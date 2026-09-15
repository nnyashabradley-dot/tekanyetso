#!/usr/bin/env python3
"""
Tekanyetso — data pipeline.

Pulls the Bank of Botswana daily exchange-rate table, stores every raw download
immutably, screens it for the fault modes BoB's table actually exhibits, and
maintains an append-only clean table recording when each observation was first
seen.

Screening policy: a bad row is QUARANTINED, never fatal. One mistyped row in
2019 must not stop today's prediction from being logged. Only structural
failures — no data, unparseable, a whole column inverted or rescaled — abort
the run, because those mean the table itself has changed shape.

Run daily:  python3 fetch.py
"""

import argparse
import csv
import hashlib
import io
import os
import statistics
import urllib.request
from datetime import date, datetime

BOB_CSV = "https://www.bankofbotswana.bw/export/exchange-rates.csv?page&_format=csv"

RAW_DIR = "raw"
CLEAN = "data/rates.csv"
QUARANTINE = "data/quarantine.csv"

# BoB quotes FOREIGN CURRENCY PER PULA.
# NOTE: the CSV export orders columns Date,CHN,EUR,GBP,USD,SDR,YEN,ZAR while the
# HTML table on the site orders them Date,CHN,EUR,GBP,USD,ZAR,SDR,YEN. Both label
# their headers correctly, so reading by NAME rather than position is
# load-bearing. Never index these columns positionally.
COLS = ["CHN", "EUR", "GBP", "USD", "ZAR", "SDR", "YEN"]

# Only these two feed the model. Everything else is carried for completeness and
# must never be able to block a row. BoB did not quote the yuan for most of the
# history and wrote 0.0000 into those cells rather than leaving them blank, so
# more than half the CHN column is zeros. A zero exchange rate is not a price;
# it is a missing value wearing a number, and it is treated as missing here.
MODEL_COLS = ["ZAR", "SDR"]
MIN_COVERAGE = 200   # non-missing observations before a column's median is judged

# The export serves the entire history in one response (~6,400 rows, 255 pages
# of 25 on the site, back to roughly 2001). A short response means the endpoint
# changed, not that the history shrank.
MIN_ROWS = 4000

# Per-row bounds. Deliberately loose: these exist to catch a single wild cell,
# not to judge borderline observations.
RANGES = {
    "ZAR": (0.5, 3.0),
    "SDR": (0.02, 0.30),
    "USD": (0.03, 0.40),
    "EUR": (0.03, 0.40),
    "GBP": (0.02, 0.30),
    "CHN": (0.2, 2.0),
    "YEN": (3.0, 40.0),
}

# Column-median bounds. This is the real inversion test. A per-row range check
# cannot catch an inverted ZAR: 1/1.23 = 0.81 sits comfortably inside [0.5, 3.0]
# because a Pula and a rand are worth roughly the same. The median of the whole
# column cannot hide like that.
MEDIAN_RANGES = {
    "ZAR": (1.0, 2.0),
    "SDR": (0.04, 0.09),
    "USD": (0.06, 0.12),
    "EUR": (0.05, 0.11),
    "GBP": (0.04, 0.10),
    "CHN": (0.4, 0.8),
    "YEN": (8.0, 14.0),
}

# Relative deviation from the local median that marks a row as suspect. The real
# fault in BoB's table is column contamination — one cell carrying another
# column's value for that day — which lands well inside RANGES but far outside
# its own neighbours. SDR is a basket and barely moves day to day, so it gets the
# tightest tolerance; ZAR is a single volatile currency and gets the loosest.
SPIKE_TOL = {"SDR": 0.03, "ZAR": 0.06, "USD": 0.05,
             "EUR": 0.05, "GBP": 0.05, "CHN": 0.06, "YEN": 0.06}
SPIKE_WINDOW = 5   # observations either side
SPIKE_MAX_GAP = 21  # ...but only those within this many CALENDAR days.
# The series has real holes — BoB skipped 26 Sep to 02 Oct 2025, and 17 to 22
# Jul 2026. Without the calendar bound, a row beside a hole gets compared
# against observations weeks away and is flagged for ordinary drift.


def fetch(url):
    req = urllib.request.Request(url, headers={"User-Agent": "tekanyetso/0.2"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.status, r.read()


def parse(blob):
    """BoB's export -> list of dicts. Tolerant of column reordering."""
    text = blob.decode("utf-8-sig", errors="replace")
    rows = list(csv.DictReader(io.StringIO(text)))
    out = []
    for row in rows:
        norm = {k.strip().upper(): (v or "").strip() for k, v in row.items() if k}
        raw_date = next((norm[k] for k in norm if "DATE" in k), None)
        if not raw_date:
            continue
        d = None
        for fmt in ("%d %b %Y", "%Y-%m-%d", "%d/%m/%Y", "%d %B %Y"):
            try:
                d = datetime.strptime(raw_date, fmt).date()
                break
            except ValueError:
                continue
        if d is None:
            continue
        rec = {"date": d.isoformat()}
        for c in COLS:
            try:
                v = float(norm.get(c, ""))
            except ValueError:
                v = None
            rec[c] = None if (v is None or v == 0.0) else v
        out.append(rec)
    return out


def check_structure(rows):
    """Fatal checks. These mean the table changed shape, not that a row is bad.

    Medians are taken over PRESENT values only. Judging a column on a median
    that includes its missing-data zeros says nothing about whether the column
    is inverted — it only says the column is sparse.
    """
    if len(rows) < MIN_ROWS:
        raise SystemExit(
            f"only {len(rows)} rows (expected >= {MIN_ROWS}); export endpoint "
            f"may have changed — raw kept, clean not touched")
    for c in MODEL_COLS:
        n = sum(1 for r in rows if r[c] is not None)
        if n < MIN_ROWS:
            raise SystemExit(
                f"column {c} has only {n} present values of {len(rows)} rows, "
                f"and the model cannot run without it — raw kept, clean not touched")
    for c, (lo, hi) in MEDIAN_RANGES.items():
        vals = [r[c] for r in rows if r[c] is not None]
        if len(vals) < MIN_COVERAGE and c not in MODEL_COLS:
            continue          # too sparse to judge, and not load-bearing
        med = statistics.median(vals)
        if not (lo <= med <= hi):
            raise SystemExit(
                f"column {c} has median {med:.4f} over {len(vals)} present "
                f"values, outside [{lo}, {hi}] — the series looks inverted or "
                f"rescaled; raw kept, clean not touched")


def screen(rows):
    """Split rows into (clean, quarantined). Never raises.

    Three faults, all present in the live table:
      1. duplicate dates carrying different values — a mistyped date elsewhere
         in the table landing on a date that already exists;
      2. column contamination — one cell holding another column's value for
         that day, e.g. 02 Nov 2023 SDR = 0.0605, which is that day's GBP;
      3. isolated typos, e.g. 05 Apr 2019 CHN = 0.0939 in a 0.63 neighbourhood.
    All three survive a per-row range check. None survives comparison against
    the row's own neighbours.
    """
    rows = sorted(rows, key=lambda r: r["date"])
    reasons = {i: [] for i in range(len(rows))}

    # 1. duplicate dates -> quarantine every copy; we cannot tell which is real
    seen = {}
    for i, r in enumerate(rows):
        seen.setdefault(r["date"], []).append(i)
    for idxs in seen.values():
        if len(idxs) > 1:
            for i in idxs:
                reasons[i].append(f"duplicate date ({len(idxs)} copies)")

    # 2. a row is unusable only if a column the MODEL needs is missing.
    #    A missing CHN is fine — nothing reads it.
    for i, r in enumerate(rows):
        for c in MODEL_COLS:
            if r[c] is None:
                reasons[i].append(f"{c} missing")

    # 3. per-row range, present values only
    for i, r in enumerate(rows):
        for c, (lo, hi) in RANGES.items():
            if r[c] is not None and not (lo <= r[c] <= hi):
                reasons[i].append(f"{c}={r[c]} outside [{lo},{hi}]")

    # 4. spike against local median, with column-identity noted as corroboration
    days = [date.fromisoformat(r["date"]) for r in rows]
    for c in COLS:
        vals = [r[c] for r in rows]
        for i, v in enumerate(vals):
            if v is None:
                continue
            lo = max(0, i - SPIKE_WINDOW)
            hi = min(len(vals), i + SPIKE_WINDOW + 1)
            nbrs = [vals[j] for j in range(lo, hi)
                    if j != i and vals[j] is not None
                    and abs((days[j] - days[i]).days) <= SPIKE_MAX_GAP]
            if len(nbrs) < 4:
                continue
            med = statistics.median(nbrs)
            if med <= 0:
                continue
            dev = abs(v - med) / med
            if dev > SPIKE_TOL[c]:
                twin = next((o for o in COLS
                             if o != c and rows[i][o] is not None
                             and rows[i][o] == v), None)
                note = f" (== {twin} that day)" if twin else ""
                reasons[i].append(
                    f"{c}={v} is {dev * 100:.1f}% off local median {med:.4f}{note}")

    clean = [r for i, r in enumerate(rows) if not reasons[i]]
    bad = [dict(r, reason="; ".join(reasons[i]))
           for i, r in enumerate(rows) if reasons[i]]
    return clean, bad


def load_clean():
    if not os.path.exists(CLEAN):
        return {}
    with open(CLEAN) as f:
        return {r["date"]: r for r in csv.DictReader(f)}


def write_quarantine(bad):
    if not bad:
        return
    fields = ["date"] + COLS + ["reason"]
    with open(QUARANTINE, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in bad:
            w.writerow({k: ("" if r[k] is None else r[k]) for k in fields})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--from-file", help="screen a local CSV instead of fetching")
    a = ap.parse_args()

    os.makedirs(RAW_DIR, exist_ok=True)
    os.makedirs(os.path.dirname(CLEAN), exist_ok=True)
    today = date.today().isoformat()

    if a.from_file:
        blob, status = open(a.from_file, "rb").read(), 0
    else:
        status, blob = fetch(BOB_CSV)
    digest = hashlib.sha256(blob).hexdigest()

    rows = parse(blob)
    if not a.from_file:
        check_structure(rows)
    clean_rows, bad = screen(rows)

    # Raw is written once per (day, content-hash) and never overwritten.
    raw_path = os.path.join(RAW_DIR, f"{today}_{digest[:12]}.csv")
    if not os.path.exists(raw_path):
        with open(raw_path, "wb") as f:
            f.write(blob)

    # Clean table is append-only. first_seen is never rewritten; a value that
    # changes under us is reported, not silently overwritten.
    existing = load_clean()
    added, revised = 0, []
    for r in clean_rows:
        old = existing.get(r["date"])
        if old is None:
            r["first_seen"] = today
            existing[r["date"]] = r
            added += 1
        else:
            for c in COLS:
                o = old[c] if old[c] not in ("", None) else None
                n = r[c]
                if o is None and n is None:
                    continue
                if o is None or n is None or abs(float(o) - n) > 1e-12:
                    revised.append(f"{r['date']} {c}: {old[c]!r} -> {n!r}")

    fields = ["date"] + COLS + ["first_seen"]
    with open(CLEAN, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for d in sorted(existing):
            w.writerow({k: ("" if existing[d][k] is None else existing[d][k])
                        for k in fields})

    write_quarantine(bad)

    with open("fetch.log", "a") as f:
        f.write(f"{datetime.now().isoformat()}\tHTTP {status}\t{digest[:16]}\t"
                f"{len(rows)} rows\t{len(bad)} quarantined\t+{added} new\t"
                f"{len(revised)} revisions\n")

    print(f"HTTP {status}  sha256 {digest[:16]}  {len(rows)} parsed  "
          f"{len(clean_rows)} clean  {len(bad)} quarantined  +{added} new")
    if bad:
        print(f"\nquarantined (written to {QUARANTINE}):")
        for r in bad[:15]:
            print(f"   {r['date']}  {r['reason']}")
        if len(bad) > 15:
            print(f"   ... and {len(bad) - 15} more")
    if revised:
        print(f"\n!! {len(revised)} REVISIONS to already-seen values — investigate:")
        for line in revised[:10]:
            print("   ", line)
    if existing:
        print(f"\ncoverage: {min(existing)} .. {max(existing)}  ({len(existing)} days)")
        print("\nper-column coverage (blank cells are BoB's, not ours):")
        for c in COLS:
            present = [d for d in sorted(existing)
                       if existing[d][c] not in ("", None)]
            mark = "  <-- used by the model" if c in MODEL_COLS else ""
            if present:
                print(f"  {c:<5}{len(present):>6} values   from {present[0]}{mark}")
            else:
                print(f"  {c:<5}{0:>6} values{mark}")


if __name__ == "__main__":
    main()
