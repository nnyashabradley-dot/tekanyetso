#!/usr/bin/env python3
"""
Tekanyetso -- live scoreboard.

Scores every row of predictions.csv against what the Bank of Botswana then
published, and writes the result into README.md between two marker lines, so
anyone landing on the repository can read the track record without working it
out from a CSV.

Each prediction is scored the way it was made: the log basket index
w*log(ZAR) + (1-w)*log(SDR) on the target date, using the weight logged with
the prediction.

The coverage line is the important one. Timestamps stop a prediction being
backdated; nothing stops someone running the system daily and only keeping the
days that went well. Publishing predictions-logged against publication-days-
elapsed makes a missing day visible.

    python3 scoreboard.py            # print, and update README.md
    python3 scoreboard.py --dry-run  # print only
"""

import argparse
import csv
from datetime import date, datetime

import numpy as np

RATES = "data/rates.csv"
PREDS = "predictions.csv"
README = "README.md"
START = "<!-- SCOREBOARD:START -->"
END = "<!-- SCOREBOARD:END -->"

# Out-of-sample error of this exact model over 2011-2026, from backtest.py.
# Quoted so a reader can judge the live errors against the long-run record.
BACKTEST_RMSE_BP = 3.56


def load_rates(path=RATES):
    out = {}
    with open(path) as f:
        for r in csv.DictReader(f):
            out[date.fromisoformat(r["date"])] = (float(r["ZAR"]), float(r["SDR"]))
    return out


def score(preds, rates):
    latest = max(rates)
    rows = []
    for p in preds:
        fd = date.fromisoformat(p["for_date"])
        made = datetime.fromisoformat(p["made_at"])
        w = float(p["w"])
        row = {"for_date": fd, "made_at": made, "pred": float(p["pred_log_index"]),
               "w": w, "same_day": made.date() == fd}
        if fd in rates:
            zar, sdr = rates[fd]
            actual = w * np.log(zar) + (1 - w) * np.log(sdr)
            row["err_bp"] = 1e4 * (actual - row["pred"])
            row["status"] = "scored"
        elif fd > latest:
            row["status"] = "pending"
        else:
            row["status"] = "no publication"       # a holiday
        rows.append(row)
    return rows, latest


def coverage(rows, rates):
    """Publication days from the first prediction's target to the latest
    published date, against how many of them had a prediction."""
    if not rows:
        return 0, 0, []
    first = min(r["for_date"] for r in rows)
    latest = max(rates)
    pub = sorted(x for x in rates if first <= x <= latest)
    have = {r["for_date"] for r in rows}
    missed = [x for x in pub if x not in have]
    return len(pub) - len(missed), len(pub), missed


def render(rows, latest, cov):
    got, total, missed = cov
    scored = [r for r in rows if r["status"] == "scored"]
    e = np.array([r["err_bp"] for r in scored]) if scored else np.array([])
    lines = [START, "", "## Live track record", "",
             f"*Updated {datetime.now():%Y-%m-%d %H:%M}. Latest Bank publication: "
             f"{latest}.*", ""]
    if len(e):
        lines += [f"- **{len(e)} predictions scored.** Mean absolute error "
                  f"{np.abs(e).mean():.2f} bp, RMS {np.sqrt((e**2).mean()):.2f} bp, "
                  f"mean {e.mean():+.2f} bp.",
                  f"- For scale: the same model's out-of-sample RMS over 2011-2026 "
                  f"in the backtest is {BACKTEST_RMSE_BP:.2f} bp, and the Bank's own "
                  f"four-decimal rounding puts a floor of about 3 bp under any "
                  f"prediction anchored on a published rate."]
    lines += [f"- **Coverage: {got} of {total} publication days predicted.**"
              + (f" Missed: {', '.join(str(x) for x in missed)}." if missed else ""),
              "- Scheduled predictions run at 09:00 Gaborone time; the Made column "
              "shows the actual time. A *same day* prediction was made when the "
              "Bank's file did not yet contain that day's rate.", "",
              "| Target date | Made | Predicted | Error (bp) |",
              "|---|---|---:|---:|"]
    for r in sorted(rows, key=lambda r: r["for_date"], reverse=True)[:10]:
        made = r["made_at"].strftime("%Y-%m-%d %H:%M")
        tag = " (same day)" if r["same_day"] else ""
        err = f"{r['err_bp']:+.2f}" if r["status"] == "scored" else r["status"]
        lines.append(f"| {r['for_date']} | {made}{tag} | {r['pred']:.6f} | {err} |")
    lines += ["", "Every prediction ever made is in `predictions.csv`, which is "
              "append-only; its git history shows when each row was written.",
              "", END]
    return "\n".join(lines)


def update_readme(block, path=README):
    try:
        text = open(path, encoding="utf-8").read()
    except FileNotFoundError:
        text = "# Tekanyetso\n"
    if START in text and END in text:
        a, rest = text.split(START, 1)
        _, b = rest.split(END, 1)
        new = a + block + b
    else:
        # after the title and first paragraph
        parts = text.split("\n\n", 2)
        head = "\n\n".join(parts[:2]) if len(parts) >= 2 else text
        tail = parts[2] if len(parts) == 3 else ""
        new = head + "\n\n" + block + "\n\n" + tail
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(new)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    rates = load_rates()
    preds = list(csv.DictReader(open(PREDS)))
    rows, latest = score(preds, rates)
    block = render(rows, latest, coverage(rows, rates))
    print(block)
    if not a.dry_run:
        update_readme(block)
        print(f"\nupdated {README}")


if __name__ == "__main__":
    main()
