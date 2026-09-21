#!/usr/bin/env python3
"""
Tekanyetso — estimator.

The peg holds a weighted log basket of the Pula's value on a linear path:

    w*log(ZAR_per_Pula) + (1-w)*log(SDR_per_Pula) = a + b*t

Writing u = log(ZAR_pP) - log(SDR_pP), this rearranges to a plain OLS with
one regressor of interest:

    log(SDR_pP) = a + b*t - w*u + e

so w = -coef(u) and the annual rate of crawl is 365*b. The sum-to-one
constraint is imposed by construction, not estimated, which removes the
degree of freedom that would otherwise make w and the crawl collinear.

t is in CALENDAR days, which recovers the annual rate correctly. The Bank in
fact steps the parity once per published day, not per calendar day (FINDINGS
section 3), so this model over-predicts the crawl across weekends; backtest.py
measures the cost. Annualising a per-publication slope by 365 rather than by
publications per year would put the rate out by about 45 percent (365/252).

Usage:  python3 estimate.py [--window 90] [--log]
"""

import argparse
import csv
import os
from datetime import date, datetime, timedelta

import numpy as np

CLEAN = "data/rates.csv"
PRED_LOG = "predictions.csv"


def load():
    with open(CLEAN) as f:
        rows = sorted(csv.DictReader(f), key=lambda r: r["date"])
    # Blank cells are BoB's missing data, not ours. Only ZAR and SDR matter; a
    # row missing either cannot be used, a row missing anything else is fine.
    usable = [r for r in rows if r.get("ZAR") and r.get("SDR")]
    dropped = len(rows) - len(usable)
    if dropped:
        print(f"  ({dropped} rows skipped: ZAR or SDR missing)")
    rows = usable
    if not rows:
        raise SystemExit("no rows with both ZAR and SDR — check data/rates.csv")
    d = np.array([date.fromisoformat(r["date"]) for r in rows])
    t = np.array([(x - d[0]).days for x in d], dtype=float)
    lz = np.log(np.array([float(r["ZAR"]) for r in rows]))
    ls = np.log(np.array([float(r["SDR"]) for r in rows]))
    return d, t, lz, ls


def fit(t, lz, ls):
    """Return weight on rand, annual crawl %, residual RMSE in bp, condition number."""
    u = lz - ls
    X = np.column_stack([np.ones_like(t), t, u])
    beta, *_ = np.linalg.lstsq(X, ls, rcond=None)
    resid = ls - X @ beta

    # Standard errors. Without these the headline "recovered -2.76 against an
    # announced -2.76" is an anecdote; with them it is a measurement. Classical
    # OLS covariance: sigma^2 (X'X)^-1, sigma^2 = RSS/(n-k).
    n, k = X.shape
    se_w = se_crawl = float("nan")
    if n > k:
        s2 = float(resid @ resid) / (n - k)
        try:
            cov = s2 * np.linalg.inv(X.T @ X)
            se_w = float(np.sqrt(cov[2, 2]))
            se_crawl = float(np.sqrt(cov[1, 1])) * 365 * 100
        except np.linalg.LinAlgError:
            pass

    # Identification diagnostic (build notes §5). The crawl and the weight are
    # separable only if the rand/SDR spread u moves independently of time
    # within the window. Condition the two non-constant regressors after
    # demeaning and unit-scaling; including the intercept would make this
    # singular by construction and report inf.
    Z = np.column_stack([t - t.mean(), u - u.mean()])
    scale = np.linalg.norm(Z, axis=0)
    cond = float(np.linalg.cond(Z / np.where(scale > 0, scale, 1.0)))

    return {
        "w": -beta[2],
        "se_w": se_w,
        "crawl_pct": 365 * beta[1] * 100,
        "se_crawl": se_crawl,
        "rmse_bp": 1e4 * float(np.sqrt((resid**2).mean())),
        "cond": cond,
        "n": len(t),
    }


def rolling(d, t, lz, ls, window):
    out = []
    for i in range(window, len(t) + 1):
        s = slice(i - window, i)
        f = fit(t[s] - t[s][0], lz[s], ls[s])
        f["date"] = d[s][-1]
        out.append(f)
    return out


def quantisation_floor(rates_path=CLEAN):
    """BoB rounds every column to 4dp. On SDR (~0.056) that is coarse.
    Returns the resulting noise floor on the basket index, in bp."""
    with open(rates_path) as f:
        rows = [r for r in csv.DictReader(f) if r.get("ZAR") and r.get("SDR")]
    last = rows[-1]
    sdr, zar = float(last["SDR"]), float(last["ZAR"])
    # uniform rounding to 1e-4 -> sd = 1e-4/sqrt(12), relative, weighted 0.5 each
    sd = lambda v: (1e-4 / np.sqrt(12)) / v
    return 1e4 * 0.5 * np.hypot(sd(sdr), sd(zar))


def log_prediction(pred_index, for_date, params, latest_obs):
    """Append-only. A row is written once and never rewritten.

    latest_obs is the last date present in rates.csv when this ran. Recording it
    is what makes the row self-describing: the same code predicts today's rate
    if BoB has not published yet, and tomorrow's if it has, and those are
    different claims. lead_days is the gap; bob_published_today says which case
    this was.
    """
    seen = set()
    if os.path.exists(PRED_LOG):
        with open(PRED_LOG) as f:
            seen = {r["for_date"] for r in csv.DictReader(f)}
    if for_date.isoformat() in seen:
        print(f"prediction for {for_date} already logged - not overwriting")
        return

    now = datetime.now()
    lead = (for_date - latest_obs).days
    new = not os.path.exists(PRED_LOG)
    with open(PRED_LOG, "a", newline="") as f:
        w = csv.writer(f)
        if new:
            w.writerow(["made_at", "for_date", "pred_log_index", "w",
                        "crawl_pct", "window", "model",
                        "latest_obs", "lead_days", "bob_published_today"])
        w.writerow([now.isoformat(timespec="seconds"), for_date.isoformat(),
                    f"{pred_index:.8f}", f"{params['w']:.4f}",
                    f"{params['crawl_pct']:.3f}", params["n"], "rolling-ols-v0",
                    latest_obs.isoformat(), lead,
                    "yes" if latest_obs == now.date() else "no"])
    print(f"logged prediction for {for_date} "
          f"({lead}d after last observation {latest_obs})")
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--window", type=int, default=90)
    ap.add_argument("--log", action="store_true", help="append tomorrow's prediction")
    a = ap.parse_args()

    d, t, lz, ls = load()
    print(f"{len(t)} observations, {d[0]} .. {d[-1]}")

    floor = quantisation_floor()
    print(f"quantisation floor from BoB's 4dp rounding: {floor:.2f} bp")
    print("  (a residual at or below this means the peg holds to measurement precision)\n")

    full = fit(t, lz, ls)
    print(f"full sample : w={full['w']:.3f}+/-{full['se_w']:.3f}  "
          f"crawl={full['crawl_pct']:+.2f}+/-{full['se_crawl']:.2f} %/yr  "
          f"resid={full['rmse_bp']:.2f} bp")
    print("  (meaningless if the sample spans changepoints — it is a smoke test only)\n")

    r = rolling(d, t, lz, ls, a.window)
    print(f"rolling {a.window}d, last 10 windows:")
    print(f"  {'date':<12}{'w':>8}{'se':>8}{'crawl %/yr':>13}{'se':>8}"
          f"{'resid bp':>11}{'cond':>7}")
    for f in r[-10:]:
        flag = "  <-- ill-conditioned" if f["cond"] > 10 else ""
        print(f"  {str(f['date']):<12}{f['w']:8.3f}{f['se_w']:8.3f}"
              f"{f['crawl_pct']:13.2f}{f['se_crawl']:8.2f}"
              f"{f['rmse_bp']:11.2f}{f['cond']:7.0f}{flag}")

    with open("rolling.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["date", "w", "se_w", "crawl_pct", "se_crawl", "rmse_bp", "cond"])
        for f in r:
            w.writerow([f["date"], f"{f['w']:.5f}", f"{f['se_w']:.5f}",
                        f"{f['crawl_pct']:.4f}", f"{f['se_crawl']:.4f}",
                        f"{f['rmse_bp']:.3f}", f"{f['cond']:.1f}"])
    print("\nwrote rolling.csv")

    # --- prediction ---------------------------------------------------------
    # The falsifiable claim is about the BANK, not about the FX market: the log
    # basket index tomorrow equals today's plus one day of crawl. This needs no
    # forecast of the rand or the SDR, so it is testable from the BoB file alone.
    cur = r[-1]
    idx_today = cur["w"] * lz[-1] + (1 - cur["w"]) * ls[-1]
    nxt = d[-1] + timedelta(days=1)
    while nxt.weekday() >= 5:
        nxt += timedelta(days=1)
    step = (nxt - d[-1]).days
    pred = idx_today + (cur["crawl_pct"] / 100 / 365) * step
    print(f"\nprediction for {nxt} ({step}d ahead)")
    print(f"  log basket index: {pred:.6f}   (today {idx_today:.6f})")
    print(f"  using w={cur['w']:.3f}, crawl={cur['crawl_pct']:+.2f} %/yr")

    if a.log:
        log_prediction(pred, nxt, cur, d[-1])
    else:
        print("  (not logged; pass --log to append to predictions.csv)")


if __name__ == "__main__":
    main()
