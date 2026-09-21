#!/usr/bin/env python3
"""
Tekanyetso -- no-lookahead backtest.

Replays the history day by day. At each published day i, using only data
published on or before d[i], each model predicts the log basket index for the
next weekday -- exactly the claim the live log makes. The prediction is scored
against what the Bank then published for that date, using the weight the model
itself used, which is how predictions.csv is scored too.

Scoring against the model's own weight is the same as predicting the SDR rate
given the rand/SDR spread on the target day: the two errors are identical,
because I = ls + w*u and the w*u term cancels. So every model is scored on one
common target, the Bank's published SDR-per-Pula, and no model forecasts the
currency market.

The models
----------
  live         rolling 90-observation OLS, calendar-day crawl, anchored on
               today's published index. This is rolling-ols-v0, the model that
               writes predictions.csv, reimplemented from estimate.fit.
  trading      same fit, but the crawl is one step per published day
               (FINDINGS section 3), still anchored on today's published index.
  trading-fit  same, anchored on the regression's fitted index instead of the
               published one, which carries today's rounding noise.
  kalman       kalman.py, trading steps, process noise fitted on
               2005-06 .. 2010-12 only and then frozen.
  no-crawl     today's index, unchanged. The baseline.

Proof of blindness
------------------
Structure is not proof, so the backtest checks itself. It re-runs every model
on the history truncated at several cutoffs, and again with every value after
the cutoff replaced by random numbers, and requires every prediction made on or
before the cutoff to be bit-identical to the full run. If any model had peeked,
this fails.

The Kalman hyperparameters are the one place future data could leak: fitted on
2005-2010, they would be lookahead for predictions made inside 2005-2010. The
headline table therefore starts on 2011-01-01, where they are genuinely
out of sample.

    python3 backtest.py                  # run, check blindness, write backtest.csv
    python3 backtest.py --skip-poison    # faster, no blindness check
    python3 backtest.py --check-live     # compare against predictions.csv
"""

import argparse
import csv
from datetime import date, timedelta

import numpy as np

from estimate import fit as live_fit
from kalman import FIT_TO, INIT_N, fit_q, init_state, kfilter, steps
from segment import load, noise_sd

START = date(2005, 6, 1)
TEST_FROM = date(2011, 1, 1)
W = 90
MODELS = ["live", "trading", "trading-fit", "kalman", "no-crawl"]


def next_weekday(x):
    y = x + timedelta(days=1)
    while y.weekday() >= 5:
        y += timedelta(days=1)
    return y


def load_series(path, end=None):
    d, t, u, ls, sg, sdr, zar = load(path, START, end)
    lz = ls + u
    return {"d": d, "t": t, "u": u, "ls": ls, "lz": lz, "sdr": sdr, "zar": zar}


def predict_all(S, q=None):
    """Every model's prediction made at each index i, for next_weekday(d[i]).

    Returns dict model -> (pred_index[i], w_used[i]); NaN where a model has no
    prediction yet. Uses only S's arrays; nothing here looks at an index > i
    except the Kalman filter's own sequential pass, which the poison test
    covers."""
    d, ls, lz, u = S["d"], S["ls"], S["lz"], S["u"]
    n = len(d)
    out = {m: (np.full(n, np.nan), np.full(n, np.nan)) for m in MODELS}
    tdays = np.array([(x - d[0]).days for x in d], dtype=float)
    k = np.arange(W, dtype=float)

    for i in range(W - 1, n):
        s = slice(i - W + 1, i + 1)
        step_days = (next_weekday(d[i]) - d[i]).days

        f = live_fit(tdays[s] - tdays[s][0], lz[s], ls[s])
        w = f["w"]
        idx = w * lz[i] + (1 - w) * ls[i]
        out["live"][0][i] = idx + f["crawl_pct"] / 100 / 365 * step_days
        out["live"][1][i] = w
        out["no-crawl"][0][i] = idx
        out["no-crawl"][1][i] = w

        X = np.column_stack([np.ones(W), k, -u[s]])
        a, b, w2 = np.linalg.lstsq(X, ls[s], rcond=None)[0]
        idx2 = w2 * lz[i] + (1 - w2) * ls[i]
        out["trading"][0][i] = idx2 + b
        out["trading"][1][i] = w2
        out["trading-fit"][0][i] = a + b * W       # fitted index one step on
        out["trading-fit"][1][i] = w2

    if q is not None:
        dt = steps(d, S["t"], "trading")
        x0, P0 = init_state(ls, u, dt, S["sdr"], S["zar"])
        # Run the filter, then turn its state after i into a prediction of the
        # index for i+1. The index uses the predicted weight, and the index
        # prediction does not need u[i+1] (see the next comment).
        pred_ls, pred_var, xf, _ = kfilter(ls, u, dt, S["sdr"], S["zar"], q,
                                           x0, P0, INIT_N)
        for i in range(INIT_N, n):
            m_, b_, w_ = xf[i]
            # m_{i+1} = m_i + b - w_i (u_i - u_{i-1}); index at i+1 = m_{i+1} + w u_i
            m_next = m_ + b_ - w_ * (u[i] - u[i - 1])
            out["kalman"][0][i] = m_next + w_ * u[i]
            out["kalman"][1][i] = w_
    return out


def score(S, preds):
    """Error in bp for each model at each prediction index, NaN if unscorable
    (the target weekday was a holiday with no publication)."""
    d, ls, lz = S["d"], S["ls"], S["lz"]
    pos = {x: j for j, x in enumerate(d)}
    n = len(d)
    target = np.full(n, -1)
    for i in range(n):
        target[i] = pos.get(next_weekday(d[i]), -1)
    errs = {}
    for m, (p, w) in preds.items():
        e = np.full(n, np.nan)
        ok = (target >= 0) & ~np.isnan(p)
        j = target[ok]
        e[ok] = 1e4 * ((w[ok] * lz[j] + (1 - w[ok]) * ls[j]) - p[ok])
        errs[m] = e
    return errs, target


def summary(S, errs, target, lo, hi, label):
    d = S["d"]
    sel = np.array([lo <= x <= hi for x in d])
    j = target[sel & (target >= 0)]
    fl = noise_sd(S["sdr"][j], S["zar"][j], 0.5) * 1e4
    print(f"\n{label}: predictions made {lo} .. {hi}")
    print(f"  rounding floor on the target day: mean {fl.mean():.2f} bp; an anchor on"
          f"\n  today's published value adds today's rounding too, so ~{np.sqrt(2)*fl.mean():.2f} bp")
    print(f"  {'model':<13}{'n':>6}{'bias':>8}{'RMSE':>8}{'MAE':>8}{'p95':>8}"
          f"{'Mon bias':>10}{'Tue-Fri':>9}")
    dow = np.array([next_weekday(x).weekday() for x in d])
    rows = {}
    for m in MODELS:
        e = errs[m][sel]
        ok = ~np.isnan(e)
        if ok.sum() == 0:
            continue
        v = e[ok]
        mon = e[ok & (dow[sel] == 0)]
        rest = e[ok & (dow[sel] != 0)]
        rows[m] = dict(n=int(ok.sum()), bias=v.mean(), rmse=np.sqrt((v**2).mean()),
                       mae=np.abs(v).mean(), p95=np.percentile(np.abs(v), 95),
                       mon=mon.mean(), rest=rest.mean())
        r = rows[m]
        print(f"  {m:<13}{r['n']:>6}{r['bias']:>8.2f}{r['rmse']:>8.2f}{r['mae']:>8.2f}"
              f"{r['p95']:>8.2f}{r['mon']:>10.2f}{r['rest']:>9.2f}")
    return rows


def dm_test(e1, e2, lag=5):
    """Diebold-Mariano test on squared errors, Newey-West variance.
    Positive t: model 2 is more accurate. Errors from models anchored on
    today's published value are autocorrelated (today's rounding appears in
    two consecutive errors), which is why the long-run variance is used."""
    ok = ~np.isnan(e1) & ~np.isnan(e2)
    dd = e1[ok] ** 2 - e2[ok] ** 2
    n = len(dd)
    x = dd - dd.mean()
    lrv = (x @ x) / n
    for k in range(1, lag + 1):
        lrv += 2 * (1 - k / (lag + 1)) * (x[k:] @ x[:-k]) / n
    return dd.mean() / np.sqrt(lrv / n), n


def after_breaks(S, errs, lo, horizon=20):
    """Error just after a policy change versus everywhere else. Break dates
    come from PELT on the FULL history -- used only to sort the errors for
    this table, never by any model."""
    from segment import P, SegmentCost, segment_pelt
    d = S["d"]
    sg = noise_sd(S["sdr"], S["zar"], 0.5)
    c = SegmentCost(S["t"], S["u"], S["ls"], sg, S["sdr"], S["zar"])
    cuts, _ = segment_pelt(c, 10 * P * np.log(c.n), 90)
    near = np.zeros(len(d), bool)
    for k in cuts:
        near[max(0, k - 1):k + horizon] = True
    sel = np.array([x >= lo for x in d])
    print(f"\nerrors in the {horizon} publications after each of the {len(cuts)} detected"
          f"\npolicy changes, against all other days (predictions made from {lo})")
    print(f"  {'model':<13}{'RMSE near':>11}{'RMSE other':>12}{'worst near':>12}")
    for m in MODELS:
        e = errs[m]
        a = e[sel & near & ~np.isnan(e)]
        b = e[sel & ~near & ~np.isnan(e)]
        print(f"  {m:<13}{np.sqrt((a**2).mean()):>11.2f}{np.sqrt((b**2).mean()):>12.2f}"
              f"{np.abs(a).max():>12.1f}")


def poison_check(path, full_preds, S_full, q, cutoffs):
    """Truncate and corrupt at each cutoff; require bit-identical predictions."""
    print("\nblindness check: predictions made on or before each cutoff must not"
          "\nchange when the data after it is deleted, or replaced with noise")
    rng = np.random.default_rng(0)
    allok = True
    for c in cutoffs:
        S_tr = load_series(path, end=c)
        k = len(S_tr["d"]) - 1
        P_tr = predict_all(S_tr, q)

        S_bad = {key: (v.copy() if isinstance(v, np.ndarray) else v)
                 for key, v in S_full.items()}
        nb = len(S_bad["d"]) - k - 1
        S_bad["ls"][k + 1:] = rng.normal(-3, 1, nb)
        S_bad["lz"][k + 1:] = rng.normal(0, 1, nb)
        S_bad["u"][k + 1:] = S_bad["lz"][k + 1:] - S_bad["ls"][k + 1:]
        S_bad["sdr"][k + 1:] = np.exp(S_bad["ls"][k + 1:])
        S_bad["zar"][k + 1:] = np.exp(S_bad["lz"][k + 1:])
        P_bad = predict_all(S_bad, q)

        line = []
        for m in MODELS:
            a = full_preds[m][0][:k + 1]
            same = (np.array_equal(a, P_tr[m][0][:k + 1], equal_nan=True) and
                    np.array_equal(a, P_bad[m][0][:k + 1], equal_nan=True))
            allok &= same
            line.append(f"{m} {'ok' if same else 'CHANGED'}")
        print(f"  cutoff {c}: " + ", ".join(line))
    print("  PASS: no model uses data from after the day it predicts from"
          if allok else "  FAIL: a model is looking ahead")
    return allok


def check_live(S, preds, path_log="predictions.csv"):
    """Does the backtest reproduce what the live system actually logged?"""
    try:
        rows = list(csv.DictReader(open(path_log)))
    except FileNotFoundError:
        print("predictions.csv not found")
        return
    d = S["d"]
    print("\nlive log against the backtest's reconstruction of rolling-ols-v0")
    print(f"  {'for_date':<12}{'logged':>14}{'backtest':>14}{'diff':>10}")
    for r in rows:
        fd = date.fromisoformat(r["for_date"])
        i = int(np.searchsorted(d, fd)) - 1       # last publication before target
        if i < 0 or next_weekday(d[i]) != fd:
            print(f"  {r['for_date']:<12}  cannot reconstruct")
            continue
        b = preds["live"][0][i]
        print(f"  {r['for_date']:<12}{float(r['pred_log_index']):>14.8f}{b:>14.8f}"
              f"{abs(float(r['pred_log_index']) - b):>10.1e}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/rates.csv")
    ap.add_argument("--skip-poison", action="store_true")
    ap.add_argument("--check-live", action="store_true")
    a = ap.parse_args()

    S = load_series(a.data)
    d = S["d"]
    nfit = int(np.searchsorted(d, FIT_TO, side="right"))
    dt = steps(d, S["t"], "trading")
    print(f"{len(d)} observations {d[0]} .. {d[-1]}")
    print(f"Kalman process noise fitted on {d[0]} .. {d[nfit-1]} only")
    q, _ = fit_q(S["ls"][:nfit], S["u"][:nfit], dt[:nfit], S["sdr"][:nfit],
                 S["zar"][:nfit])

    preds = predict_all(S, q)
    errs, target = score(S, preds)

    rows = summary(S, errs, target, TEST_FROM, d[-1],
                   "HEADLINE, out of sample for every model")
    summary(S, errs, target, date(2005, 10, 1), date(2010, 12, 31),
            "early period (Kalman hyperparameters in sample here)")

    sel = np.array([x >= TEST_FROM for x in d])
    print("\nDiebold-Mariano, 2011 onward (t > 2: second model significantly better)")
    for m1, m2 in [("live", "trading"), ("live", "kalman"), ("trading", "kalman"),
                   ("no-crawl", "live")]:
        tstat, n = dm_test(np.where(sel, errs[m1], np.nan), np.where(sel, errs[m2], np.nan))
        print(f"  {m1:>9} vs {m2:<9} t = {tstat:+6.2f}   (n = {n})")

    after_breaks(S, errs, TEST_FROM)

    if a.check_live:
        check_live(S, preds)

    ok = True
    if not a.skip_poison:
        cuts = [date(2011, 3, 15), date(2017, 8, 1), date(2024, 12, 31)]
        ok = poison_check(a.data, preds, S, q, cuts)

    with open("backtest.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["made_on", "for_date"] + [f"err_{m}_bp" for m in MODELS])
        for i in range(len(d)):
            if target[i] < 0:
                continue
            vals = [errs[m][i] for m in MODELS]
            if all(np.isnan(v) for v in vals):
                continue
            w.writerow([d[i], d[target[i]]] +
                       ["" if np.isnan(v) else f"{v:.3f}" for v in vals])
    print("\nwrote backtest.csv")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
