#!/usr/bin/env python3
"""
Tekanyetso — residual diagnostic.

`segment.py` established that the model

    w*log(ZAR_pP) + (1-w)*log(SDR_pP) = a + b*t

fits to measurement precision from 2009 onward (chi2/dof ~ 1.0) and fails badly
before that (chi2/dof ~ 26 in 2005-06). This asks WHY, by fitting one line
across an interval with no announced policy change and examining what is left.

If the model is right, the residual is rounding noise: no autocorrelation, no
weekday pattern, RMS equal to the derived floor. Any structure is the Bank
doing something the straight line does not describe.

Leading hypothesis. The model assumes the crawl is applied continuously, every
calendar day. If instead the parity is reset every T days, the residual is a
sawtooth of period T and RMS (c/100)/365*T/sqrt(12). At the 2005-06 crawl of
-4.80 %/yr that is 0.38 bp daily, 2.7 weekly, 5.3 fortnightly, 11.4 monthly.
The observed excess of ~4 bp therefore rules out monthly on magnitude alone and
points at weekly or fortnightly. A related variant: the crawl might be applied
on TRADING days only, which would show up as a weekday pattern rather than a
pure sawtooth.

Usage
-----
    python3 diagnose.py --from 2005-06-01 --to 2006-06-29
    python3 diagnose.py --eras            # every announced-crawl interval
    python3 diagnose.py --self-test       # inject a known step, recover it
"""

import argparse
import sys
from datetime import date

import numpy as np

from segment import load, noise_sd, make_synthetic, P


# ---------------------------------------------------------------------------

def fit_window(t, u, ls, sdr, zar):
    """One weighted least squares fit. Returns residual in bp plus summary."""
    sg = noise_sd(sdr, zar, 0.5)
    X = np.column_stack([np.ones_like(t), t, u])
    b, *_ = np.linalg.lstsq(X / sg[:, None], ls / sg, rcond=None)
    w_hat = float(np.clip(-b[2], 0.0, 1.0))
    sg = noise_sd(sdr, zar, w_hat)          # refit the floor to the fitted w
    b, *_ = np.linalg.lstsq(X / sg[:, None], ls / sg, rcond=None)
    resid = ls - X @ b
    dof = max(len(t) - P, 1)
    return {
        "resid_bp": 1e4 * resid,
        "w": float(-b[2]),
        "crawl": float(b[1] * 100),
        "rms_bp": float(1e4 * np.sqrt((resid ** 2).mean())),
        "floor_bp": float(1e4 * sg.mean()),
        "chi2_dof": float(((resid / sg) ** 2).sum() / dof),
    }


def excess(rms, floor):
    """Residual not explained by rounding, removed in quadrature."""
    return float(np.sqrt(max(rms ** 2 - floor ** 2, 0.0)))


def autocorr(x, maxlag=21):
    x = x - x.mean()
    denom = (x * x).sum()
    return [1.0 if k == 0 else float((x[:-k] * x[k:]).sum() / denom)
            for k in range(maxlag)]


def sawtooth_rms(crawl_pct, T):
    """Expected residual RMS in bp if the parity steps every T calendar days."""
    return 1e4 * (abs(crawl_pct) / 100) / 365 * T / np.sqrt(12)



def period_from_acf(ac, n, maxlag=19):
    """Smallest lag that is a local maximum and significant.

    A sawtooth correlates with itself at EVERY multiple of its period, so the
    largest peak is not necessarily the period -- with noise, a harmonic often
    edges it out. The fundamental is the smallest significant peak.
    """
    thresh = 2 / np.sqrt(n)
    peaks = [k for k in range(2, maxlag)
             if ac[k] > thresh and ac[k] >= ac[k - 1] and ac[k] >= ac[k + 1]]
    return (peaks[0] if peaks else None), peaks, thresh


def fold(resid, key, labels, title):
    """Mean residual grouped by some periodic key, with standard errors."""
    print(f"\n  {title}")
    print(f"    {'group':<14}{'n':>6}{'mean bp':>10}{'se':>8}")
    flagged = False
    for k, lab in labels:
        m = key == k
        if m.sum() < 3:
            continue
        v = resid[m]
        se = v.std(ddof=1) / np.sqrt(len(v))
        star = "  *" if abs(v.mean()) > 2 * se else ""
        flagged |= bool(star)
        print(f"    {lab:<14}{len(v):>6}{v.mean():>10.2f}{se:>8.2f}{star}")
    print("    (* = mean differs from zero by more than 2 standard errors)")
    return flagged


def plot(resid, dates, height=15, width=78):
    """ASCII scatter. numpy-only project, so no matplotlib."""
    n = len(resid)
    lo, hi = resid.min(), resid.max()
    if hi - lo < 1e-9:
        return
    grid = [[" "] * width for _ in range(height)]
    for i, v in enumerate(resid):
        c = int(i * (width - 1) / max(n - 1, 1))
        r = int((hi - v) * (height - 1) / (hi - lo))
        grid[r][c] = "#" if grid[r][c] != " " else "."
    zero = int(hi * (height - 1) / (hi - lo))
    print()
    for r, row in enumerate(grid):
        axis = f"{hi - r*(hi-lo)/(height-1):>8.1f} |"
        line = "".join(row)
        if r == zero:
            line = "".join(ch if ch != " " else "-" for ch in line)
        print(f"  {axis}{line}")
    print(f"  {'':>8} +{'-'*width}")
    print(f"  {'':>8}  {dates[0]}{' '*(width-22)}{dates[-1]}   (bp)")


# ---------------------------------------------------------------------------

def diagnose(d, t, u, ls, sdr, zar, label="", quiet=False):
    f = fit_window(t, u, ls, sdr, zar)
    r = f["resid_bp"]
    exc = excess(f["rms_bp"], f["floor_bp"])

    print(f"\n{'='*84}")
    print(f"{label}   n={len(t)}   {d[0]} .. {d[-1]}")
    print(f"{'='*84}")
    print(f"  w = {f['w']:.3f}   crawl = {f['crawl']:+.2f} %/yr")
    print(f"  residual RMS {f['rms_bp']:.2f} bp   rounding floor {f['floor_bp']:.2f} bp"
          f"   chi2/dof {f['chi2_dof']:.2f}")
    print(f"  UNEXPLAINED EXCESS: {exc:.2f} bp")

    if quiet:
        return exc, f

    # -- what step period would produce this much excess? -------------------
    print("\n  If the parity steps every T days instead of crawling continuously,"
          "\n  the residual is a sawtooth of this size:")
    print(f"    {'period':<14}{'expected bp':>13}{'':>4}")
    best, bestd = None, 1e9
    for T, name in [(1, "daily"), (5, "weekly (5 obs)"), (7, "weekly"),
                    (10, "fortnightly*"), (14, "fortnightly"), (30, "monthly"),
                    (91, "quarterly")]:
        e = sawtooth_rms(f["crawl"], T)
        mark = ""
        if abs(e - exc) < bestd:
            bestd, best = abs(e - exc), name
        if abs(e - exc) / max(exc, 1e-9) < 0.35:
            mark = "  <-- consistent"
        print(f"    {name:<14}{e:>13.2f}{mark}")
    print(f"  closest match to the observed {exc:.2f} bp: {best}")

    # -- is there actually periodic structure? ------------------------------
    ac = autocorr(r)
    print("\n  autocorrelation of the residual (rounding noise is ~0 at every lag)")
    if exc < 1.0:
        print("  NOTE: the excess is under 1 bp, i.e. the residual is essentially")
        print("  at the measurement floor. Any peaks below are noise -- on")
        print("  model-conforming test data spurious peaks appear about half the")
        print("  time. Only read the periodicity when there is an excess to explain.")
    print("    lag " + "".join(f"{k:>6}" for k in range(1, 11)))
    print("    ac  " + "".join(f"{ac[k]:>6.2f}" for k in range(1, 11)))
    print("    lag " + "".join(f"{k:>6}" for k in range(11, 21)))
    print("    ac  " + "".join(f"{ac[k]:>6.2f}" for k in range(11, 21)))
    thresh = 2 / np.sqrt(len(r))
    strong = [k for k in range(1, 21) if abs(ac[k]) > thresh]
    print(f"    lags beyond 2/sqrt(n) = {thresh:.3f}: {strong if strong else 'none'}")
    peak, peaks, _ = period_from_acf(ac, len(r))
    if peak is not None:
        print(f"    -> significant peaks at lags {peaks}, fundamental {peak}"
              f" (+{ac[peak]:.2f}).")
        print(f"       If the excess is a sawtooth, its period is {peak}"
              " observations.")
        print("       Harmonics at multiples of the period are expected, and one"
              "\n       of them may be larger than the fundamental.")
        print("       Lag 1 can be low even so: consecutive steps are similar in"
              "\n       size to the rounding noise, which dilutes it.")
    if not strong:
        print("    -> no periodic structure at all. The excess is NOT a sawtooth;"
              "\n       look instead at the weekday fold and the plot below.")

    # -- fold by weekday and by day of month --------------------------------
    dow = np.array([x.weekday() for x in d])
    fold(r, dow, list(enumerate(["Mon", "Tue", "Wed", "Thu", "Fri"])),
         "by weekday (a pattern here means the crawl is applied on trading days,"
         "\n   not calendar days -- or that the model's t should not be calendar)")

    dom = np.array([(x.day - 1) // 7 for x in d])
    fold(r, dom, [(0, "days 1-7"), (1, "days 8-14"),
                  (2, "days 15-21"), (3, "days 22-28")],
         "by position in month (a ramp here means monthly stepping)")

    plot(r, [str(d[0]), str(d[-1])])
    return exc, f


# ---------------------------------------------------------------------------

def self_test():
    """Inject a known stepping period into synthetic data and try to recover it."""
    print("Injecting a KNOWN 7-day step into synthetic data.")
    print("The generator crawls continuously; we re-quantise the parity onto a")
    print("7-day grid, which is exactly the hypothesis being tested.\n")

    t, u, ls, sg, sdr, zar, _ = make_synthetic([(0.65, -4.80)], n_per=600, seed=9,
                                               quantise=False)
    # Step the crawl: hold the drift constant within each 7-day block.
    step = 7
    drift = (-4.80 / 100) * t
    blocks = (np.arange(len(t)) // step) * step
    stepped = (-4.80 / 100) * t[blocks]
    ls_step = ls - drift + stepped
    sdr2 = np.round(np.exp(ls_step), 4)
    zar2 = np.round(np.exp(ls_step + u), 4)
    ls2, u2 = np.log(sdr2), np.log(zar2) - np.log(sdr2)
    d = [date(2005, 6, 1)] * len(t)   # dates unused by the parts we check

    f = fit_window(t, u2, ls2, sdr2, zar2)
    exc = excess(f["rms_bp"], f["floor_bp"])
    # A block of 7 OBSERVATIONS is 7 weekdays, which is ~9.8 CALENDAR days, and
    # sawtooth_rms is in calendar days. Getting this wrong understates the
    # prediction by 40 percent -- the same calendar-vs-trading-day trap that
    # PROGRESS.md 4.4 flags for the crawl itself.
    T_cal = step * float(np.diff(t).mean()) * 365
    pred = sawtooth_rms(f["crawl"], T_cal)
    ac = autocorr(f["resid_bp"])
    peak, peaks, _ = period_from_acf(ac, len(f["resid_bp"]))
    print(f"  excess {exc:.2f} bp   predicted for a {T_cal:.1f}-calendar-day "
          f"step {pred:.2f} bp   ({100*abs(exc-pred)/pred:.0f}% apart)")
    print(f"  significant autocorrelation peaks at lags {peaks}"
          f" -> fundamental {peak}; lag 1 is only {ac[1]:+.2f}")
    print("  Lag 1 stays low by construction: consecutive sawtooth steps are"
          "\n  about the same size as the rounding noise, so the signal shows"
          "\n  at the PERIOD lag, not at lag 1. Look for the peak, not lag 1.")
    ok = abs(exc - pred) / pred < 0.25 and peak == step
    print("\n  PASS - a real stepping artefact is detectable this way"
          if ok else "\n  FAIL")
    return ok


def eras(args):
    """Every interval between announced crawl changes, oldest first."""
    bounds = [
        ("2005-06-01", "2006-06-29", "2005-06 .. 2006-06   announced -4.80"),
        ("2006-07-03", "2007-06-29", "2006-07 .. 2007-06   announced -3.90"),
        ("2007-07-02", "2009-02-27", "2007-07 .. 2009-02   announced -2.30"),
        ("2009-03-02", "2010-03-31", "2009-03 .. 2010-03   announced -2.91"),
        ("2010-04-01", "2012-05-31", "2010-04 .. 2012-05   announced -2.61"),
        ("2012-06-01", "2014-12-31", "2012-06 .. 2014-12   announced -0.16"),
        ("2020-05-04", "2022-12-30", "2020-05 .. 2022-12   announced -2.87"),
        ("2023-01-02", "2025-06-30", "2023-01 .. 2025-06   announced -1.51"),
        ("2025-07-31", "2026-09-15", "2025-07 .. now       announced -2.76"),
    ]
    print(f"  {'interval':<38}{'n':>6}{'crawl':>8}{'RMS':>8}"
          f"{'floor':>8}{'excess':>9}{'chi2/dof':>10}")
    for lo, hi, lab in bounds:
        try:
            d, t, u, ls, sg, sdr, zar = load(args.data, date.fromisoformat(lo),
                                             date.fromisoformat(hi))
        except Exception as e:
            print(f"  {lab:<38}  skipped ({e})")
            continue
        if len(t) < 60:
            print(f"  {lab:<38}  skipped (only {len(t)} obs)")
            continue
        f = fit_window(t, u, ls, sdr, zar)
        print(f"  {lab:<38}{len(t):>6}{f['crawl']:>8.2f}{f['rms_bp']:>8.2f}"
              f"{f['floor_bp']:>8.2f}{excess(f['rms_bp'], f['floor_bp']):>9.2f}"
              f"{f['chi2_dof']:>10.2f}")
    print("\n  Each interval has NO announced change inside it, so the excess is"
          "\n  what the straight-line model fails to explain within one policy regime.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--from", dest="start", type=date.fromisoformat)
    ap.add_argument("--to", dest="end", type=date.fromisoformat)
    ap.add_argument("--eras", action="store_true")
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--data", default="data/rates.csv")
    a = ap.parse_args()

    if a.self_test:
        sys.exit(0 if self_test() else 1)
    if a.eras:
        return eras(a)

    d, t, u, ls, sg, sdr, zar = load(a.data, a.start, a.end)
    diagnose(d, t, u, ls, sdr, zar, label="residual diagnostic")


if __name__ == "__main__":
    main()
