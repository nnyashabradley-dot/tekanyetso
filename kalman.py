#!/usr/bin/env python3
"""
Tekanyetso -- Kalman filter for a time-varying peg.

segment.py says the weight w and the crawl b are constant inside a regime and
jump at changepoints. This says something different: they drift, and the best
estimate of them is updated every day as each new rate arrives. Same
quantities, different assumption about how they move.

State
-----
    x_t = [ m_t, b, w_t ]

    m_t  the parity: the log SDR-per-Pula the rule implies at the previous
         day's rand/SDR spread
    b    crawl per step (a trading day or a calendar day, see --dt)
    w_t  weight on the rand

Why m and not the basket index L = ls + w*u.  When the Bank reweights, the
Pula itself does not jump; the index does, by dw*u ~ 0.05*3.1 ~ 15 percent,
because it has changed definition. With L in the state, every reweighting
would need a 1,500 bp jump in L exactly correlated with the jump in w. With m,
a reweighting leaves m continuous and only w moves.

Transition (linear, time-varying):
    m_t = m_{t-1} + b*dt_t - w_{t-1}*(u_{t-1} - u_{t-2}) + eta_m
    b   = b
    w_t = w_{t-1} + eta_w
    (b also gets a small random walk eta_b)

Observation:
    ls_t = m_t - w_t*(u_t - u_{t-1}) + e_t,   Var(e_t) = floor_t^2

floor_t is the rounding noise derived in PROGRESS.md section 5, per
observation and at the current weight, so it is known rather than estimated.
Only the three process variances are fitted, by maximum likelihood.

Numerics: the covariance update is in Joseph form,
    P = (I-KH) P (I-KH)' + K R K'
which stays symmetric and positive semi-definite under rounding, unlike the
textbook P = (I-KH)P.

Usage
-----
    python3 kalman.py --self-test
    python3 kalman.py                    # fit on 2005-06..2010-12, filter all
    python3 kalman.py --compare-dt       # trading-day vs calendar-day steps
"""

import argparse
import csv
import sys
from datetime import date

import numpy as np

from segment import GRID, load, make_synthetic, noise_sd

FIT_FROM = date(2005, 6, 1)
FIT_TO = date(2010, 12, 31)       # hyperparameters see only this, then freeze
INIT_N = 90                       # observations used to initialise the state


# ---------------------------------------------------------------------------

def steps(d, t_years, mode):
    """dt for each observation. 'trading': 1 per published day. 'calendar':
    calendar days elapsed, so b is per calendar day."""
    if mode == "trading":
        dt = np.ones(len(t_years))
    else:
        dt = np.diff(t_years, prepend=t_years[0]) * 365.0
    dt[0] = 0.0
    return dt


def init_state(ls, u, dt, sdr, zar, n0=INIT_N):
    """OLS on the first n0 observations. Returns x0, P0 at index n0-1."""
    tt = np.cumsum(dt[:n0])
    X = np.column_stack([np.ones(n0), tt, -u[:n0]])
    beta, *_ = np.linalg.lstsq(X, ls[:n0], rcond=None)
    resid = ls[:n0] - X @ beta
    s2 = max(float(resid.var()), 1e-12)
    cov = s2 * np.linalg.inv(X.T @ X)
    a, b, w = beta
    L = a + b * tt[-1]
    m = L - w * u[n0 - 2]                       # m is defined at u_{t-1}
    x0 = np.array([m, b, w])
    # Jacobian of (m, b, w) wrt (a, b, w)
    J = np.array([[1.0, tt[-1], -u[n0 - 2]], [0, 1, 0], [0, 0, 1]])
    P0 = J @ cov @ J.T * 25.0                   # inflate: OLS is overconfident
    return x0, P0


def kfilter(ls, u, dt, sdr, zar, q, x0, P0, start):
    """Run the filter from index start (state x0, P0 given at start-1).

    Returns per-index arrays over the whole series (NaN before start):
      pred_ls[i]   one-step prediction of ls[i] made at i-1, using u[i]
      pred_var[i]  its variance
      xf[i]        filtered state after seeing i
      loglik       sum over i >= start
    Everything at index i depends only on data at indices <= i.
    """
    n = len(ls)
    qm, qb, qw = q
    Q = np.diag([qm, qb, qw])
    x, P = x0.copy(), P0.copy()
    I3 = np.eye(3)
    pred_ls = np.full(n, np.nan)
    pred_var = np.full(n, np.nan)
    xf = np.full((n, 3), np.nan)
    xf[start - 1] = x
    ll = 0.0
    for i in range(start, n):
        du_prev = u[i - 1] - u[i - 2]
        F = np.array([[1.0, dt[i], -du_prev], [0, 1, 0], [0, 0, 1]])
        x = F @ x
        P = F @ P @ F.T + Q
        du = u[i] - u[i - 1]
        H = np.array([1.0, 0.0, -du])
        w_now = min(max(x[2], 0.0), 1.0)
        r = noise_sd(sdr[i:i + 1], zar[i:i + 1], w_now)[0] ** 2
        yhat = H @ x
        S = H @ P @ H + r
        pred_ls[i], pred_var[i] = yhat, S
        v = ls[i] - yhat
        ll += -0.5 * (np.log(2 * np.pi * S) + v * v / S)
        K = P @ H / S
        x = x + K * v
        A = I3 - np.outer(K, H)
        P = A @ P @ A.T + r * np.outer(K, K)     # Joseph form
        xf[i] = x
    return pred_ls, pred_var, xf, ll


def nelder_mead(f, x0, step=1.0, iters=300, tol=1e-6):
    """Minimal Nelder-Mead, so the project keeps numpy as its only dependency."""
    n = len(x0)
    pts = [np.array(x0, float)] + [np.array(x0, float) + step * np.eye(n)[k]
                                     for k in range(n)]
    vals = [f(p) for p in pts]
    for _ in range(iters):
        order = np.argsort(vals)
        pts = [pts[k] for k in order]
        vals = [vals[k] for k in order]
        if abs(vals[-1] - vals[0]) < tol * (1 + abs(vals[0])):
            break
        c = np.mean(pts[:-1], axis=0)
        xr = c + (c - pts[-1]); fr = f(xr)
        if fr < vals[0]:
            xe = c + 2 * (c - pts[-1]); fe = f(xe)
            pts[-1], vals[-1] = (xe, fe) if fe < fr else (xr, fr)
        elif fr < vals[-2]:
            pts[-1], vals[-1] = xr, fr
        else:
            xc = c + 0.5 * (pts[-1] - c); fc = f(xc)
            if fc < vals[-1]:
                pts[-1], vals[-1] = xc, fc
            else:
                for k in range(1, n + 1):
                    pts[k] = pts[0] + 0.5 * (pts[k] - pts[0])
                    vals[k] = f(pts[k])
    k = int(np.argmin(vals))
    return pts[k], vals[k]


def fit_q(ls, u, dt, sdr, zar, n0=INIT_N, verbose=True):
    """Maximum likelihood for the three process variances (log scale)."""
    x0, P0 = init_state(ls, u, dt, sdr, zar, n0)

    def nll(theta):
        q = np.exp(theta)
        return -kfilter(ls, u, dt, sdr, zar, q, x0, P0, n0)[3]

    b_scale = abs(x0[1]) if abs(x0[1]) > 0 else 1e-4
    start = np.log([(0.5e-4) ** 2, (0.01 * b_scale) ** 2, 0.002 ** 2])
    theta, val = nelder_mead(nll, start, step=2.0, iters=250)
    theta, val = nelder_mead(nll, theta, step=0.5, iters=250)
    if verbose:
        q = np.exp(theta)
        print(f"  fitted: sd(eta_m) {1e4*np.sqrt(q[0]):.3f} bp/step   "
              f"sd(eta_w) {np.sqrt(q[2]):.5f}/step   "
              f"sd(eta_b) {np.sqrt(q[1]):.2e}/step   loglik {-val:.1f}")
    return np.exp(theta), -val


# ---------------------------------------------------------------------------

def prepare(path, start, end=None):
    d, t, u, ls, sg, sdr, zar = load(path, start, end)
    return d, t, u, ls, sdr, zar


def run(path, mode="trading", fit_to=FIT_TO, verbose=True):
    """Fit q on [FIT_FROM, fit_to] only, then filter the whole post-2005 series."""
    d, t, u, ls, sdr, zar = prepare(path, FIT_FROM)
    dt = steps(d, t, mode)
    nfit = int(np.searchsorted(d, fit_to, side="right"))
    if verbose:
        print(f"fitting process noise on {d[0]} .. {d[nfit-1]} ({nfit} obs), "
              f"steps = {mode}")
    q, _ = fit_q(ls[:nfit], u[:nfit], dt[:nfit], sdr[:nfit], zar[:nfit],
                 verbose=verbose)
    x0, P0 = init_state(ls, u, dt, sdr, zar)
    pred_ls, pred_var, xf, ll = kfilter(ls, u, dt, sdr, zar, q, x0, P0, INIT_N)
    return {"d": d, "t": t, "u": u, "ls": ls, "sdr": sdr, "zar": zar, "dt": dt,
            "q": q, "pred_ls": pred_ls, "pred_var": pred_var, "xf": xf,
            "loglik": ll, "mode": mode}


def annual_crawl(xf, d, mode):
    """Filtered crawl in %/yr. Trading steps are annualised by the observed
    number of publications per year in the trailing 252 observations."""
    b = xf[:, 1]
    if mode == "calendar":
        return b * 365 * 100
    t = np.array([(x - d[0]).days for x in d], float)
    out = np.full(len(b), np.nan)
    for i in range(len(b)):
        j = max(0, i - 252)
        span = (t[i] - t[j]) / 365.0
        out[i] = b[i] * ((i - j) / span if span > 0 else 250.0) * 100
    return out


def compare_dt(path):
    """Same model, same data, same number of parameters; only the clock
    differs. The higher likelihood says which clock the Bank steps on."""
    print("Trading-day steps against calendar-day steps, full post-2005 history.")
    print("Process noise is re-fitted by maximum likelihood for each.\n")
    out = {}
    d, t, u, ls, sdr, zar = prepare(path, FIT_FROM)
    for mode in ("trading", "calendar"):
        dt = steps(d, t, mode)
        print(f"  {mode}:")
        q, ll = fit_q(ls, u, dt, sdr, zar)
        out[mode] = ll
    diff = out["trading"] - out["calendar"]
    print(f"\n  log-likelihood, trading minus calendar: {diff:+.1f}")
    print(f"  likelihood ratio: e^{diff:.1f}. The models have the same number of"
          "\n  parameters, so this is a direct comparison; a difference above ~5"
          "\n  is decisive by the usual Bayes-factor conventions.")
    return diff


def self_test():
    ok = True
    print("1. recovers a known weight and crawl, and tracks a reweighting")
    segs = [(0.65, -4.80), (0.60, -2.30)]
    t, u, ls, sg, sdr, zar, _ = make_synthetic(segs, n_per=600, seed=4)
    dt = np.ones(len(t)) * (7 / 5)    # synthetic data has weekday spacing in days
    dt[0] = 0
    q, _ = fit_q(ls, u, dt, sdr, zar, verbose=False)
    x0, P0 = init_state(ls, u, dt, sdr, zar)
    _, _, xf, _ = kfilter(ls, u, dt, sdr, zar, q, x0, P0, INIT_N)
    w1, w2 = xf[550, 2], xf[1150, 2]
    c1, c2 = xf[550, 1] * 365 * 100, xf[1150, 1] * 365 * 100
    print(f"   w  before/after the break: {w1:.3f} / {w2:.3f}   (true 0.65 / 0.60)")
    print(f"   crawl before/after:        {c1:+.2f} / {c2:+.2f}   (true -4.80 / -2.30)")
    good = abs(w1 - 0.65) < 0.01 and abs(w2 - 0.60) < 0.01 \
        and abs(c1 + 4.8) < 0.5 and abs(c2 + 2.3) < 0.5
    ok &= good
    print("   PASS\n" if good else "   FAIL\n")

    print("2. Joseph form keeps the covariance symmetric positive definite")
    x0, P0 = init_state(ls, u, dt, sdr, zar)
    # run with tiny noise, where the textbook update is most fragile
    x, P = x0.copy(), P0.copy()
    worst = np.inf
    for i in range(INIT_N, len(ls)):
        F = np.array([[1.0, dt[i], -(u[i-1] - u[i-2])], [0, 1, 0], [0, 0, 1]])
        P = F @ P @ F.T + np.diag(q * 1e-6)
        H = np.array([1.0, 0.0, -(u[i] - u[i-1])])
        r = noise_sd(sdr[i:i+1], zar[i:i+1], 0.6)[0] ** 2
        K = P @ H / (H @ P @ H + r)
        A = np.eye(3) - np.outer(K, H)
        P = A @ P @ A.T + r * np.outer(K, K)
        Pn = P / np.sqrt(np.outer(np.diag(P), np.diag(P)))   # correlation form
        worst = min(worst, float(np.linalg.eigvalsh((Pn + Pn.T) / 2).min()))
    print(f"   smallest eigenvalue of the correlation matrix over the run: {worst:.2e}")
    good = worst > -1e-9
    ok &= good
    print("   PASS\n" if good else "   FAIL\n")

    print("3. causal: prediction at i unchanged when data after i is corrupted")
    x0, P0 = init_state(ls, u, dt, sdr, zar)
    p1 = kfilter(ls, u, dt, sdr, zar, q, x0, P0, INIT_N)[0]
    k = 700
    ls2, u2 = ls.copy(), u.copy()
    rng = np.random.default_rng(0)
    ls2[k + 1:] = rng.normal(0, 1, len(ls) - k - 1)
    u2[k + 1:] = rng.normal(0, 1, len(ls) - k - 1)
    p2 = kfilter(ls2, u2, dt, sdr, zar, q, x0, P0, INIT_N)[0]
    good = np.array_equal(p1[:k + 1], p2[:k + 1], equal_nan=True)
    print(f"   predictions up to index {k} bit-identical: {good}")
    ok &= good
    print("   PASS\n" if good else "   FAIL\n")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--compare-dt", action="store_true")
    ap.add_argument("--data", default="data/rates.csv")
    a = ap.parse_args()
    if a.self_test:
        sys.exit(0 if self_test() else 1)
    if a.compare_dt:
        return compare_dt(a.data)

    r = run(a.data)
    d, xf = r["d"], r["xf"]
    crawl = annual_crawl(xf, d, r["mode"])
    with open("kalman.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["date", "w", "crawl_pct", "parity_m", "pred_err_bp", "pred_sd_bp"])
        for i in range(len(d)):
            if np.isnan(xf[i, 0]):
                continue
            e = 1e4 * (r["ls"][i] - r["pred_ls"][i]) if not np.isnan(r["pred_ls"][i]) else ""
            s = 1e4 * np.sqrt(r["pred_var"][i]) if not np.isnan(r["pred_var"][i]) else ""
            w.writerow([d[i], f"{xf[i,2]:.5f}", f"{crawl[i]:.4f}", f"{xf[i,0]:.8f}",
                        f"{e:.3f}" if e != "" else "", f"{s:.3f}" if s != "" else ""])
    print("\n  filtered weight and crawl, one row per year-end:")
    print(f"    {'date':<12}{'w':>8}{'crawl %/yr':>12}")
    last_year = None
    for i in range(len(d)):
        if np.isnan(xf[i, 0]):
            continue
        if last_year is not None and d[i].year != last_year:
            j = i - 1
            print(f"    {str(d[j]):<12}{xf[j,2]:>8.3f}{crawl[j]:>12.2f}")
        last_year = d[i].year
    print(f"    {str(d[-1]):<12}{xf[-1,2]:>8.3f}{crawl[-1]:>12.2f}")
    err = 1e4 * (r["ls"] - r["pred_ls"])
    m = ~np.isnan(err)
    z = err[m] / (1e4 * np.sqrt(r["pred_var"][m]))
    print(f"\n  one-step prediction error: RMS {np.sqrt(np.mean(err[m]**2)):.2f} bp"
          f"   standardised RMS {np.sqrt(np.mean(z**2)):.2f} (1.00 = calibrated)")
    print("wrote kalman.csv")


if __name__ == "__main__":
    main()
