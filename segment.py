#!/usr/bin/env python3
"""
Tekanyetso — segmentation.  Week 3: naive optimal partitioning (exact, O(n^2)).

Splits the history into segments minimising total within-segment misfit, with a
fixed penalty per extra segment.  Within a segment the model is the one
`estimate.py` fits:

    ls = a + b*t - w*u + e            u = lz - ls,   t in YEARS

A changepoint is a date where (w, b) change -- a basket reweighting, a change in
the rate of crawl, or both at once.  STATUS.md §4 found the Bank does both at the
same review, so "both at once" is the common case, not the exotic one.

    python3 segment.py --self-test      # O(1) cost == lstsq; DP == exhaustive search
    python3 segment.py --synthetic      # known breaks, recovered?
    python3 segment.py --sweep          # penalty sweep on real data
    python3 segment.py                  # segment data/rates.csv
    python3 segment.py --from 2025-01-01 --to 2026-03-01 --min-seg 20   # zoom

---------------------------------------------------------------------------
DESIGN -- four decisions doing the real work
---------------------------------------------------------------------------

1.  The segment cost is O(1), not O(m).
    A DP that refits an OLS per candidate segment is O(n^3) and will not run on
    6,324 points.  Instead precompute cumulative sums of the 3x3 matrix X'X, the
    3-vector X'y and the scalar y'y.  Any segment's normal equations are a
    difference of two cumulative sums, solved in closed form.  The DP touches the
    raw data exactly once, at construction.

2.  Sigma is known in advance -- and it is NOT constant.
    PROGRESS.md §5 derives the noise on the basket index from BoB's own 4-decimal
    rounding, from first principles, before seeing a residual.  Having a
    physically-derived noise floor is a real luxury; most changepoint problems
    must estimate it from the same residuals they are testing.
    But §5 evaluates it at today's SDR of 0.0557 and treats the resulting 2.59 bp
    as a constant.  It is not.  The rounding grid is fixed at 1e-4 in ABSOLUTE
    terms while the printed SDR level moves, so the noise in log terms goes as
    1/SDR: about 1.6 bp in 2001 against 2.6 bp today.  A homoscedastic cost would
    therefore over-detect breaks in recent data and under-detect in early data.
    Each observation is weighted by its own 1/sigma_i instead.  The cost is then
    a true chi-square, so the BIC penalty is a clean 3*log(n) and chi2/dof has a
    direct reading: 1.0 means the peg is holding exactly to measurement precision.

3.  Every segment gets its own free intercept.
    Structural, not lazy.  At a weight change the modelled quantity changes
    definition -- w*lz + (1-w)*ls at w = 0.45 is a different index from the same
    expression at w = 0.50 -- so its level jumps mechanically even though the
    Pula itself does not.  Forcing continuity would model an artefact of the
    parameterisation.  It also keeps segments statistically independent, which is
    what PELT's pruning argument will require in week 4.

4.  DP and PELT will share one code path.
    Both reduce to "for each tau, evaluate the cost over an array of candidate
    starts".  Only the array differs: DP uses all of them, PELT uses a pruned
    set.  Writing a vectorised DP and a scalar PELT would make next week's
    runtime plot meaningless, so the primary metric is cost EVALUATIONS, with
    wall clock reported alongside.
"""

import argparse
import csv
import itertools
import sys
import time
from datetime import date

import numpy as np

P = 3                # parameters per segment: intercept, crawl, weight
GRID = 1e-4          # BoB prints every column to 4 decimals
CLEAN = "data/rates.csv"


def noise_sd(sdr, zar, w=0.5):
    """Per-observation sd of the basket-index residual, in log units.

    PROGRESS.md §5 made per-observation.  A value rounded to a grid of width h
    carries error of sd h/sqrt(12); in logs what matters is the RELATIVE size,
    h/sqrt(12)/level.  The two columns enter with weights (1-w) and w and their
    rounding errors are independent, so they add in quadrature.
    """
    s = GRID / np.sqrt(12)
    return np.hypot((1.0 - w) * s / sdr, w * s / zar)


# ---------------------------------------------------------------------------
# Segment cost
# ---------------------------------------------------------------------------

class SegmentCost:
    """Weighted RSS of `ls = a + b*t - w*u` on any interval [i, j), in O(1)."""

    def __init__(self, t_years, u, y, sigma):
        self.n = n = len(y)
        # Centring is free: shifting t, u or y only moves the per-segment
        # intercept, which is estimated separately in every segment, so w and b
        # are untouched.  It buys back several digits in the cumulative sums,
        # which is the one place this scheme is fragile.
        t = t_years - t_years.mean()
        u = u - u.mean()
        y = y - y.mean()
        self.t, self.u, self.y, self.sigma = t, u, y, sigma
        self.X = np.column_stack([np.ones(n), t, u])

        # Whitening: divide through by sigma_i so ordinary least squares on the
        # transformed data is weighted least squares on the original, and the
        # resulting RSS is a chi-square with no further scaling.
        Xw = self.X / sigma[:, None]
        yw = y / sigma

        # Preconditioning -- exact, not approximate.
        #
        # A segment's RSS is recovered as y'y - beta'X'y, and y'y runs ~6e4 times
        # larger than the RSS being extracted from it, so five digits are lost to
        # cancellation before the DP sees the number.  Fix: subtract one global
        # fit and segment the residual.  Every segment fit already spans the
        # columns of X, so subtracting X @ b_global shifts that segment's beta by
        # exactly -b_global and leaves its residual, hence its RSS, bit-for-bit
        # identical.  Same problem, two orders of magnitude less cancellation.
        # b_global is added back on the way out, so reported w and crawl are
        # unchanged.
        self.b_global, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
        yw = yw - Xw @ self.b_global

        cxx = np.zeros((n + 1, P, P))
        cxx[1:] = np.cumsum(Xw[:, :, None] * Xw[:, None, :], axis=0)
        cxy = np.zeros((n + 1, P))
        cxy[1:] = np.cumsum(Xw * yw[:, None], axis=0)
        cyy = np.zeros(n + 1)
        cyy[1:] = np.cumsum(yw * yw)
        self.cxx, self.cxy, self.cyy = cxx, cxy, cyy

        self.n_evals = 0          # the fair DP-vs-PELT metric

    def cost_many(self, starts, end):
        """chi-square of each segment [s, end).  Returns (chi2, beta)."""
        starts = np.asarray(starts)
        self.n_evals += starts.size

        M = self.cxx[end] - self.cxx[starts]
        v = self.cxy[end] - self.cxy[starts]
        yy = self.cyy[end] - self.cyy[starts]

        m00, m01, m02 = M[:, 0, 0], M[:, 0, 1], M[:, 0, 2]
        m11, m12, m22 = M[:, 1, 1], M[:, 1, 2], M[:, 2, 2]

        # Adjugate of a symmetric 3x3, closed form -- no batched LAPACK call.
        a00 = m11 * m22 - m12 * m12
        a01 = m02 * m12 - m01 * m22
        a02 = m01 * m12 - m02 * m11
        a11 = m00 * m22 - m02 * m02
        a12 = m01 * m02 - m00 * m12
        a22 = m00 * m11 - m01 * m01
        det = m00 * a00 + m01 * a01 + m02 * a02

        scale = np.abs(m00 * m11 * m22)          # m00 is the segment's weight mass
        bad = ~(np.abs(det) > 1e-10 * np.maximum(scale, 1e-300))
        safe = np.where(bad, 1.0, det)

        b0 = (a00 * v[:, 0] + a01 * v[:, 1] + a02 * v[:, 2]) / safe
        b1 = (a01 * v[:, 0] + a11 * v[:, 1] + a12 * v[:, 2]) / safe
        b2 = (a02 * v[:, 0] + a12 * v[:, 1] + a22 * v[:, 2]) / safe

        chi2 = yy - (b0 * v[:, 0] + b1 * v[:, 1] + b2 * v[:, 2])
        chi2 = np.maximum(chi2, 0.0)
        chi2[bad] = np.inf
        beta = np.column_stack([b0, b1, b2]) + self.b_global
        return chi2, beta

    def fit(self, i, j):
        """Readable single-segment summary.  Direct WLS -- O(m), reporting only."""
        X, y, sg = self.X[i:j], self.y[i:j], self.sigma[i:j]
        b, *_ = np.linalg.lstsq(X / sg[:, None], y / sg, rcond=None)
        resid = y - X @ b
        dof = max((j - i) - P, 1)
        chi2 = float(((resid / sg) ** 2).sum())
        return {
            "i": i, "j": j, "n": j - i,
            "w": float(-b[2]),
            "crawl_pct": float(b[1] * 100),          # t already in years
            "rmse_bp": float(1e4 * np.sqrt((resid ** 2).mean())),
            "floor_bp": float(1e4 * sg.mean()),
            "chi2_dof": chi2 / dof,
        }


# ---------------------------------------------------------------------------
# Naive optimal partitioning -- exact, O(n^2) cost evaluations
# ---------------------------------------------------------------------------

def segment_dp(cost, penalty, min_seg):
    """Exact segmentation.  F[tau] = cost of the best partition of [0, tau)."""
    n = cost.n
    F = np.full(n + 1, np.inf)
    F[0] = -penalty                        # so the first segment is not charged
    prev = np.zeros(n + 1, dtype=int)

    for tau in range(min_seg, n + 1):
        starts = np.arange(0, tau - min_seg + 1)
        live = np.isfinite(F[starts])
        if not live.any():
            continue
        starts = starts[live]
        chi2, _ = cost.cost_many(starts, tau)
        total = F[starts] + chi2 + penalty
        k = int(np.argmin(total))
        F[tau] = total[k]
        prev[tau] = starts[k]

    cuts, tau = [], n
    while tau > 0:
        cuts.append(tau)
        tau = int(prev[tau])
    return sorted(cuts)[:-1], float(F[n])


def segment_brute(cost, penalty, min_seg, max_k=3):
    """Exhaustive search over partitions.  Tiny n only -- a correctness oracle."""
    n = cost.n
    best, best_cuts = np.inf, None
    for k in range(0, max_k + 1):
        for cuts in itertools.combinations(range(min_seg, n - min_seg + 1), k):
            bnds = (0,) + cuts + (n,)
            if any(bnds[i + 1] - bnds[i] < min_seg for i in range(len(bnds) - 1)):
                continue
            tot = -penalty
            for i in range(len(bnds) - 1):
                c, _ = cost.cost_many(np.array([bnds[i]]), bnds[i + 1])
                tot += c[0] + penalty
            if tot < best:
                best, best_cuts = tot, list(cuts)
    return best_cuts, best


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

def load(path=CLEAN, start=None, end=None):
    with open(path) as f:
        rows = sorted(csv.DictReader(f), key=lambda r: r["date"])
    d = [date.fromisoformat(r["date"]) for r in rows]
    keep = [k for k, x in enumerate(d)
            if (start is None or x >= start) and (end is None or x <= end)]
    d = np.array([d[k] for k in keep])
    zar = np.array([float(rows[k]["ZAR"]) for k in keep])
    sdr = np.array([float(rows[k]["SDR"]) for k in keep])
    t = np.array([(x - d[0]).days for x in d], dtype=float) / 365.0
    lz, ls = np.log(zar), np.log(sdr)
    return d, t, lz - ls, ls, noise_sd(sdr, zar)


def make_synthetic(segments, n_per=400, seed=0, quantise=True):
    """Data obeying the peg exactly, then rounded the way BoB rounds it.

    `segments` is a list of (w, crawl_pct_pa).  The spread u is a random walk:
    u = log(ZAR/SDR) has the Pula cancelled out of it (PROGRESS.md §4.3), so it
    is a pure world-market quantity and simulating it as exogenous is faithful.
    """
    rng = np.random.default_rng(seed)
    n = n_per * len(segments)
    t = np.arange(n) / 365.0 * (7 / 5)                 # weekdays -> calendar years
    u = 3.09 + np.cumsum(rng.normal(0, 0.006, n))

    ls = np.empty(n)
    a = np.log(0.0557) + 0.45 * 3.09
    truth = []
    for k, (w, crawl) in enumerate(segments):
        s, e = k * n_per, (k + 1) * n_per
        b = crawl / 100.0
        if k > 0:
            # Choose `a` so the Pula's own value does not jump at the boundary.
            # The index jumps -- it has changed definition -- the currency does not.
            a = ls[s - 1] - (b * t[s] - w * u[s])
        ls[s:e] = a + b * t[s:e] - w * u[s:e]
        truth.append({"start": s, "w": w, "crawl": crawl})
    lz = ls + u

    sdr, zar = np.exp(ls), np.exp(lz)
    if quantise:
        # The actual noise mechanism, not an additive Gaussian stand-in.
        sdr, zar = np.round(sdr, 4), np.round(zar, 4)
        ls, lz = np.log(sdr), np.log(zar)
    return t, lz - ls, ls, noise_sd(sdr, zar), truth


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def report(cost, cuts, dates=None):
    bnds = [0] + list(cuts) + [cost.n]
    print(f"  {'segment':<26}{'n':>6}{'w':>8}{'crawl %/yr':>12}"
          f"{'resid bp':>10}{'floor bp':>10}{'chi2/dof':>10}")
    for i in range(len(bnds) - 1):
        f = cost.fit(bnds[i], bnds[i + 1])
        lab = (f"{dates[bnds[i]]} .. {dates[bnds[i+1]-1]}" if dates is not None
               else f"[{bnds[i]}, {bnds[i+1]})")
        print(f"  {lab:<26}{f['n']:>6}{f['w']:>8.3f}{f['crawl_pct']:>12.2f}"
              f"{f['rmse_bp']:>10.2f}{f['floor_bp']:>10.2f}{f['chi2_dof']:>10.2f}")


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def self_test():
    ok = True

    print("1. O(1) cumulative-sum cost == direct weighted least squares")
    print("   Tolerance is in COST units, not relative.  The DP only ever compares")
    print("   costs against the penalty, so an error is harmless until it is a")
    print("   material fraction of one penalty term.")
    t, u, ls, sg, _ = make_synthetic([(0.45, -1.51), (0.50, -2.76)], n_per=500, seed=1)
    c = SegmentCost(t, u, ls, sg)
    pen = P * np.log(c.n)
    rng = np.random.default_rng(7)
    wc = wb = 0.0
    for _ in range(300):
        i = rng.integers(0, c.n - 40)
        j = rng.integers(i + 40, c.n + 1)
        fast, beta = c.cost_many(np.array([i]), j)
        Xw = c.X[i:j] / c.sigma[i:j, None]
        yw = c.y[i:j] / c.sigma[i:j]
        bb, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
        slow = float(((yw - Xw @ bb) ** 2).sum())
        wc = max(wc, abs(fast[0] - slow))
        wb = max(wb, float(np.abs(beta[0] - bb).max()))
    print(f"   worst cost error over 300 random segments: {wc:.2e}"
          f"  ({wc / pen:.1e} of one penalty term)")
    print(f"   worst coefficient error:                   {wb:.2e}")
    ok &= wc < 1e-3 * pen and wb < 1e-9
    print("   PASS\n" if wc < 1e-3 * pen and wb < 1e-9 else "   FAIL\n")

    print("2. dynamic program == exhaustive search on a small series")
    t, u, ls, sg, _ = make_synthetic([(0.45, -1.51), (0.50, -2.76)], n_per=60, seed=3)
    c = SegmentCost(t, u, ls, sg)
    pen = P * np.log(c.n)
    dc, dv = segment_dp(c, pen, min_seg=20)
    bc, bv = segment_brute(c, pen, min_seg=20, max_k=3)
    print(f"   DP        : cuts={dc} cost={dv:.4f}")
    print(f"   exhaustive: cuts={bc} cost={bv:.4f}")
    good = dc == bc and abs(dv - bv) < 1e-6
    ok &= good
    print("   PASS\n" if good else "   FAIL -- the DP is not exact\n")

    print("3. a series with no break gets no break (5 seeds)")
    bad = []
    for seed in range(5):
        t, u, ls, sg, _ = make_synthetic([(0.50, -2.76)], n_per=600, seed=seed)
        c = SegmentCost(t, u, ls, sg)
        cuts, _ = segment_dp(c, P * np.log(c.n), min_seg=90)
        if cuts:
            bad.append((seed, cuts))
    print(f"   spurious breaks: {bad if bad else 'none'}")
    ok &= not bad
    print("   PASS\n" if not bad else "   FAIL\n")

    print("4. chi2/dof == 1 when sigma is right")
    t, u, ls, sg, _ = make_synthetic([(0.50, -2.76)], n_per=1500, seed=2)
    c = SegmentCost(t, u, ls, sg)
    f = c.fit(0, c.n)
    print(f"   chi2/dof = {f['chi2_dof']:.3f}   resid {f['rmse_bp']:.2f} bp"
          f"   vs floor {f['floor_bp']:.2f} bp")
    good = 0.85 < f["chi2_dof"] < 1.15
    ok &= good
    print("   PASS\n" if good else "   FAIL -- sigma is misspecified\n")
    return ok


def synthetic_run(min_seg, mult):
    segs = [(0.55, -0.16), (0.55, 0.00), (0.50, 0.26), (0.45, -1.51), (0.50, -2.76)]
    t, u, ls, sg, truth = make_synthetic(segs, n_per=400, seed=11)
    c = SegmentCost(t, u, ls, sg)
    pen = mult * P * np.log(c.n)
    print(f"n={c.n}, true breaks at {[x['start'] for x in truth[1:]]}, "
          f"penalty={pen:.1f} ({mult}x BIC), min_seg={min_seg}")
    print("note: break 1 is crawl-only, and only 0.16 pp of it -- 0.04 bp/day under\n"
          "      ~3.6 bp of rounding noise.  That is the hard case.\n")

    t0 = time.perf_counter()
    cuts, _ = segment_dp(c, pen, min_seg)
    dt = time.perf_counter() - t0
    print(f"found {len(cuts)} breaks in {dt:.1f}s, {c.n_evals:,} cost evaluations")
    print(f"cuts: {list(map(int, cuts))}\n")
    report(c, cuts)

    print("\n  true break -> nearest detected (error in observations)")
    for x in truth[1:]:
        if cuts:
            k = int(min(cuts, key=lambda z: abs(z - x["start"])))
            print(f"  {x['start']:>6} -> {k:>6}   ({k - x['start']:+d})"
                  f"   w {x['w']}, crawl {x['crawl']}")
        else:
            print(f"  {x['start']:>6} -> none")
    extra = [int(k) for k in cuts
             if min(abs(k - x["start"]) for x in truth[1:]) > 10]
    print(f"  false positives (>10 obs from any true break): {extra or 'none'}")


def sweep(cost, min_seg):
    print(f"  {'mult':>6}{'penalty':>10}{'breaks':>8}")
    for mult in (1, 2, 5, 10, 20, 50, 100, 200, 500):
        cuts, _ = segment_dp(cost, mult * P * np.log(cost.n), min_seg)
        print(f"  {mult:>6}{mult * P * np.log(cost.n):>10.0f}{len(cuts):>8}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--synthetic", action="store_true")
    ap.add_argument("--sweep", action="store_true")
    ap.add_argument("--min-seg", type=int, default=90,
                    help="minimum segment length.  The crawl moves ~0.8 bp/day "
                         "under ~2.6 bp of rounding noise, so short segments "
                         "cannot identify it and will fit artefacts instead.")
    ap.add_argument("--mult", type=float, default=1.0,
                    help="penalty as a multiple of BIC, 3*log(n)")
    ap.add_argument("--from", dest="start", type=date.fromisoformat)
    ap.add_argument("--to", dest="end", type=date.fromisoformat)
    ap.add_argument("--data", default=CLEAN)
    a = ap.parse_args()

    if a.self_test:
        sys.exit(0 if self_test() else 1)
    if a.synthetic:
        return synthetic_run(a.min_seg, a.mult)

    d, t, u, ls, sg = load(a.data, a.start, a.end)
    c = SegmentCost(t, u, ls, sg)
    print(f"{c.n} observations, {d[0]} .. {d[-1]}")
    print(f"rounding floor: {1e4*sg[0]:.2f} bp at the start, "
          f"{1e4*sg[-1]:.2f} bp at the end\n")

    if a.sweep:
        return sweep(c, a.min_seg)

    pen = a.mult * P * np.log(c.n)
    t0 = time.perf_counter()
    cuts, _ = segment_dp(c, pen, a.min_seg)
    dt = time.perf_counter() - t0
    print(f"penalty {pen:.1f} ({a.mult}x BIC), min_seg {a.min_seg}")
    print(f"{len(cuts)} breaks in {dt:.1f}s, {c.n_evals:,} cost evaluations\n")
    report(c, cuts, d)

    with open("segments.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["start_date", "end_date", "n", "w", "crawl_pct",
                    "rmse_bp", "floor_bp", "chi2_dof"])
        bnds = [0] + list(cuts) + [c.n]
        for i in range(len(bnds) - 1):
            f = c.fit(bnds[i], bnds[i + 1])
            w.writerow([d[bnds[i]], d[bnds[i + 1] - 1], f["n"],
                        f"{f['w']:.4f}", f"{f['crawl_pct']:.3f}",
                        f"{f['rmse_bp']:.2f}", f"{f['floor_bp']:.2f}",
                        f"{f['chi2_dof']:.3f}"])
    print("\nwrote segments.csv")


if __name__ == "__main__":
    main()
