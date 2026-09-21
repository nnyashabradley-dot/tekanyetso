#!/usr/bin/env python3
"""
Tekanyetso -- PELT versus the naive dynamic program.

Two claims, both checked here:

  1. Identical answers.  Same breakpoints, same optimal cost, on synthetic data
     with known breaks and on the real Bank of Botswana history.
  2. Near-linear cost.  The primary metric is COST EVALUATIONS, because both
     methods share the same vectorised cost function and differ only in the
     array of candidate starts.  Wall clock is reported alongside.

    python3 bench_pelt.py                 # synthetic scaling + real data
    python3 bench_pelt.py --quick         # smaller sizes
    python3 bench_pelt.py --counterexample

Writes bench_pelt.csv for the runtime figure.
"""

import argparse
import csv
import time
from datetime import date

import numpy as np

from segment import (P, SegmentCost, load, make_synthetic, segment_dp,
                     segment_pelt)


def one(cost_factory, pen, m):
    """Run DP and PELT on fresh cost objects; return a result row."""
    c1 = cost_factory()
    t0 = time.perf_counter()
    cuts_dp, f_dp = segment_dp(c1, pen, m)
    t_dp = time.perf_counter() - t0

    c2 = cost_factory()
    t0 = time.perf_counter()
    cuts_pe, f_pe = segment_pelt(c2, pen, m)
    t_pe = time.perf_counter() - t0

    same = cuts_dp == cuts_pe and abs(f_dp - f_pe) < 1e-6
    return {
        "n": c1.n, "breaks": len(cuts_dp), "identical": same,
        "evals_dp": c1.n_evals, "evals_pelt": c2.n_evals,
        "sec_dp": t_dp, "sec_pelt": t_pe,
        "max_candidates": max(segment_pelt.last_sizes or [0]),
    }


def synthetic_scaling(sizes, mult, m, seed=3):
    """Breaks roughly every 300 observations -- the real data's density."""
    rows = []
    rng = np.random.default_rng(seed)
    for n in sizes:
        k = max(1, n // 300)
        segs = [(float(rng.choice([0.45, 0.5, 0.55, 0.6, 0.65])),
                 float(rng.uniform(-5, 0.5))) for _ in range(k)]
        t, u, ls, sg, sdr, zar, _ = make_synthetic(segs, n_per=n // k, seed=seed + n)
        f = lambda: SegmentCost(t, u, ls, sg, sdr, zar)
        pen = mult * P * np.log(len(t))
        r = one(f, pen, m)
        r.update(kind="synthetic", mult=mult, min_seg=m)
        rows.append(r)
        print(f"  n={r['n']:>6}  breaks={r['breaks']:>3}  "
              f"DP {r['evals_dp']:>11,} evals {r['sec_dp']:6.2f}s   "
              f"PELT {r['evals_pelt']:>9,} evals {r['sec_pelt']:5.2f}s   "
              f"x{r['evals_dp']/max(r['evals_pelt'],1):5.1f}   "
              f"identical={r['identical']}")
    return rows


def real(path, start, mult, m):
    d, t, u, ls, sg, sdr, zar = load(path, start)
    f = lambda: SegmentCost(t, u, ls, sg, sdr, zar)
    pen = mult * P * np.log(len(t))
    r = one(f, pen, m)
    r.update(kind="real", mult=mult, min_seg=m)
    print(f"  real {d[0]}..{d[-1]} n={r['n']}  mult={mult}  breaks={r['breaks']}  "
          f"DP {r['evals_dp']:,} ({r['sec_dp']:.1f}s)  PELT {r['evals_pelt']:,} "
          f"({r['sec_pelt']:.1f}s)  x{r['evals_dp']/r['evals_pelt']:.1f}  "
          f"identical={r['identical']}")
    return r


def counterexample():
    """Reproduce a case where pruning WITHOUT the min_seg delay is wrong."""
    rng = np.random.default_rng(1)
    for i in range(100000):
        n = int(rng.integers(30, 120))
        m = int(rng.integers(4, 25))
        t = np.cumsum(rng.uniform(0.5, 2, n)) / 365
        u = np.cumsum(rng.normal(0, 0.01, n))
        y = np.zeros(n)
        pos = 0
        while pos < n:
            L = int(rng.integers(3, 40))
            a, b, w = rng.normal(0, 0.01), rng.normal(0, 0.05), rng.uniform(0, 1)
            y[pos:pos + L] = a + b * t[pos:pos + L] - w * u[pos:pos + L]
            pos += L
        sig = np.full(n, 10 ** rng.uniform(-4.5, -3))
        y = y + rng.normal(0, 1, n) * sig
        pen = float(10 ** rng.uniform(-0.5, 2))
        c = SegmentCost(t, u, y, sig)
        a_, fa = segment_dp(c, pen, m)
        d_, fd = segment_pelt(c, pen, m, delay=False)
        if fd > fa + 1e-6:
            b_, fb = segment_pelt(c, pen, m, delay=True)
            print(f"case {i}: n={n}, min_seg={m}, penalty={pen:.3f}")
            print(f"  exact DP            cuts {a_}  cost {fa:.4f}")
            print(f"  PELT, no delay      cuts {d_}  cost {fd:.4f}   <-- worse")
            print(f"  PELT, delay min_seg cuts {b_}  cost {fb:.4f}")
            print("Without the delay, a start that had been pruned turned out to be")
            print("the optimum later. Why, and why a delay of min_seg is enough, is the")
            print("correctness argument.")
            return
    print("no counterexample found")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--counterexample", action="store_true")
    ap.add_argument("--data", default="data/rates.csv")
    a = ap.parse_args()

    if a.counterexample:
        return counterexample()

    sizes = [500, 1000, 2000, 4000] if a.quick else [500, 1000, 2000, 4000, 6000, 8000]
    rows = []
    for mult in (1, 10):
        print(f"\nsynthetic, penalty {mult}x BIC, min_seg 90")
        rows += synthetic_scaling(sizes, mult, 90)

    print("\nreal data, from 2005-06-01, min_seg 90")
    for mult in (1, 10):
        try:
            rows.append(real(a.data, date(2005, 6, 1), mult, 90))
        except FileNotFoundError:
            print("  data/rates.csv not found -- skipped")
            break

    with open("bench_pelt.csv", "w", newline="") as fh:
        keys = ["kind", "mult", "min_seg", "n", "breaks", "identical",
                "evals_dp", "evals_pelt", "sec_dp", "sec_pelt", "max_candidates"]
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in keys})
    print("\nwrote bench_pelt.csv")
    if not all(r["identical"] for r in rows):
        raise SystemExit("MISMATCH between PELT and DP")


if __name__ == "__main__":
    main()
