"""Synthetic check: build a series with a KNOWN w and crawl, rounded to 4dp
exactly as BoB rounds, and see whether fit() recovers them inside its own
stated standard errors."""
import sys, os, csv, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from estimate import fit

rng = np.random.default_rng(7)
W_TRUE, CRAWL_TRUE = 0.50, -2.76
N = 25

for N in (25, 90, 250):
    t = np.arange(N, dtype=float) * 1.4          # ~1.4 calendar days per obs
    u = np.cumsum(rng.normal(0, 0.006, N))       # rand vs SDR wanders
    a = np.log(0.0557)
    b = CRAWL_TRUE / 100 / 365
    ls_true = a + b * t - W_TRUE * u
    lz_true = ls_true + u
    sdr = np.round(np.exp(ls_true), 4)           # BoB's 4dp grid
    zar = np.round(np.exp(lz_true) * 22.1, 4)    # scale into ZAR's ~1.23 range
    f = fit(t, np.log(zar), np.log(sdr))
    dw = abs(f["w"] - W_TRUE) / f["se_w"]
    dc = abs(f["crawl_pct"] - CRAWL_TRUE) / f["se_crawl"]
    print(f"n={N:4d}  w={f['w']:.4f} +/-{f['se_w']:.4f} ({dw:.1f} sd)   "
          f"crawl={f['crawl_pct']:+.2f} +/-{f['se_crawl']:.2f} ({dc:.1f} sd)   "
          f"resid={f['rmse_bp']:.2f} bp  cond={f['cond']:.1f}")
