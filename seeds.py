import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from estimate import fit
W, C = 0.50, -2.76
for N in (25, 60, 90, 180, 250):
    sw, sc, hits_w, hits_c = [], [], 0, 0
    for seed in range(200):
        rng = np.random.default_rng(seed)
        t = np.arange(N, dtype=float) * 1.4
        u = np.cumsum(rng.normal(0, 0.006, N))
        ls = np.log(0.0557) + (C/100/365)*t - W*u
        sdr = np.round(np.exp(ls), 4); zar = np.round(np.exp(ls+u)*22.1, 4)
        f = fit(t, np.log(zar), np.log(sdr))
        sw.append(f["se_w"]); sc.append(f["se_crawl"])
        hits_w += abs(f["w"]-W) <= 1.96*f["se_w"]
        hits_c += abs(f["crawl_pct"]-C) <= 1.96*f["se_crawl"]
    print(f"n={N:4d}  median se(w)={np.median(sw):.4f}  median se(crawl)={np.median(sc):.3f} %/yr"
          f"   95% CI coverage: w {hits_w/2:.0f}%  crawl {hits_c/2:.0f}%")
