#!/usr/bin/env python3
"""
Tekanyetso -- figures for the writeup.

Needs matplotlib (the rest of the project needs only numpy):
    pip install matplotlib
    python3 bench_pelt.py        # writes bench_pelt.csv
    python3 backtest.py          # writes backtest.csv
    python3 kalman.py            # writes kalman.csv
    python3 segment.py --from 2005-06-01 --mult 10   # writes segments.csv
    python3 figures.py           # writes figures/*.png

Three figures:
  fig1_pelt_runtime.png   cost evaluations against n, DP and PELT
  fig2_weight_path.png    Kalman-filtered weight against the PELT segments
  fig3_backtest.png       out-of-sample error by model, against the floor
"""

import csv
import os
from datetime import date

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#e4e3df"
S1 = "#2a78d6"     # categorical slot 1, blue
S2 = "#eb6834"     # categorical slot 2, orange

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "axes.edgecolor": GRID, "axes.labelcolor": INK2, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2, "font.size": 10,
    "axes.titlesize": 12, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "legend.frameon": False, "savefig.dpi": 200, "savefig.bbox": "tight",
})


def rows(path):
    return list(csv.DictReader(open(path)))


def fig1():
    r = [x for x in rows("bench_pelt.csv") if x["kind"] == "synthetic" and x["mult"] == "1"]
    n = np.array([int(x["n"]) for x in r])
    dp = np.array([int(x["evals_dp"]) for x in r])
    pe = np.array([int(x["evals_pelt"]) for x in r])
    real = [x for x in rows("bench_pelt.csv") if x["kind"] == "real" and x["mult"] == "1"]

    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.plot(n, dp, color=S1, lw=2, marker="o", ms=5, label="Naive DP")
    ax.plot(n, pe, color=S2, lw=2, marker="o", ms=5, label="PELT")
    ax.set_yscale("log")
    ax.set_xlabel("Observations, n")
    ax.set_ylabel("Cost evaluations (log scale)")
    ax.set_title("PELT returns the same breaks for a fraction of the work")
    ax.text(n[-1] * 1.01, dp[-1], "  Naive DP", color=INK, va="center")
    ax.text(n[-1] * 1.01, pe[-1], "  PELT", color=INK, va="center")
    if real:
        x = real[0]
        ax.scatter([int(x["n"])] * 2, [int(x["evals_dp"]), int(x["evals_pelt"])],
                   s=60, facecolor=SURFACE, edgecolor=INK2, zorder=5, lw=1.2)
        ax.annotate("Bank of Botswana data,\n2005-2026 (hollow)", (int(x["n"]), int(x["evals_pelt"])),
                    xytext=(10, -28), textcoords="offset points", color=INK2, fontsize=9)
    ax.legend(loc="upper left")
    ax.set_xlim(0, n[-1] * 1.18)
    fig.text(0, -0.10, "Synthetic data with a policy change roughly every 300 observations, "
             "penalty 1x BIC, minimum segment 90.\nIdentical breakpoints and optimal cost "
             "in every run.", color=INK2, fontsize=8.5)
    fig.savefig("figures/fig1_pelt_runtime.png")
    plt.close(fig)


def fig2():
    k = rows("kalman.csv")
    kd = [date.fromisoformat(x["date"]) for x in k]
    kw = np.array([float(x["w"]) for x in k])
    seg = rows("segments.csv")
    fig, ax = plt.subplots(figsize=(8, 4.2))
    ax.plot(kd, kw, color=S1, lw=1.2, label="Kalman filter (daily)")
    first = True
    for s in seg:
        a, b = date.fromisoformat(s["start_date"]), date.fromisoformat(s["end_date"])
        ax.plot([a, b], [float(s["w"])] * 2, color=S2, lw=2.4, solid_capstyle="butt",
                label="PELT segments" if first else None)
        first = False
    ax.set_ylabel("Weight on the rand, w")
    ax.set_title("The basket weight, two ways")
    ax.set_ylim(0.40, 0.70)
    ax.xaxis.set_major_locator(mdates.YearLocator(3))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax.annotate("January 2025:\nunannounced", (date(2025, 1, 6), 0.50),
                xytext=(date(2019, 6, 1), 0.60), color=INK2, fontsize=9,
                arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8))
    ax.legend(loc="upper right")
    fig.text(0, -0.04, "Kalman process noise fitted on 2005-2010 only. PELT: penalty "
             "10x BIC, minimum segment 90.", color=INK2, fontsize=8.5)
    fig.savefig("figures/fig2_weight_path.png")
    plt.close(fig)


def fig3():
    b = rows("backtest.csv")
    names = {"live": "Live model (calendar-day crawl)",
             "trading": "Trading-day crawl",
             "trading-fit": "Trading-day, fitted anchor",
             "kalman": "Kalman filter",
             "no-crawl": "No crawl (baseline)"}
    sel = [x for x in b if date.fromisoformat(x["made_on"]) >= date(2011, 1, 1)]
    vals = {}
    for m in names:
        e = np.array([float(x[f"err_{m}_bp"]) for x in sel if x[f"err_{m}_bp"] != ""])
        vals[m] = np.sqrt((e ** 2).mean())
    order = sorted(vals, key=lambda m: vals[m])
    fig, ax = plt.subplots(figsize=(7, 3.6))
    y = np.arange(len(order))
    ax.set_axisbelow(True)
    ax.barh(y, [vals[m] for m in order], color=S1, height=0.55)
    for i, m in enumerate(order):
        ax.text(vals[m] + 0.05, i, f"{vals[m]:.2f}", va="center", color=INK, fontsize=9)
    ax.set_yticks(y, [names[m] for m in order])
    ax.invert_yaxis()
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("Out-of-sample RMS error, basis points (2011-2026)")
    ax.set_title("Next-day prediction error, no lookahead")
    ax.set_xlim(0, 4.4)
    fig.text(0, -0.13, "3,725 predictions per model. Rounding floor on the target day "
             "averages 2.18 bp; a prediction anchored on\ntoday's published value "
             "carries today's rounding as well, about 3.08 bp.", color=INK2, fontsize=8.5)
    fig.savefig("figures/fig3_backtest.png")
    plt.close(fig)


if __name__ == "__main__":
    os.makedirs("figures", exist_ok=True)
    fig1(); fig2(); fig3()
    print("wrote figures/fig1_pelt_runtime.png, fig2_weight_path.png, fig3_backtest.png")
