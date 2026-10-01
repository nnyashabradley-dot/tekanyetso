# Tekanyetso

Recovering the Bank of Botswana's crawling-peg rule from its own published
exchange rates, detecting when the Bank changes the rule's settings without
being told, and predicting the next published rate in public before it appears.

<!-- SCOREBOARD:START -->

## Live track record

*Updated 2026-10-01 09:00. Latest Bank publication: 2026-09-29.*

- **9 predictions scored.** Mean absolute error 2.50 bp, RMS 2.95 bp, mean -0.93 bp.
- For scale: the same model's out-of-sample RMS over 2011-2026 in the backtest is 3.56 bp, and the Bank's own four-decimal rounding puts a floor of about 3 bp under any prediction anchored on a published rate.
- **Coverage: 9 of 13 publication days predicted.** Missed: 2026-09-17, 2026-09-23, 2026-09-25, 2026-09-29.
- Scheduled predictions run at 09:00 Gaborone time; the Made column shows the actual time. A *same day* prediction was made when the Bank's file did not yet contain that day's rate.

| Target date | Made | Predicted | Error (bp) |
|---|---|---:|---:|
| 2026-09-30 | 2026-10-01 09:00 | -1.347494 | pending |
| 2026-09-28 | 2026-09-25 09:12 | -1.347959 | +3.34 |
| 2026-09-24 | 2026-09-24 09:09 (same day) | -1.355184 | -4.83 |
| 2026-09-22 | 2026-09-22 09:00 (same day) | -1.357726 | -5.06 |
| 2026-09-21 | 2026-09-19 09:00 | -1.357263 | +1.59 |
| 2026-09-18 | 2026-09-18 09:00 (same day) | -1.356616 | -0.17 |
| 2026-09-16 | 2026-09-15 10:20 | -1.354681 | -1.16 |
| 2026-09-15 | 2026-09-15 09:00 (same day) | -1.354924 | -2.78 |
| 2026-09-14 | 2026-09-12 09:05 | -1.353263 | +2.12 |
| 2026-09-11 | 2026-09-10 22:28 | -1.345392 | -1.41 |

Every prediction ever made is in `predictions.csv`, which is append-only; its git history shows when each row was written.

<!-- SCOREBOARD:END -->

## What this is

Since 2005 the Pula has been tied to a basket of the South African rand and the
IMF's SDR, and the peg is moved down a little every day at an announced annual
rate of crawl. The Bank announces the weights and the crawl rate in general
terms. It does not publish the daily arithmetic, the exact dates it changes its
settings, or how tightly it holds to its own formula.

The rule turns out to be

    w * log(ZAR per Pula) + (1 - w) * log(SDR per Pula)  =  a + b * t

and it holds to within the Bank's own four-decimal rounding. Everything else
follows from estimating `w` and `b` and finding the days they change.

## Results

Details, evidence and caveats are in [`FINDINGS.md`](FINDINGS.md) and
[`RESULTS.md`](RESULTS.md).

- **Every announced change of crawl rate since 2005 is found from the data
  alone**, 14 out of 14: eight within a week of the start of the announced month,
  thirteen within a month.
- **The basket was reweighted in January 2025 without an announcement**, six
  months before the July 2025 crawl change it had previously been attributed to.
- **The Bank steps the peg once per published day, not once per calendar day.**
  Over the whole post-2005 history a trading-day model beats a calendar-day model
  by 284 log-likelihood units with the same number of parameters.
- **PELT returns exactly the same breakpoints as the exhaustive dynamic program**
  with 10 to 30 times fewer cost evaluations, and its minimum-segment handling is
  shown to be necessary by a concrete counterexample.
- **In a backtest that provably never sees the future**, a Kalman filter predicts
  the next published rate with 2.68 bp RMS error, against 3.56 bp for the model
  that currently writes the live log.

## Run it

Python 3 and numpy. `figures.py` also needs matplotlib.

```
pip install numpy
python fetch.py                          # pull the Bank's table, validate, store
python estimate.py                       # rolling estimate of w and the crawl
python estimate.py --log                 # also append tomorrow's prediction
python segment.py --from 2005-06-01 --mult 10   # find the policy changes (PELT)
python kalman.py                         # time-varying w and crawl
python backtest.py --check-live          # no-lookahead backtest, all models
python bench_pelt.py                     # PELT against the exact DP
python scoreboard.py                     # refresh the track record above
python figures.py                        # figures for the writeup
```

Every analysis script has a `--self-test` or a built-in check:
`segment.py --self-test`, `kalman.py --self-test`, `diagnose.py --self-test`,
and `backtest.py` checks its own blindness on every run.

## Checking the live log

`predictions.csv` is append-only and committed by the scheduled daily run, so
its git history shows when each row was written. `backtest.py --check-live`
re-derives every logged prediction from the published data and reproduces each
one to eight decimal places, which shows the logged numbers came from the model
described here. The coverage line in the track record above shows any day on
which no prediction was made.

Git commit times are self-reported, so this is evidence that is hard to fake
quietly, not a cryptographic guarantee.

## Data

One table: `https://www.bankofbotswana.bw/content/exchange-rates`, CSV export
`/export/exchange-rates.csv`, quoted as foreign currency per Pula. The Bank
publishes its own SDR column alongside the rand, so no IMF or ECB data is
needed. `data/rates.csv` is append-only; `fetch.py` reports revisions rather than
overwriting. `data/crawl_changes.csv` is the Bank's table of announced crawl
changes, used only for scoring.

## Files

| File | Purpose |
|---|---|
| `fetch.py` | Download, screen and store the Bank's table |
| `estimate.py` | Rolling OLS estimate; writes the daily prediction |
| `segment.py` | Changepoint detection: exact DP and PELT |
| `bench_pelt.py` | PELT against the DP: identical answers, cost curves |
| `kalman.py` | Kalman filter with time-varying weight and crawl |
| `backtest.py` | No-lookahead backtest of every model, with blindness check |
| `diagnose.py` | Residual structure within a period |
| `scoreboard.py` | Scores the live log and writes the track record above |
| `figures.py` | Figures for the writeup |
| `FINDINGS.md`, `RESULTS.md` | What was found, with evidence |
