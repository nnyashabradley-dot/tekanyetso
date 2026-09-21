# Tekanyetso: Technical Results

**21 September 2026.** Everything built after `FINDINGS.md`: PELT, the Kalman
filter, the no-lookahead backtest and the live scoreboard. Every number here was
produced from `data/rates.csv` as fetched on 20 September (6,330 observations,
4 January 2001 to 18 September 2026) and can be regenerated with the command
shown beside it.

This is material for the writeup, not the writeup. The PELT correctness argument
is deliberately left out; section 1.4 says what it has to cover.

---

## 1. PELT

### 1.1 What was built

`segment_pelt` in `segment.py`. Same recursion and same cost object as
`segment_dp`; the only difference is the array of candidate starts passed to
`cost.cost_many` at each step. PELT is now the default in `segment.py`, and the
DP stays available as `--method dp` as the reference.

Pruning rule: after computing `F[tau]`, a start `s` that was evaluated at `tau`
is marked dead if `F[s] + C(s, tau) > F[tau] + PRUNE_TOL`. A dead start is not
removed immediately. It remains a candidate until `tau + min_seg` and is dropped
after that.

`PRUNE_TOL = 1e-5` cost units. The cumulative-sum cost is accurate to about
1e-6, so the pruning inequality can be violated by rounding alone; the tolerance
only makes pruning more conservative.

### 1.2 Identical answers

| Test | Result | Command |
|---|---|---|
| PELT against exhaustive search, small series | identical | `segment.py --self-test` (test 4) |
| PELT against DP, 48 synthetic configurations | 0 mismatches | `segment.py --self-test` (test 4) |
| PELT against DP, 360 synthetic configurations, min_seg 20/60/90, penalty 0.5x to 20x BIC | 0 mismatches | development run |
| PELT against DP, 16,987 adversarial random cases | 0 mismatches | development run |
| PELT against DP, real data 2005-2026, 1x and 10x BIC | identical breaks and cost | `bench_pelt.py` |

"Identical" means the same breakpoints and the same optimal cost to 1e-6.

### 1.3 The minimum-segment delay is necessary

Pruning without the delay, removing a dead start at once as in the textbook
algorithm, lost the true optimum in **546 of 16,987 adversarial cases (3.2%)**.
The delayed version was exact in all of them.

Smallest reproducible case (`bench_pelt.py --counterexample`):

```
case 40: n=63, min_seg=11, penalty=1.189
  exact DP            cuts [19, 32, 48]  cost 60.8531
  PELT, no delay      cuts [19, 36, 48]  cost 60.9792   <-- worse
  PELT, delay min_seg cuts [19, 32, 48]  cost 60.8531
```

On the peg data and the synthetic peg data, the undelayed version happened to
give the right answer in every case tried; the failures appeared only in the
adversarial cases, which have short segments and low penalties. So the bug did
not show up on the real data, and only the proof or the adversarial search
reveals it.

### 1.4 What the correctness argument has to establish

Left for you, per the standing decision that writing it is how you find out
whether you understand it. The argument needs to answer:

1. Why `C(s,t) + C(t,T) <= C(s,T)` holds for this cost, so that `K = 0`.
2. For which later `T` a start pruned at `tau` can be shown never to be optimal,
   and which step of the argument decides that.
3. Why the removal therefore has to wait, and why `min_seg` is exactly enough.
   The counterexample above is what it looks like when it does not wait.
4. Why `PRUNE_TOL` cannot cost exactness.

### 1.5 Cost

Cost evaluations, the primary metric, since both methods share one vectorised
cost function. From `bench_pelt.py`; figure: `figures/fig1_pelt_runtime.png`.

Synthetic, a policy change every ~300 observations, min_seg 90, 1x BIC:

| n | breaks | DP evaluations | PELT evaluations | ratio | PELT per observation |
|---|---|---|---|---|---|
| 500 | 0 | 52,092 | 52,092 | 1.0 | 104 |
| 999 | 2 | 337,520 | 97,700 | 3.5 | 98 |
| 1,998 | 5 | 1,657,199 | 231,457 | 7.2 | 116 |
| 3,991 | 12 | 7,271,480 | 483,912 | 15.0 | 121 |
| 6,000 | 19 | 16,950,842 | 724,340 | 23.4 | 121 |
| 7,982 | 25 | 30,455,199 | 1,036,623 | 29.4 | 130 |

DP grows as n squared over 2. PELT's cost per observation stays roughly flat at
100 to 130, which is the near-linear behaviour. At 10x BIC the ratio at n = 7,982
is 21.5.

Real data, 1 June 2005 to 18 September 2026, n = 5,254:

| Penalty | breaks | DP | PELT | ratio | wall clock DP / PELT |
|---|---|---|---|---|---|
| 1x BIC | 22 | 12,885,515 | 946,985 | 13.6 | 1.5 s / 0.4 s |
| 10x BIC | 16 | 12,885,515 | 1,302,545 | 9.9 | 1.5 s / 0.4 s |

Why PELT's cost per observation sits near 100 rather than near the number of
live candidates: the delay keeps every pruned start alive for another `min_seg`
= 90 steps. The minimum segment length therefore costs roughly a constant factor
of `min_seg`. That is the price of correctness, and it is worth one sentence in
the writeup.

Wall clock understates the difference because the cost function is already O(1)
and vectorised.

### 1.6 Segments, now from PELT

`segment.py --from 2005-06-01 --mult 10` gives the same 16 breaks and the same
table as `FINDINGS.md` section 1, with three more days of data.

**Correction to `FINDINGS.md` section 1**, now fixed in that file: measured from
the first day of the announced month, 8 of the 14 detected dates are within a
week and 13 within a month. April 2010 is the outlier, detected on 12 May, 41
days late. The earlier text said eleven were within about a week.

| Announced | Detected | Days |
|---|---|---|
| 2006-07 | 2006-06-30 | -1 |
| 2007-07 | 2007-07-03 | +2 |
| 2009-03 | 2009-02-27 | -2 |
| 2010-04 | 2010-05-12 | +41 |
| 2012-06 | 2012-06-04 | +3 |
| 2015-01 | 2015-01-28 | +27 |
| 2016-01 | 2015-12-14 | -18 |
| 2017-01 | 2017-01-04 | +3 |
| 2018-01 | 2017-12-29 | -3 |
| 2019-01 | 2018-12-27 | -5 |
| 2020-01 | 2019-12-20 | -12 |
| 2020-05 | 2020-05-06 | +5 |
| 2023-01 | 2022-12-21 | -11 |
| 2025-07 | 2025-07-30 | +29 |

---

## 2. Kalman filter

### 2.1 Model

`kalman.py`. State `[m, b, w]`: the parity (log SDR per Pula implied at the
previous day's rand/SDR spread), the crawl per step, and the weight. All three
follow random walks. The observation noise is the rounding floor, per
observation and at the current weight, so it is known rather than estimated.
Only the three process variances are fitted, by maximum likelihood, on
**1 June 2005 to 31 December 2010 only**, then frozen.

Design choice worth explaining in the writeup: the state holds the parity `m`
rather than the basket index. A reweighting changes the index by about
`0.05 x 3.1`, roughly 15 percent, because the index changes definition, while
the Pula itself does not move. With the index in the state, every reweighting
would need a 1,500 bp jump exactly correlated with the jump in `w`. With `m`, a
reweighting leaves `m` continuous and only `w` moves.

Covariance update in Joseph form. Optimiser: a small Nelder-Mead written in
numpy, so the project still needs only numpy.

### 2.2 Fitted values (2005-2010)

| Process noise | Fitted standard deviation per published day |
|---|---|
| Parity `m` | 0.107 bp |
| Weight `w` | 0.00128 |
| Crawl `b` | 2.76e-06 |

The parity noise says the implied rate wanders by about 0.1 bp a day beyond the
crawl, on top of the known rounding. It absorbs both any looseness in how the
Bank applies its formula and any imperfection in this model, so it is an upper
bound on the Bank's looseness, not a measurement of it. It is also a daily
increment rather than a deviation, so do not compare it directly with the 1.1 bp
bound in `STATUS.md` section 3.

### 2.3 Checks

| Check | Result | Command |
|---|---|---|
| Recovers known w and crawl across a synthetic reweighting | 0.658 / 0.598 vs 0.65 / 0.60; -4.99 / -2.37 vs -4.80 / -2.30 | `kalman.py --self-test` |
| Joseph form keeps covariance positive definite | smallest eigenvalue 9.7e-03 | `kalman.py --self-test` |
| Causal: corrupting later data leaves earlier predictions unchanged | bit-identical | `kalman.py --self-test` |
| Calibration on real data, 2005-2026 | standardised one-step RMS 1.05 (1.00 is ideal) | `kalman.py` |
| One-step RMS error on real data | 2.41 bp | `kalman.py` |

### 2.4 What it shows

`figures/fig2_weight_path.png`. The filtered weight follows the PELT segments
closely, including the unannounced January 2025 reweighting. It lags each break
by a few weeks, which is what a random walk does with a jump. The filtered crawl
is noisier: it wanders by a few tenths of a percent inside regimes. That
matches README finding 3, that the crawl is invisible day to day and only
emerges over long windows.

### 2.5 The pooled test of trading-day stepping

`FINDINGS.md` section 3 left open whether the trading-day mechanism holds over
the whole history. `kalman.py --compare-dt` fits the same model twice on all
post-2005 data, once with one crawl step per published day and once with one per
calendar day. Same data, same number of parameters.

| Clock | Log-likelihood | Fitted parity noise |
|---|---|---|
| Trading day | 36,114.1 | 0.231 bp |
| Calendar day | 35,829.8 | 0.443 bp |

**Trading-day wins by 284.3 log-likelihood units.** The calendar model also needs
twice the parity noise to absorb the weekend misfit.

By era (development script, same fitted values):

| Era | Observations | Gain for trading day |
|---|---|---|
| 2005-2007 | 552 | +167.9 |
| 2008-2012 | 1,236 | +78.9 |
| 2013-2019 | 1,722 | -0.6 |
| 2020-2026 | 1,654 | +38.1 |

This is the pattern the mechanism predicts. From 2013 to 2019 the crawl was
between -0.30 and +0.38 percent a year, so there was almost no crawl to step and
the two clocks cannot be told apart. Where the crawl is large the trading-day
model wins, including decisively in the modern era. This answers the question
behind open question 3 in `FINDINGS.md` by a different and stronger route than
the per-segment weekday test proposed there.

---

## 3. No-lookahead backtest

### 3.1 Design

`backtest.py`. At each published day, every model predicts the log basket index
for the next weekday using only data published on or before that day, exactly as
the live log does, and is scored against what the Bank then published, using the
model's own weight.

Scoring with the model's own weight is the same as predicting the published SDR
rate given the rand/SDR spread on the target day, because `I = ls + w*u` and the
`w*u` term cancels. So all models are scored on one target, and none of them
forecasts the currency market. The claim being tested is about the Bank.

Models:

| Model | What it is |
|---|---|
| live | rolling-ols-v0: 90-observation OLS, calendar-day crawl, anchored on today's published index. The model that writes `predictions.csv`. |
| trading | same fit, one crawl step per published day |
| trading-fit | same, anchored on the regression's fitted value instead of today's published value |
| kalman | section 2, hyperparameters from 2005-2010 only |
| no-crawl | today's index unchanged; the baseline |

### 3.2 Proof of blindness

Checked on every run. At three cutoffs (15 March 2011, 1 August 2017,
31 December 2024) the backtest re-runs every model on the history deleted after
the cutoff, and again with every value after the cutoff replaced by random
numbers. Every prediction made on or before the cutoff must be bit-identical to
the full run. **All five models pass at all three cutoffs.**

The one place future data could leak is the Kalman hyperparameters, fitted on
2005-2010. The headline table therefore starts on 1 January 2011, where they are
out of sample.

### 3.3 The live system is the backtested system

`backtest.py --check-live` re-derives all six logged predictions from the
published data. Each matches the logged value to within 4e-9, which is the
eight-decimal rounding in the log. The numbers in `predictions.csv` came from
exactly this model.

### 3.4 Headline: out of sample, predictions made 2011-01-01 to 2026-09-18

3,725 scored predictions per model. Rounding floor on the target day averages
2.18 bp; a prediction anchored on today's published value carries today's
rounding too, so about 3.08 bp is the floor for those models.
Figure: `figures/fig3_backtest.png`.

| Model | Mean error | RMS | Mean abs | 95th pct | Monday mean | Tue-Fri mean |
|---|---|---|---|---|---|---|
| live | -0.03 | 3.56 | 2.78 | 6.82 | +0.69 | -0.20 |
| trading | +0.01 | 3.54 | 2.76 | 6.76 | +0.20 | -0.04 |
| trading-fit | +0.05 | 3.75 | 2.62 | 6.60 | +0.11 | +0.04 |
| kalman | -0.02 | **2.68** | **2.17** | **4.79** | +0.06 | -0.04 |
| no-crawl | -0.48 | 3.61 | 2.83 | 6.89 | -0.27 | -0.53 |

All in basis points.

Diebold-Mariano tests on squared error, Newey-West variance (t above 2 means
the second model is significantly more accurate):

| Comparison | t |
|---|---|
| live vs trading | +2.83 |
| live vs kalman | +15.33 |
| trading vs kalman | +15.19 |
| no-crawl vs live | +5.40 |

### 3.5 Early period, 2005-10-01 to 2010-12-31

Kalman hyperparameters are in sample here, so treat its row as optimistic.

| Model | Mean error | RMS | Monday mean | Tue-Fri mean |
|---|---|---|---|---|
| live | -0.03 | 1.91 | +1.26 | -0.34 |
| trading | +0.04 | 1.77 | -0.08 | +0.07 |
| trading-fit | +0.58 | 2.93 | +0.53 | +0.59 |
| kalman | +0.01 | 1.37 | -0.06 | +0.03 |
| no-crawl | -1.26 | 2.22 | -1.38 | -1.23 |

### 3.6 What the backtest says

1. **The live model has a weekend bias, and it is the trading-day mechanism.**
   It applies three days of crawl over a weekend where the Bank applies one, so
   its Monday predictions are too low: +0.69 bp on average from 2011, +1.26 bp
   before. The trading-day model removes it. This is `FINDINGS.md` section 3
   showing up as a forecasting error.
2. **The Kalman filter is 25 percent more accurate than the live model**, 2.68
   against 3.56 bp RMS, and its worst days are much better (95th percentile 4.79
   against 6.82). The gain comes mainly from where the prediction starts. Any
   model that starts from today's published value carries today's rounding as
   well as tomorrow's, so it cannot get below about 3.08 bp however good its
   crawl and weight are. The filter starts from its estimate of the true parity,
   and its 2.68 bp is below that limit, closer to the 2.18 bp floor that every
   prediction faces.
3. **The crawl is detectable in the forecasts.** Ignoring it leaves a mean error
   of -0.48 bp from 2011 and -1.26 bp before, and loses to the live model with
   t = 5.4.
4. **Anchoring on the fitted value alone is worse** (what did not work). The
   90-observation regression lags every policy change, so its fitted value is
   wrong for weeks afterwards. Its RMS within 20 publications of a detected
   break is 5.20 bp against 3.64 elsewhere, and its early-period bias is
   +0.58 bp. The Kalman filter gets the benefit of a fitted anchor without the
   lag because it updates every day.

Errors within 20 publications after each detected policy change, from 2011
(break dates are used only to sort errors for this table, never by a model):

| Model | RMS near a break | RMS elsewhere | Worst near a break |
|---|---|---|---|
| live | 4.06 | 3.53 | 12.5 |
| trading | 4.07 | 3.51 | 12.5 |
| trading-fit | 5.20 | 3.64 | 19.1 |
| kalman | 3.67 | 2.60 | 12.0 |
| no-crawl | 4.10 | 3.58 | 12.5 |

---

## 4. Live scoreboard

`scoreboard.py`, now called by `run_daily.bat` after the prediction, writes a
track record into `README.md` between two marker lines, and the daily commit now
includes `README.md`.

As of 21 September: 5 predictions scored, mean absolute error 1.53 bp, RMS
1.76 bp. **Coverage 5 of 6 publication days; 17 September missed.** The coverage
line is what makes a skipped day visible.

---

## 5. Caveats to carry into the writeup

- The backtest tests the rule claim: given the rand/SDR spread, what will the
  Bank publish. It is not a forecast of the exchange rate.
- The Kalman comparison is out of sample from 2011, but its model structure was
  chosen after seeing the data. Only its fitted values are clean.
- The trading-day evidence carries no information for 2013 to 2019, when the
  crawl was near zero.
- The live log has 5 scored predictions. It becomes evidence of anything only
  over months. The backtest is what supports the accuracy claims until then.
- The 2005-2006 first-year anomaly (`FINDINGS.md` section 4) and the January
  2025 reweighting still need the documentary checks listed there.

## 6. Decisions left for you

1. **Whether to change the live model.** The Kalman filter is clearly better. I
   have not switched it, because the live log is a public record and changing
   what it measures is your call. The cleanest option is to log a second model
   alongside rolling-ols-v0 rather than replacing it; `predictions.csv` already
   has a `model` column, but `log_prediction` currently refuses a second row for
   the same date and would need to key on date and model.
2. **The scheduler.** One missed day in the first six.
3. **OpenTimestamps**, if you want the log's timestamps anchored independently of
   git. Not built.

## 7. Files added or changed

| File | Change |
|---|---|
| `segment.py` | PELT added and made the default; `--method`; self-test 4; refine counts evaluations correctly |
| `bench_pelt.py` | new: PELT against DP, cost curves, counterexample |
| `kalman.py` | new |
| `backtest.py` | new |
| `scoreboard.py` | new |
| `figures.py` | new, needs matplotlib |
| `figures/*.png` | three figures |
| `estimate.py` | docstring only: the old calendar-day reasoning corrected |
| `run_daily.bat` | runs the scoreboard and commits README.md |
| `README.md` | rewritten for a stranger; live track record at the top |
| `FINDINGS.md` | section 1 count corrected |
| `.gitignore` | derived CSVs from the new scripts |
