# Tekanyetso: Progress and Handover

**15 September 2026. Week 3 complete. The system is live and logging.**

This document is the single entry point. It assumes no memory of any previous
conversation. It replaces the earlier version of `PROGRESS.md`, which was written
in week 1 before most of the data had been looked at. The old version is in git
history if it is ever needed.

---

## 1. What the project is

The Botswana Pula does not float and is not pegged to one currency. Since May
2005 the Bank of Botswana has run a crawling band. The Pula is tied to a basket,
currently 50 percent South African rand and 50 percent IMF SDR, and that peg is
adjusted downward gradually, currently at 2.76 percent per year, to keep the real
effective exchange rate stable.

The Bank announces the basket weights and the crawl rate at a high level. It does
not publish the exact daily arithmetic, the precise dates it changes its settings,
or how tightly it holds the published rate to its own formula.

The job is to recover all of that from published numbers alone, detect the dates
the parameters changed without being told, forecast the next published rate before
it appears, and log every forecast publicly so the claim cannot be adjusted after
the fact.

### Why it is a computer science project

The centrepiece is changepoint detection. Two versions of the same algorithm get
built: a naive dynamic program that is exact but quadratic, and PELT, which prunes
candidate breakpoints that provably cannot be optimal and returns the identical
answer in near-linear time. The correctness argument for the pruning, plus the
runtime curves diverging, is the technical core of the writeup.

Unusually for this kind of problem, there is real ground truth to score against.
The Bank publishes a table of every rate-of-crawl change since 2005.

## 2. Where it stands

| Component | State |
|---|---|
| Data pipeline | Done. Runs daily against the live endpoint. |
| Estimator | Done. Recovers weight and crawl across the full history. |
| Changepoint detector | Naive DP done and scored. PELT not started. |
| Forecaster and log | Live since 10 September. Publishes to GitHub daily. |

Repository: `github.com/nnyashabradley-dot/tekanyetso`, public.

Local folder: `Documents\Developing stuff\Pula peg` on Windows, Python via the
`py` launcher. The only dependency is numpy. Everything else is standard library.

### What runs daily

`run_daily.bat`, scheduled at 09:00 through Task Scheduler. It fetches, estimates,
appends one prediction, commits and pushes. Output goes to `run.log`.

The 09:00 time is deliberate and should never be moved. Evidence from
`predictions.csv` shows the Bank published on 15 September somewhere between 09:00
and 10:20, so a 09:00 run predicts a rate that does not yet exist anywhere. That
is the strongest form of the claim.

Two runs were missed last week, on 11 and 12 September. Task Scheduler history has
not been checked. Do not enable "run task as soon as possible after a scheduled
start is missed", because a run firing at 14:00 would predict a rate the Bank
published hours earlier and the row would look identical to a legitimate one. A
missed day should stay a visible gap.

Maintenance is one weekly check: `predictions.csv` should gain about five rows a
week. If it did not, read the tail of `run.log`.

## 3. Document map

Read in this order. Where documents conflict, the later one wins.

| Document | What it holds |
|---|---|
| `PROGRESS.md` | This file. Current state and what is left. |
| `FINDINGS.md` | Results from week 3, with evidence and caveats. |
| `STATUS.md` | Week 1 and 2 results. Sections 3, 4 and 6 are partly superseded by `FINDINGS.md`. |
| `WEEK1_ADDENDUM.md` | Data pipeline detail. |
| `docs/method-changepoint.md` | Method note for the segmentation algorithm. |
| `README.md` | Public-facing summary and how to run. |
| `tekanyetso-build-notes.md`, `tekanyetso-plan.md` | Original plans. Largely superseded, kept for the schedule. |

## 4. The model

The peg holds a weighted geometric basket on a straight line in log space:

    w * log(ZAR_per_Pula) + (1 - w) * log(SDR_per_Pula) = a + b * t

Geometric rather than arithmetic because an arithmetic average changes value if
you quote the rates the other way round, so it is not measuring anything real.

Writing `u = log(ZAR_pP) - log(SDR_pP)`, this rearranges to ordinary least squares
with one regressor of interest:

    log(SDR_pP) = a + b * t - w * u + e

so `w` is minus the coefficient on `u`, and the annual crawl is `365 * b` when `t`
is in calendar days. Sum-to-one is imposed algebraically rather than estimated,
which removes the degree of freedom that would otherwise make the weight and the
crawl collinear.

The single regressor is exogenous. `u = log(ZAR / SDR)` has the Pula cancelled out
of it, so it is a pure world-market quantity with no Botswana content. This is why
identification is clean rather than lucky, and it belongs in the writeup.

Note that `estimate.py` uses `t` in calendar days and annualises by 365, while
`segment.py` uses `t` in years directly. Both are correct. Do not mix them.

## 5. Facts worth not re-deriving

**Data source.** One file. `https://www.bankofbotswana.bw/content/exchange-rates`,
CSV export at `/export/exchange-rates.csv?page&_format=csv`. It carries CHN, EUR,
GBP, USD, ZAR, SDR and YEN, all quoted as foreign currency per Pula. ZAR is about
1.23, SDR about 0.056. The IMF and ECB are not needed, because the Bank publishes
its own SDR column alongside its rand column on the same dates.

Current coverage: 6,327 observations, 4 January 2001 to 15 September 2026.

**Column order differs between the CSV export and the HTML table.** Read by name.
Reading by position silently swaps the two columns the model uses.

**The measurement floor is not constant.** The Bank prints to four decimals, and
in log terms the resulting noise scales as one over the SDR level. It was 1.02
basis points in 2001 and is 2.61 today. It also depends on the weight. Any cost
function that assumes a constant floor will over-detect in recent data and
under-detect in early data.

**The crawl is invisible day to day.** At 2.76 percent per year it moves 0.76
basis points per day, well under the rounding noise. It only emerges once
cumulative drift beats the rounding, which takes tens of days. Use long windows
for the crawl. Treat short-window crawl estimates as unreliable.

**The prediction target is the rule, not the level.** Once the formula is known,
tomorrow's Pula rate is fully determined by tomorrow's rand and SDR, so predicting
it reduces to predicting world currency markets. The log currently claims
something different and testable: tomorrow's log basket index equals today's plus
one day of crawl. That is a claim about the Bank, needs no market forecast, and is
checkable from the Bank's own file.

## 6. Findings so far

Full evidence and caveats are in `FINDINGS.md`. In brief:

1. **The detector finds every announced crawl change.** Fourteen out of fourteen
   in the post-2005 window, eleven of them within about a week of the announced
   month.
2. **The basket was reweighted in January 2025 and the Bank did not announce it.**
   `STATUS.md` had placed this six months late and bundled it with the July 2025
   crawl change. They are two separate events.
3. **The crawl is applied on trading days, not calendar days.** The predicted
   weekday residual ramp matches observation with chi-squared 0.41 on 4 degrees of
   freedom, and the residual autocorrelation has its fundamental at lag 5. This
   does not bias the annual rate estimates, because a calendar-day fit averages
   the steps correctly.
4. **The framework's first year ran about half a percentage point fast.** From
   July 2006 onward the recovered rate matches the announced rate to within a few
   hundredths and stays that way for twenty years. Three mundane explanations have
   not been ruled out, and the next step there is documentary.
5. **The results do not depend on the noise weighting.** Every segment from July
   2007 onward is identical with and without iterative reweighting.

## 7. What is left

The schedule from `tekanyetso-plan.md`, updated for where things actually are.

### Week 4: PELT

This is the next thing to build and the most important remaining piece.

Three parts:

- **The pruning rule.** Discard a candidate start `s` when
  `F[s] + C(s, tau) + K > F[tau]`, where `K` is a constant depending on the cost
  function. For this cost `K = 0`, because splitting a segment can only reduce the
  residual sum of squares. The minimum segment length adds a wrinkle: candidates
  can only be pruned once they are older than that minimum.
- **The correctness argument.** An induction showing a pruned candidate cannot be
  optimal at any later `tau`. This is where the additivity of the cost function
  earns its keep. Writing this argument is how you find out whether you understand
  it, and it cannot be faked.
- **The empirical side.** Identical breaks to the naive DP on real and synthetic
  data, plus runtime curves. The DP baseline is already measured at 0.208, 0.338,
  0.415 and 0.457 cost evaluations per n squared as n grows through 500, 1000,
  2000 and 4000, converging on one half. PELT should flatten toward a constant
  multiple of n.

`segment.py` was written for this. Both algorithms reduce to "for each tau,
evaluate the cost over an array of candidate starts", and only the array differs.
The existing `--self-test` already compares the DP against exhaustive search, so
the same harness verifies PELT.

Keep cost evaluations as the primary metric rather than wall clock, since a
vectorised DP against a scalar-loop PELT would produce a meaningless speedup.

Scoring against announcements is largely done already, in `FINDINGS.md` section 1.

### Week 5: Kalman filter

Time-varying weight and crawl as hidden states, fitted by maximum likelihood.
Joseph form or a square-root filter if time allows, plain clipping if not. This is
a sophistication upgrade, not the spine of the project, and is the first thing to
cut if the schedule slips.

### Week 6: backtest, site, README

A backtest provably blind to future data. A live scoreboard. A repository a
stranger can run.

### Week 7: writeup

Eight to twelve pages including what did not work.

### Not yet started, from earlier plans

- **OpenTimestamps anchoring** of the prediction log. Git commit timestamps are
  self-reported and can be backdated, and GitHub's public push-event record only
  retains about ninety days. Anchoring each day's file hash to Bitcoin closes that
  permanently. About five minutes of setup, and it is what makes the log's central
  claim defensible rather than merely stated.
- **A scoreboard in the README**: last ten predictions against what the Bank
  published, running mean absolute error, and predictions logged versus trading
  days elapsed. That last number is what makes selective omission visible, which
  is the one weakness timestamping cannot address.
- **`docs/method.md`** beyond the changepoint section.

## 8. Immediate next actions

1. Check Task Scheduler history for the two missed runs on 11 and 12 September.
2. Add a pointer at the top of `STATUS.md` noting that sections 3, 4 and 6 are
   superseded by `FINDINGS.md`.
3. Start PELT.

Cheap and worth doing when convenient:

4. Check the Bank's press releases for late 2024 and early 2025 for any
   announcement of a basket reweighting.
5. Check the 2005 and 2006 Annual Reports and Monetary Policy Statements for the
   operational rate of crawl.
6. Review the 40 quarantined rows, all between 2001 and 2004. Some are clear
   errors. Others sit only 3 to 5 percent from their local median, which is inside
   a normal week of currency movement, so real observations may be being
   discarded.
7. Document what `seeds.py`, `synth.py` and `fixture_bob.csv` do. They are in the
   repository and referenced nowhere.

## 9. Standing decisions

- **Never cut the daily prediction log.** No lookahead. Never overwrite, never
  backfill silently.
- **Never move the 09:00 run time.**
- **Documentation happens alongside the code, not after it.** The point is
  understanding the work, not possessing it.
- **The writeup is mandatory.** Writing the PELT correctness argument is how you
  find out whether you understand it.
- **No vanity publication.** The project's credibility rests on a public error log
  that is hard to fake. Do not dilute it with a credential that is easy to fake.
- **The real paper, if it happens, is the comparative extension**: the same
  recovery across several African basket pegs, testing whether drift from the
  stated rule predicts reserve loss or devaluation. That is twelve to eighteen
  months out and needs the track record the log accumulates by itself.

## 10. Files

| File | State |
|---|---|
| `fetch.py` | Runs daily against the live endpoint. |
| `estimate.py` | Full history. Logs predictions with lead time recorded. |
| `segment.py` | Naive DP changepoint detector. Four self-tests, all passing. |
| `diagnose.py` | Residual structure diagnostic. Self-test passing. |
| `run_daily.bat` | Fetch, estimate, log, commit, push. Scheduled 09:00. |
| `data/rates.csv` | 6,327 rows, append-only, with `first_seen`. |
| `data/crawl_changes.csv` | Ground truth. 15 announced changes, month precision. |
| `data/quarantine.csv` | Screened-out rows with reasons. 40 rows, all 2001 to 2004. |
| `predictions.csv` | Live and public. |
| `rolling.csv`, `segments.csv` | Derived output, regenerated each run, gitignored. |
| `raw/` | Every download verbatim. Local only, gitignored. |
| `seeds.py`, `synth.py`, `fixture_bob.csv` | Present, undocumented. |

### How to run

```
py fetch.py                                  # pull, screen, store
py estimate.py --window 90 --log             # estimate, predict, append
py segment.py --self-test                    # verify the detector
py segment.py --from 2005-06-01 --mult 10    # segment the history
py diagnose.py --eras                        # residual structure by era
```
