# Tekanyetso: Findings

**15 September 2026.** Results from the first run of the changepoint detector and
the residual diagnostic against the full Bank of Botswana history.

This document records what was found, what the evidence is, and what has not been
ruled out. Where it contradicts `STATUS.md`, `PROGRESS.md` or the two original
planning documents, this one is later and should be taken as correct.

Everything here comes from `segment.py` and `diagnose.py` run on
`data/rates.csv` (6,327 observations, 4 January 2001 to 15 September 2026).

---

## 1. The detector finds every announced crawl change

Naive dynamic programming, BIC penalty at 10x, minimum segment 90 observations,
restricted to 1 June 2005 onward (5,251 observations). Sixteen breaks found in
1.2 seconds over 12.87 million cost evaluations.

Fourteen announced rate-of-crawl changes fall inside that window. The detector
found a boundary for all fourteen, and the recovered rate in each segment matches
the announced rate closely.

| Announced | Detected | Recovered | Announced rate |
|---|---|---|---|
| July 2006 | 2006-06-30 | -3.92 | -3.90 |
| July 2007 | 2007-07-03 | -2.26 | -2.30 |
| March 2009 | 2009-02-27 | -2.92 | -2.91 |
| April 2010 | 2010-05-12 | -2.58 | -2.61 |
| June 2012 | 2012-06-04 | -0.16 | -0.16 |
| January 2015 | 2015-01-28 | -0.01 | 0.00 |
| January 2016 | 2015-12-14 | 0.38 | 0.38 |
| January 2017 | 2017-01-04 | 0.25 | 0.26 |
| January 2018 | 2017-12-29 | -0.31 | -0.30 |
| January 2019 | 2018-12-27 | 0.29 | 0.30 |
| January 2020 | 2019-12-20 | -1.32 | -1.51 |
| May 2020 | 2020-05-06 | -2.86 | -2.87 |
| January 2023 | 2022-12-21 | -1.52 | -1.51 |
| July 2025 | 2025-07-30 | -2.77 | -2.76 |

Announced dates are month precision, so a detected date within a month of the
announcement counts as a match. Eleven of the fourteen land within about a week.

Recall is therefore 14 out of 14. Precision is harder to state, because two of the
sixteen breaks are not announced crawl changes. One of them is a real finding
(section 2) and one is an artefact of model misfit (section 4).

The January 2020 segment is the weakest match. The recovered -1.32 sits 0.19
percentage points from the announced -1.51, and the segment is only 90
observations long, which is the minimum allowed. It is bounded by the January
2020 and May 2020 announcements, which are close together.

## 2. The basket was reweighted in January 2025, and the Bank did not announce it

`STATUS.md` section 4 reported five basket weight changes recovered from rolling
windows, bracketed to 90 days each, and concluded that all five coincided with an
announced crawl change. It inferred from this that the Bank adjusts weight and
crawl at the same review and announces only the crawl.

The dynamic program dates the weight changes directly rather than bracketing
them. Four of the five do coincide with an announced crawl change. The fifth does
not.

| Weight change | Detected | Crawl before | Crawl after | Announced crawl change |
|---|---|---|---|---|
| 0.65 to 0.60 | 2007-07-03 | -3.92 | -2.26 | July 2007 |
| 0.60 to 0.55 | 2012-06-04 | -2.58 | -0.16 | June 2012 |
| 0.55 to 0.50 | 2015-01-28 | -0.16 | -0.01 | January 2015 |
| 0.50 to 0.45 | 2017-01-04 | 0.38 | 0.25 | January 2017 |
| **0.45 to 0.50** | **2025-01-06** | **-1.52** | **-1.55** | **none** |

The January 2025 change moves the weight on the rand from 0.450 to 0.502 while
the crawl stays where it was. The announced rate of -1.51 percent ran from
January 2023 to June 2025 with no change recorded in between.

`STATUS.md` had placed this change in the window 8 July to 6 October 2025 and
matched it to the July 2025 crawl change. That is six months late. The July 2025
break is real and separate: the detector finds it at 2025-07-30, with the weight
already at 0.500 and staying there while the crawl moves from -1.55 to -2.77.

So these are two distinct events, and the Bank announced only the second one.

This also explains the observation in `STATUS.md` section 6 that the 2025
transition looked messier than the other four, wandering through 0.49, 0.48 and
0.51 before settling. A 90 day rolling window spanning two separate policy
changes six months apart will produce exactly that.

**Independent confirmation.** Fitting a single line across 2 January 2023 to
30 June 2025, which spans the January 2025 change, gives chi2/dof of 3.23 and an
unexplained residual of 4.33 basis points. Stopping the same window at
31 December 2024, just before the detected change, drops the excess below 1 basis
point and flattens the residual autocorrelation to near zero at every lag. The
misfit is entirely attributable to the break.

**What has not been done.** The press releases for late 2024 and early 2025 have
not been checked. If one of them announces a reweighting, then the finding is that
the crawl-change table is incomplete rather than that the Bank was silent. Either
way the date recovered from data stands.

## 3. The crawl is applied on trading days, not calendar days

`PROGRESS.md` section 4.4 argues that `t` must be in calendar days because the
crawl represents an inflation differential accumulating in real time, and weekends
still happen. The practical advice is correct. The reasoning is not.

If the Bank steps the parity once per trading day while the model assumes a
continuous calendar-day drift, the residual acquires a predictable weekly shape.
Over a weekend the model adds three days of drift and the Bank adds one, so the
mismatch builds up from Monday to Friday and resets. The size of the ramp is fixed
entirely by the fitted crawl rate, with nothing free to tune. The predicted
residual relative to Monday is 0, -2r/7, -4r/7, -6r/7, -8r/7 across the week,
where r is the crawl per trading day.

Tested on 1 September 2005 to 29 June 2006, where the crawl of -5.23 percent per
year makes the effect largest:

| | Mon | Tue | Wed | Thu | Fri |
|---|---|---|---|---|---|
| Predicted | +1.19 | +0.59 | 0.00 | -0.59 | -1.19 |
| Observed | +1.07 | +0.61 | +0.01 | -0.45 | -1.16 |

Chi-squared is 0.41 on 4 degrees of freedom. Three of the five weekday means
differ from zero by more than two standard errors.

Two further signatures agree. The residual autocorrelation has its fundamental at
lag 5, which is one week in trading days, with harmonics at 9 and 15. And the
observed excess of 1.87 basis points is closest to the 2.07 predicted for a
five-observation stepping period, against 12.42 for monthly and 37.67 for
quarterly, both of which are ruled out by size alone.

The result replicates on the two halves of that window separately, on 92 and 115
observations, with the fundamental at lag 5 in each.

**Why it does not change the crawl estimates.** With roughly 250 observations per
year, a calendar-day fit averages the trading-day steps correctly and recovers the
annual rate without bias. This is why the modern estimates land at -2.79 against
an announced -2.76. The mechanism shows up only in the within-week residual, which
in the modern era is about 0.24 basis points against a 2.83 basis point
measurement floor.

**Strength of the claim.** The evidence is strong in the 2005 and 2006 window
where the crawl is large. In the modern era the effect is too small to confirm
from any single interval. Testing 2 January 2023 to 30 June 2025 gives a weekday
slope of -0.172 against a predicted -0.171, which looks like an exact match, but
the shape is flat from Monday to Thursday with a Friday drop rather than a ramp,
and chi-squared against the ramp is 3.58 on 4 degrees of freedom. That window
cannot distinguish the two shapes. A single matching slope is not evidence, since
a slope is one number and can agree while the shape does not.

The test that would settle it is to fit every detected segment separately, fold
the residuals by weekday within each, and compare the observed ramp against the
prediction across all seventeen segments at a range of crawl rates. That has not
been done.

## 4. The framework's first year ran faster than announced

`STATUS.md` section 6 reported that the 2005 to 2007 crawl did not match, with
recovered values of about -5.1, -5.2 and -3.8 against announced values of -4.80,
-3.90 and -2.30, and called it the best open research question in the project.

Using interval boundaries taken from the announcement dates rather than from
rolling windows, the anomaly is much narrower than that. It is confined to the
first thirteen months.

| Interval | Recovered | Announced | Excess |
|---|---|---|---|
| 2005-09 to 2006-02 | -5.30 | -4.76 | -0.54 |
| 2006-02 to 2006-06 | -5.19 | -4.76 | -0.43 |
| 2006-07 to 2007-06 | -3.92 | -3.87 | -0.05 |
| 2007-07 to 2009-02 | -2.26 | -2.28 | +0.02 |
| 2009-03 to 2010-03 | -2.92 | -2.88 | -0.04 |

The announced column is adjusted for the 249.8 observations per year that a
calendar-day fit actually sees, which is a 0.9 percent correction.

From July 2006 onward the recovered rate matches the announced rate to within a
few hundredths of a percentage point and stays that way for twenty years. Before
that it runs about half a percentage point fast. Cumulatively this is roughly half
a percent of extra depreciation across the first year.

**The confounds that were ruled out.** Splitting the window in half gives -5.30
and -5.19, a spread of 0.11 percentage points against a gap to announced of 0.47.
The recovered weight is stable across the two halves at 0.654 and 0.649, so
weight drift aliasing into the crawl coefficient is not the explanation. Dropping
the first three months, which contained residual outliers as large as 21 basis
points, cuts the unexplained residual from 3.51 to 1.87 basis points but leaves
the crawl at -5.23. The settling period and the rate discrepancy are separate
things.

**The confounds that have not been ruled out.** The ground truth table gives month
precision only, so the entry for June 2005 may record the rate as announced rather
than the rate as operated, or may be dated approximately. A one-off realignment of
the parity, amortised across the first year, would look identical in the price
data. And the Bank may simply have taken time to calibrate a new mechanism. None
of these can be distinguished from published rates alone.

**Next step.** This one is documentary rather than computational. The Bank of
Botswana Annual Reports for 2005 and 2006, and the Monetary Policy Statements from
that period, describe the operation of the exchange rate framework and may state
the operational rate of crawl. A document stating -4.80 for a year when the data
says -5.2 would make this a finding. A document describing a realignment would
make it an explanation, which is also worth reporting.

**A related artefact.** The detector places a break at 2005-11-16 where the weight
goes from 0.662 to 0.651 and the crawl from -5.23 to -5.20. Neither is a
meaningful change. The segment lengths on either side are 116 and 154, close to
the 90 observation minimum, which is what a detector does when it is trying to
patch a misfit it cannot fix. With iterative reweighting the early period gets
three more such breaks at lengths 122, 90, 91 and 215. These are not policy
changes. They should be excluded from any precision score, and they are the reason
precision cannot be quoted cleanly for the early period.

## 5. The measurement floor is not constant

`PROGRESS.md` section 5 derives the rounding noise floor from first principles.
The Bank prints every column to four decimals, a uniform error of width h has
standard deviation h over the square root of twelve, and in logs the relative size
is what matters. Evaluated at an SDR of 0.0557 this gives 2.59 basis points.

That derivation is correct. Treating the result as a constant across the history
is not. The rounding grid is fixed at 0.0001 in absolute terms while the printed
SDR level moves, so the noise in log terms scales as one over SDR. The Pula has
depreciated for twenty years, so the floor has risen.

Measured across the file: **1.02 basis points in January 2001, 2.61 basis points
today.** A factor of 2.6.

This matters for the detector. A cost function assuming a constant floor judges
recent data against too tight a standard and early data against too loose a one,
which produces breaks in the recent era that are not there and hides breaks in the
early era that are. `segment.py` therefore weights each observation by its own
floor.

The floor also depends on the weight, since the noise is
hypot((1-w) s / SDR, w s / ZAR) and the SDR term dominates because SDR per Pula is
about 22 times smaller than ZAR per Pula. Assuming w = 0.5 where the truth is 0.66
overstates the floor by about 40 percent and understates chi2/dof by about a
factor of two, in exactly the early segments most worth diagnosing. The reported
chi2/dof is now computed from each segment's own fitted weight. The first segment
reads 26.91 rather than 12.69 as a result.

**Consequence for an earlier claim.** `STATUS.md` section 3 puts an upper bound of
about 1.1 basis points on how far the Bank lets the published rate drift from its
own formula. That figure was computed against a constant 2.59 basis point floor
over recent windows, so it holds for the recent era. It should not be extended
backwards. The per-segment chi2/dof column in `segments.csv` is the right basis
for restating it by era.

## 6. Robustness

Running the detector with and without iterative reweighting changes nothing from
3 July 2007 onward. Every segment boundary, weight, crawl rate and residual is
identical. Only the 2005 to 2007 block differs, and it differs by adding breaks
rather than removing them.

The results in sections 1 and 2 therefore do not depend on the choice of noise
weighting, which was the modelling decision most open to question.

## 7. Corrections to earlier documents

- `STATUS.md` section 4: the claim that all five weight changes coincide with an
  announced crawl change is wrong. Four do. The January 2025 change stands alone
  and was not announced. The bracket given for it, 8 July to 6 October 2025, is
  six months late.
- `STATUS.md` section 6: the 2005 to 2007 crawl mismatch is confined to the first
  thirteen months. July 2006 onward matches to within 0.05 percentage points.
- `STATUS.md` section 3: the 1.1 basis point bound on policy tolerance applies to
  the recent era only.
- `PROGRESS.md` section 4.4: the reasoning that weekends count because inflation
  accrues in real time is contradicted by the data. The Bank steps on trading
  days. The practical advice to use calendar days for `t` is unaffected.
- `PROGRESS.md` section 5: the 2.59 basis point floor is the value today, not a
  constant. It was 1.02 basis points in 2001.

## 8. Open questions, in order of value

1. Check the Bank's press releases for late 2024 and early 2025 for any
   announcement of a basket reweighting (section 2).
2. Check the 2005 and 2006 Annual Reports and Monetary Policy Statements for the
   operational rate of crawl (section 4).
3. Run the weekday test across all seventeen detected segments to establish
   whether trading-day stepping holds throughout the history or only in the early
   period (section 3).
4. Settle what the remaining 1.67 basis points of unexplained residual in 2005 and
   2006 is, after the weekday ramp is removed. The autocorrelation decays smoothly
   rather than periodically, which suggests persistence rather than a second
   stepping mechanism.
5. Review the 40 quarantined rows, all of which fall between 2001 and 2004. Some
   are clear errors, such as a yen value of 0.2053 against a local median of
   20.8250. Others sit only 3 to 5 percent from their local median, which is
   inside a normal week of currency movement. If real observations are being
   discarded, that affects the pre-2005 period.
