# Changepoint detection — method

*Week 3. Draft section for `docs/method.md`. Companion to `segment.py`.*

---

## 1. The problem

The Bank changes two parameters — the basket weight `w` and the rate of crawl `b`
— at reviews it does not fully announce. `STATUS.md §4` established that it
changes both at the same review and announces only the crawl.

Formally: partition the history into segments so that total within-segment misfit
is minimised, with a penalty `β` charged per extra segment:

    minimise   Σ_segments [ C(segment) + β ]   −  β

over all partitions. Subtracting one `β` means a single segment is charged
nothing, so `β` is purely the price of a *change*.

Within a segment the model is the one `estimate.py` already fits:

    ls = a + b·t − w·u + e            u = lz − ls,   t in years

## 2. The segment cost

`C` is the weighted residual sum of squares of that regression, whitened by a
known per-observation noise level, so it is a chi-square. Three properties matter:

**It is additive.** Total cost is the sum of segment costs with nothing linking
them. This is what makes the dynamic program valid at all, and it is the
precondition for PELT's pruning argument in week 4.

**It is computable in O(1).** A naive implementation refits an OLS per candidate
segment, making the DP O(n³) — unrunnable on 6,324 points. Instead precompute
cumulative sums of the 3×3 matrix `X'X`, the 3-vector `X'y` and the scalar `y'y`.
Any segment's normal equations are a difference of two cumulative sums, solved in
closed form via the adjugate of a symmetric 3×3. The DP touches the raw data
exactly once, at construction.

**Its scale is known a priori.** See §3. This is unusual and valuable: most
changepoint problems must estimate the noise level from the same residuals they
are testing, which couples the penalty to the answer.

### 2.1 Numerical note — cancellation, and an exact fix

A segment's RSS is recovered as `y'y − β'X'y`, and on this data `y'y` runs about
6×10⁴ times larger than the RSS being extracted from it. Five significant digits
are lost to cancellation before the DP ever sees the number.

The fix is exact rather than approximate. Subtract a single global fit and
segment the *residual* instead. Because every segment fit already spans the
columns of `X`, subtracting `X·β_global` shifts that segment's coefficients by
exactly `−β_global` and leaves its residual — hence its RSS — bit-for-bit
identical. Same problem, two orders of magnitude less cancellation; `β_global` is
added back when coefficients are reported.

Measured effect: worst-case cost error over 300 random segments fell from 1.0×10⁻⁴
to 6.6×10⁻⁷, against a penalty term of 21. The right tolerance here is in *cost
units*, not relative units — the DP only ever compares costs against `β`, so an
error is harmless until it is a material fraction of one penalty.

### 2.2 Every segment gets a free intercept

Structural, not lazy. At a weight change the modelled quantity changes
*definition*: `w·lz + (1−w)·ls` at `w = 0.45` is a different index from the same
expression at `w = 0.50`, so its level jumps mechanically even though the Pula
itself does not. Forcing continuity across segments would model an artefact of
the parameterisation. It also keeps segments statistically independent, which §2
requires.

---

## 3. Correction to `PROGRESS.md §5` — σ is not constant

`PROGRESS.md §5` derives the rounding noise floor from first principles: BoB
prints to four decimals, a uniform error of width `h` has sd `h/√12`, and in logs
what matters is the *relative* size `h/√12/level`. Evaluated at today's SDR of
0.0557 this gives **2.59 bp**, confirmed empirically in `STATUS.md §3`.

That derivation is correct. Treating the answer as a constant across the history
is not. The rounding grid is fixed at 1e−4 in **absolute** terms while the
printed SDR level moves, so the noise in log terms goes as `1/SDR`:

| Year | SDR per Pula (approx) | Floor |
|---|---|---|
| 2001 | 0.092 | 1.57 bp |
| 2008 | 0.080 | 1.81 bp |
| 2015 | 0.069 | 2.08 bp |
| 2026 | 0.0557 | 2.59 bp |

The floor has risen roughly 65% across the sample, purely because the Pula has
depreciated and BoB kept printing four decimals. A homoscedastic cost function
would therefore judge recent data against too tight a floor and early data
against too loose a one — over-detecting breaks in the recent era and missing
them in the early one. This is an artefact generator, not a rounding detail.

Each observation is therefore weighted by its own `1/σᵢ`:

    σᵢ = hypot( (1−w)·(h/√12)/SDRᵢ ,  w·(h/√12)/ZARᵢ )

with `w = 0.5`. The cost becomes a true chi-square, which has two payoffs: the
BIC penalty is a clean `3·log(n)`, and **chi2/dof reads directly as "how tightly
is the peg held relative to measurement precision"** — 1.0 means exactly at the
floor, and the excess above 1.0 is the quantity `STATUS.md §3` estimates as a
≈1.1 bp upper bound on policy tolerance.

**Consequence for `STATUS.md §3`:** that 1.1 bp bound was computed against a
constant 2.59 bp floor over recent windows, so it stands for the recent era. It
should not be extended backwards. Re-derive per segment from `chi2_dof` in
`segments.csv`.

---

## 4. The penalty

BIC: `β = 3·log(n)`, three being the parameters a changepoint adds.

This is a default, not a result. The penalty's meaning depends entirely on σ
being right — mis-specifying σ by 35% changes the effective penalty by 80%, which
was visible in development as a spurious break that vanished once the weighting
of §3 was applied. Since the model is only approximately true (rounding is
assumed uniform and unbiased; the residual is assumed pure noise), report a
**penalty sweep** — breaks found against `β` — and look for a plateau, rather
than quoting one BIC number. `--sweep` does this.

---

## 5. Validation

`segment.py --self-test`, four tests:

1. **O(1) cost == direct weighted least squares.** Worst cost error 6.6×10⁻⁷
   over 300 random segments, 3×10⁻⁸ of one penalty term. Coefficients agree to
   1.3×10⁻¹².
2. **DP == exhaustive search** over all partitions on a small series. Identical
   cuts, identical cost to 10⁻⁶. This is the exactness guarantee PELT must
   reproduce in week 4 — the comparison is already wired.
3. **No break in, no break out.** Five seeds of single-regime data, zero
   spurious breaks at plain BIC.
4. **chi2/dof == 1 when σ is right.** Measured 0.991. This validates §3's noise
   model end to end.

Synthetic recovery, five regimes of 400 observations, generated by building data
that obeys the peg exactly and then **rounding it the way BoB rounds it** (not
adding Gaussian noise — the real mechanism, so the test is faithful):

| True break | Detected | Error |
|---|---|---|
| 400 (crawl only, −0.16 → 0.00) | 400 | 0 |
| 800 | 800 | 0 |
| 1200 | 1200 | 0 |
| 1600 | 1601 | +1 |

Zero false positives at plain BIC, untuned. Recovered parameters match the
generating values to within 0.001 in `w` and 0.01 %/yr in the crawl.

The first break is the hard case and the one worth quoting: a crawl change of
0.16 percentage points is 0.04 bp per day against a 3.6 bp noise floor —
invisible in any rolling window — and it is found to the exact observation,
because the DP accumulates the drift over the whole segment rather than trying
to see it locally.

## 6. Scaling

| n | cost evaluations | wall | evals / n² |
|---|---|---|---|
| 500 | 52,092 | 0.03 s | 0.208 |
| 1,000 | 338,342 | 0.09 s | 0.338 |
| 2,000 | 1,660,842 | 0.28 s | 0.415 |
| 4,000 | 7,305,842 | 0.98 s | 0.457 |

Converging on n²/2 as expected. Extrapolated to the full 6,324-row history:
~20M evaluations, ~2.5 s. This is the baseline for week 4's runtime plot.

**The comparison must be fair.** A vectorised DP against a scalar-loop PELT would
produce a meaningless speedup. Both reduce to "for each τ, evaluate the cost over
an array of candidate starts" — only the array differs, all of them for DP and a
pruned set for PELT — so the same vectorised code path serves both. The primary
metric is **cost evaluations**; wall clock is reported alongside.

---

## 7. Open for week 4

- **PELT**, sharing the cost object and the candidate-array interface above.
  The pruning condition needs care in the presence of `min_seg`: candidates must
  only be pruned once they are older than the minimum segment length.
- **Scoring** against `data/crawl_changes.csv`. Announced dates are
  month-precision, so precision/recall needs a ±1 month tolerance band; the five
  weight changes from `STATUS.md §4` are currently bracketed to 90 days and
  should be re-dated by this detector first, then used as ground truth.
- **Errors in variables.** `u` appears as a regressor and contains the same
  rounding noise as the dependent variable, which biases `w` downward. The
  attenuation is of order `var(noise)/var(u)` ≈ 27 bp² / ~10⁵ bp², i.e. about
  3×10⁻⁴ — negligible, and consistent with the observed `w = 0.499 ± 0.003`
  against an announced 0.500. Recorded as considered and dismissed, not ignored.
- **The 2005–2007 discrepancy** (`STATUS.md §6`). Run the detector restricted to
  that window; if it wants breaks where none were announced, that is evidence the
  early framework differs structurally rather than the data being noisy.
