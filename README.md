# Tekanyetso

Recovering the Bank of Botswana's crawling-band fixing rule from published data,
detecting when its parameters change, and forecasting the next published rate.

## Week 1 status

Data source settled, formula recovered, `§2.2` answered. See "Findings" below.

```
pip install numpy
python3 fetch.py            # pulls BoB, stores raw + clean, validates
python3 estimate.py         # recovers w and the crawl, writes rolling.csv
python3 estimate.py --log   # also appends tomorrow's prediction (run daily)
```

## The data

Everything comes from **one** Bank of Botswana table:
`https://www.bankofbotswana.bw/content/exchange-rates`, CSV export at
`/export/exchange-rates.csv`. It carries CHN, EUR, GBP, USD, ZAR, SDR and YEN
per Pula, on ~254 pages of 25 rows, back to roughly 2001.

The IMF and ECB are **not needed**. BoB publishes its own SDR-per-Pula column
alongside its rand column, on the same dates in the same convention, so the
calendar-alignment and unit-conversion work budgeted in build-notes §1 does not
arise. Both sides of the basket come from the same file.

Convention: **foreign currency per Pula.** ZAR ≈ 1.23 (a Pula buys 1.23 rand),
SDR ≈ 0.056.

**Column order differs between the CSV export and the HTML table.** The export
is `Date,CHN,EUR,GBP,USD,SDR,YEN,ZAR`; the web page is
`Date,CHN,EUR,GBP,USD,ZAR,SDR,YEN`. Both label their headers correctly, so
reading by name is safe and reading by position silently swaps ZAR and SDR —
the two columns the model depends on. `fetch.py` reads by name.

**Inversion is caught on the column median, not per row.** A per-row range check
cannot detect an inverted ZAR: 1/1.23 = 0.81 sits inside any plausible band,
because a Pula and a rand are worth about the same. The median of the whole
column cannot hide like that.

## Data quality

The table contains real faults, so `fetch.py` screens rather than trusts. A bad
row is **quarantined to `data/quarantine.csv`, never fatal** — one mistyped row
in 2019 must not stop today's prediction from being logged. Only structural
failures (too few rows, a whole column inverted) abort the run.

Three fault modes, all confirmed present:

| Fault | Example |
|---|---|
| Column contamination — a cell holds another column's value that day | 02 Nov 2023 `SDR = 0.0605`, which is that day's GBP; 16 Dec 2019 `SDR = 0.0929`, that day's USD |
| Isolated typo | 05 Apr 2019 `CHN = 0.0939` against a 0.63 neighbourhood |
| Duplicate date with different values | 13 Aug 2019 appears twice; the second copy's values match late Feb 2019 |

**All three survive a per-row range check** — a contaminated SDR of 0.0605 is a
perfectly plausible exchange rate. None survives comparison against the row's
own neighbours, so the screen is a relative-deviation test against a local
median, with per-column tolerances (SDR tightest, since it is a basket and
barely moves; ZAR loosest). The neighbourhood is bounded to 21 calendar days so
it never compares across the gaps in the series (BoB skipped 26 Sep – 02 Oct
2025, and 17 – 22 Jul 2026).

Contaminated SDR matters more than it looks: SDR ≈ 0.0557 enters the log basket
with a sensitivity of about 5.18 bp per 1e−4, so an 0.0605 reading is roughly
800 bp of error on the index — enough to wreck every window containing it,
against a residual the model expects to be ~2.5 bp.

Storage: `raw/` holds every download verbatim, named by date and content hash,
never overwritten. `data/rates.csv` is append-only with a `first_seen` column;
if BoB ever revises a published value, `fetch.py` reports it rather than
silently backfilling. That is most of the benefit of bitemporal storage for
almost none of the cost, because each fetch returns the full history.

`data/crawl_changes.csv` is the ground-truth scoring key: 15 announced rate-of-
crawl changes since 2005, from the Bank's Current Exchange Rate Framework page.
Dates there are month-precision only, so changepoint scoring needs a tolerance
band of roughly ±1 month.

## The model

The peg holds a weighted log basket on a linear path:

    w·log(ZAR_per_Pula) + (1−w)·log(SDR_per_Pula) = a + b·t

With `u = log(ZAR_pP) − log(SDR_pP)` this is one ordinary least squares fit:

    log(SDR_pP) = a + b·t − w·u + e

so `w = −coef(u)` and the annual crawl is `365·b`. Sum-to-one is imposed by
construction rather than estimated, which is what removes the degree of freedom
that build-notes §5 worried about.

`t` is in **calendar** days. The crawl is an annualised continuous drift applied
per calendar day; using trading days puts the estimate out by about 45 percent
(365/252 ≈ 1.45).

## Findings

**1. It is the central fixing, not a market rate.** Build-notes §2.2 flagged this
as a fork between two different projects. On the most recent 25 trading days the
recovered parameters are `w = 0.496` and crawl `−2.76 %/yr`, against announced
values of 50/50 and −2.76. Residual RMSE is **2.47 bp**.

**2. The residual is at the measurement floor.** BoB rounds every column to four
decimals. On SDR ≈ 0.0557 that is a coarse grid: uniform rounding to 1e−4
implies about 2.6 bp of noise on the basket index. The observed 2.47 bp residual
is *at* that floor, so the peg is holding as tightly as this data can resolve.
The published number is mechanical. This is the good fork — the original plan
works as written, and band width becomes a bounded-above question rather than
the centre of the project.

**3. Rounding, not collinearity, is the binding constraint.** The condition
number of the demeaned regressors runs about 2 over recent windows — very well
identified, because the rand/SDR spread moves a lot. But the crawl is only
−2.76/365 ≈ **0.76 bp per day**, well under the 2.6 bp quantisation noise. The
crawl is invisible day-to-day and only emerges once cumulative drift beats the
rounding, which takes tens of days. Expect the weight to be pinned quickly and
the crawl to need long windows. Retune §5's remedies accordingly: this argues
for longer windows and for treating short-window crawl estimates as unreliable,
not for regularising the weight.

## Open question for week 2

Given rand and SDR, the fixing is deterministic to ~2 bp — so "predict tomorrow's
Pula rate" is only a real forecast if it is made without tomorrow's rand and SDR.
Two targets, and they are different claims:

- **Rule prediction** — tomorrow's log basket index equals today's plus one day
  of crawl. Needs no FX forecast, is a claim about the *Bank*, and is testable
  from this repo alone. This is what `estimate.py --log` currently writes, and
  it can start today.
- **Level prediction** — the headline "tomorrow's USD/BWP". This needs an
  external market snapshot plus a timing story: *what* market data does BoB use,
  and *when* is the rate published relative to it? If BoB fixes off the previous
  close, this is nearly deterministic and the log becomes a precision
  demonstration. If it fixes off same-morning rates, there is a genuine lead-time
  problem worth solving.

Settling the publication timing is the first thing to do in week 2. Until then
the rule prediction is the honest one to log.

## Known limitations

- ~~The single-shot CSV export may only return the first page.~~ **Resolved
  10 Sep 2026:** the export returns the entire history in one response. The
  pagination fallback was removed — `?page=N&_format=csv` on `/content/...`
  serves HTML, not CSV, so the old fallback would have parsed zero rows and
  looked like an empty history. `fetch.py` now asserts a minimum row count
  instead.
- Basket **weight** change dates are only shown as a chart on the framework
  page, not a table. They need reading off press releases, unlike the crawl
  changes.
- The 2.47 bp result is one recent 25-day window. It needs re-running across the
  full history before it is a claim.
