# Tekanyetso — Week 1 Addendum

**10 September 2026.** Supersedes `PROGRESS.md` where the two conflict, on the
same rule `PROGRESS.md` applies to the earlier plan documents: this was written
after the live endpoint was actually hit, and `PROGRESS.md` was not.

---

## 1. The endpoint works, and returns everything

`PROGRESS.md §8.1–8.2` treats the export as untested and flags the risk that it
serves only page one. Both settled:

- `https://www.bankofbotswana.bw/export/exchange-rates.csv?page&_format=csv`
  returns the **entire history in a single response**, most recent first. The
  site paginates at 25 rows across 255 pages (last page index 254), so expect
  roughly **6,400 rows back to about 2001**.
- **The pagination fallback was broken and has been removed.** The URL it used,
  `/content/exchange-rates?page=N&_format=csv`, ignores `_format` and returns
  the HTML page. `parse()` would have found zero rows there — so had the export
  genuinely been truncated, the fallback would have quietly produced an *empty*
  history rather than a paginated one. Replaced with an explicit minimum-row
  assertion that aborts loudly.

## 2. Column order differs between the CSV and the web page

| Source | Order |
|---|---|
| CSV export | `Date, CHN, EUR, GBP, USD, SDR, YEN, ZAR` |
| HTML table | `Date, CHN, EUR, GBP, USD, ZAR, SDR, YEN` |

Both label their headers correctly. Reading by name is safe; reading by position
silently swaps **ZAR and SDR** — precisely the two columns the model consumes.
`fetch.py` reads by name and now says so in a comment, because the failure would
be silent and catastrophic.

## 3. The inversion test was aimed at the wrong statistic

`PROGRESS.md §2.1` credits the per-row range assertion as the inversion test. It
is not sufficient. An inverted ZAR gives 1/1.23 ≈ **0.81**, which sits inside the
`[0.5, 3.0]` band — because a Pula and a rand happen to be worth about the same,
the one currency in the table where inversion is invisible to a range check.

Fixed by testing the **median of each column** against a tight band. A whole
inverted column cannot hide from its own median. The per-row ranges are kept,
but their job is now catching single wild cells, not inversion.

## 4. The table has real faults, and none of them are range violations

Three fault modes, all confirmed in the live data:

| Fault | Confirmed instances |
|---|---|
| Column contamination | 02 Nov 2023 `SDR = 0.0605` (that day's GBP); 16 Dec 2019 `SDR = 0.0929` (that day's USD) |
| Isolated typo | 05 Apr 2019 `CHN = 0.0939`; 19 Mar 2019 `CHN = 0.0933`; 06 Aug 2019 `GBP = 0.0818` (that day's EUR) |
| Duplicate date, different values | 13 Aug 2019 appears twice; the second copy's values match late Feb 2019 |
| Unsorted rows | 09 Jan 2024 sits inside the Dec 2023 block; 29 May 2023 inside June |

Two consequences.

**The old `validate()` would have crashed on the first run.** It raised
`SystemExit` on any duplicate date, and duplicates exist. `PROGRESS.md §8.1`
says "expect to debug it" — this is the specific bug.

**A contaminated SDR is worse than it looks.** SDR ≈ 0.0557 enters the log basket
with sensitivity ~5.18 bp per 1e−4 (`PROGRESS.md §5`), so an 0.0605 reading is
roughly **800 bp** of error on the index, against a residual the model expects to
be ~2.5 bp. One such row poisons every window containing it.

**New policy: quarantine, never fatal.** Bad rows go to `data/quarantine.csv`
with a reason string; the run continues. One mistyped row in 2019 must not stop
today's prediction from being logged (`PROGRESS.md §9`, non-negotiable 4). Only
structural failure — too few rows, a column inverted — aborts.

The screen is a relative-deviation test against a **local median**, per-column
tolerances (SDR tightest at 3%, since it is a basket and barely moves; ZAR
loosest at 6%), neighbourhood bounded to **21 calendar days** so it never
compares across the real gaps in the series (BoB skipped 26 Sep – 02 Oct 2025,
and 17 – 22 Jul 2026). Verified on a fixture of the known-bad rows: 6/6 caught,
0 false positives.

## 5. The week-1 numbers need error bars, and one of them changes meaning

`PROGRESS.md §3` reports `w = 0.496` against 0.50 and crawl `−2.76` against
−2.76, on 25 days. `estimate.py` now returns standard errors. From 200 synthetic
draws at the true parameters, rounded on BoB's 4dp grid (95% interval coverage
measured at 90–94%, so mildly optimistic but roughly honest):

| Window | se(`w`) | se(crawl) |
|---|---|---|
| 25 obs | 0.007 | **0.25 %/yr** |
| 60 obs | 0.003 | 0.07 %/yr |
| 90 obs | 0.002 | 0.04 %/yr |
| 250 obs | 0.0008 | 0.008 %/yr |

Reading the week-1 result against these:

- **`w = 0.496 ± 0.007`.** The announced 0.500 is well inside. The 0.004 gap is
  not evidence of anything — do not report it as a deviation.
- **crawl `= −2.76 ± 0.25`.** Still consistent with the announced −2.76, but the
  95% interval spans roughly −2.3 to −3.3. **The exact digit-for-digit match was
  luck.** It is corroboration, not the precision demonstration it reads as.
  `PROGRESS.md §3` should be read with that caveat.

This is `§5`'s "use long windows for the crawl" made quantitative. Slope
precision improves as roughly n^(−3/2), so the crawl tightens fast: 25 → 90
observations buys about a 7× tighter estimate.

**It also settles the week-2 gate in advance.** 45/55 and 50/50 differ by 0.05 in
`w`; at a 90-day window se(`w`) ≈ 0.002, so the two regimes sit about **25
standard errors apart**. They will separate unmistakably. If they do not, the
problem is a data fault, not collinearity — check `quarantine.csv` first.

## 6. Independent confirmation of the rounding floor

The synthetic series is built with a known `w` and crawl and then rounded to 4dp
exactly as BoB rounds. Recovered residual RMSE: **2.46–2.69 bp** across window
lengths, against the 2.59 bp predicted analytically in `PROGRESS.md §5`. The
derivation is right, and the observed 2.47 bp on real data is the rounding floor
and nothing else.

## 7. Zeros are missing data, and only two columns are load-bearing

Found on the first live run of `fetch.py` (10 Sep 2026), which aborted with
`column CHN has median 0.0000`.

BoB did not quote the yuan for most of the history and wrote **`0.0000`** into
those cells rather than leaving them blank. More than half the CHN column is
zeros, so its median is zero. A zero exchange rate is not a price; it is a
missing value wearing a number.

Two fixes, and the second matters more than the first:

- **Zeros parse to missing**, and every median, range and spike test now
  operates on present values only. Judging a column on a median that includes
  its missing-data zeros says nothing about whether the column is inverted — it
  only says the column is sparse.
- **`MODEL_COLS = ["ZAR", "SDR"]`.** Those two are the only columns anything
  reads. A row missing CHN is perfectly usable and is now kept; a row missing
  ZAR or SDR is quarantined. The other five columns are carried for completeness
  and can no longer block anything. `estimate.py` skips unusable rows and says
  how many.

`fetch.py` now prints per-column coverage with a first-observation date, so the
real start of each series is visible rather than assumed. **Check the ZAR and SDR
lines on the first successful run.** The `README` claim of history "back to
roughly 2001" is about the *table*; if either model column starts later than
2005, the crawling-band period is shorter than the ground-truth key assumes and
some of the 15 announced changepoints are unscoreable. That would not sink the
project, but it changes what the detector can be graded on and should be
recorded before any scoring is done.

## 8. Still open

- **Publication timing** (`PROGRESS.md §7`) — untouched. Still the first real
  task of week 2, and still what decides whether "level prediction" is a
  precision demo or a genuine forecasting problem.
- **Basket weight change dates** — still only on a chart, still need reading off
  press releases.
- **Nothing has been logged.** `fetch.py` and `estimate.py --log` now need to run
  daily on a machine with network access. This remains the only thing that
  cannot be caught up later.
