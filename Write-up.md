# 1. Introduction {#introduction}

Currency runs the world, as it facilitates trade and commerce; when a currency collapses, the world feels the effects. Hence, countries have various ways of controlling their currency to protect themselves from ruin.

## 1.1 How Countries Control Their Currency {#how-countries-control-their-currency}

Control rests with the central bank: the financial institution that runs a country\'s money. This includes its currency, which means issuing notes, setting interest rates and, in this case, setting the exchange rate. There are three main systems.

### 1.1.1 Floating Rate {#floating-rate}

This is a system in which the government does not regulate the currency. Instead, it lets buyers and sellers settle on a value minute by minute, like a share price. The central bank sets no target value at all. If it buys or sells currency, it does so to calm a disorderly market, not to hold the rate at a particular level. A bank that does commit to keeping the rate inside a range or band is running a **managed float**, which sits between this system and a peg.

### 1.1.2 Fixed Peg {#fixed-peg}

A monetary policy in which a central bank sets a strict, unchanging value for its national currency.

### 1.1.3 Crawling Peg {#crawling-peg}

A managed exchange-rate system in which a central bank links the value of its domestic currency to a major foreign currency or a basket of currencies, making small, frequent and predictable adjustments over time.

*Table 1: Comparison of exchange-rate systems*

<table>
<thead>
<tr>
<th><strong>System</strong></th>
<th><strong>Pros</strong></th>
<th><strong>Cons</strong></th>
<th><strong>Countries using it</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td>Floating rate</td>
<td><p>Monetary independence: interest rates can target inflation at home</p>
<p>Shock absorber</p>
<p>No reserves spent holding a level</p></td>
<td><p>Reserve strain</p>
<p>Management cost</p>
<p>Speculative targets</p></td>
<td>USA, UK</td>
</tr>
<tr>
<td>Fixed peg</td>
<td><p>Imported price stability: inflation tracks the anchor country</p>
<p>Certainty for traders and lenders</p></td>
<td><p>No autonomy</p>
<p>Reserve drain</p></td>
<td>Saudi Arabia, Hong Kong, Eritrea</td>
</tr>
<tr>
<td>Crawling peg</td>
<td><p>An inflation anchor that still allows the rate to move</p>
<p>Trade competitiveness as costs rise</p>
<p>Market predictability: the drift is announced</p></td>
<td><p>Policy constraints</p>
<p>Speculative attacks</p></td>
<td>Botswana, China, Bolivia</td>
</tr>
</tbody>
</table>

## 1.2 Botswana\'s Crawling Peg {#botswanas-crawling-peg}

The subject of this write-up, Botswana, uses a crawling peg to help control inflation and keep its exports competitive. The crawl is a slow, deliberate weakening of the Pula, announced in advance; it is not a market movement.

To achieve this, the Bank of Botswana (the country\'s central bank) ties the Pula not to a single currency but to a basket: a mix of two things.

- **The rand (ZAR)**, because South Africa is Botswana\'s biggest trading partner.
- **The Special Drawing Right (SDR)**, which is not a country\'s currency at all. It is a unit of account created by the International Monetary Fund, and it is itself a fixed mix of five major currencies: the US dollar, the euro, the Chinese renminbi, the Japanese yen and the British pound. Using it is a compact way of saying "the rest of the world" in one number.

So "the Pula is worth a stable amount of the basket" means that the Pula\'s value is measured half against South Africa and half against the wider world.

At the time of writing (23 September 2026), one Pula buys 1.2266 rand and 0.0556 SDR.

# 2. The Problem and the Project {#the-problem-and-the-project}

## 2.1 What the Bank Does Not Publish {#what-the-bank-does-not-publish}

The Bank of Botswana controls the peg, and the crawl is currently set at −2.76 percent a year. However, the Bank does not publish everything. In particular, it does not publish:

- the exact arithmetic that turns the rule into each day\'s numbers;
- the exact dates on which it changes the weights or the crawl (it gives only the month for crawl changes, and shows weight changes only as a chart);
- how closely it actually sticks to its own rule.

## 2.2 What This Project Does {#what-this-project-does}

The aim of the project is to recover all three from the published numbers, then to predict tomorrow\'s number before it appears, and to publish every prediction so that anyone can check the claim. *(Not yet complete.)*

## 2.3 The Rounding Error and Its Size {#the-rounding-error-and-its-size}

Even when the Bank follows its rule exactly, the published points do not sit exactly on the line, because the numbers are rounded to four decimal places. Rounding introduces an error, but it is nevertheless unavoidable when the figures are published. The size of this error can be calculated using the standard deviation of a uniform distribution of width *h*:

> σ = h / √12 (1)

With four decimal places the rounding step is h = 0.0001, so the rounding error has a standard deviation of about 0.0000289.

Because the project works in logarithms (see Section 3.1), the error must also be expressed in logarithmic terms. For a published value *P* with a small error *e*:

> ln(P + e) − ln(P) ≈ e / P (2)

So the same absolute rounding error matters more for a small number:

*Table 2: Effect of the rounding error on each column*

| **Column** | **Value** | **0.0000289 ÷ value** | **In basis points (bp)** |
|------------|-----------|-----------------------|--------------------------|
| SDR        | 0.0555    | 0.00052               | 5.20                     |
| ZAR        | 1.23      | 0.0000235             | 0.23                     |

The SDR column is 22 times more sensitive, simply because its value is 22 times smaller. This is why the SDR column dominates everything about the error in this project.

## 2.4 The Noise Floor {#the-noise-floor}

The **noise floor** is the baseline level of random, meaningless variation below which economic insights or trading signals cannot be detected.

The basket uses half of each column, and the two rounding errors are independent, so they combine the way independent errors always do: in quadrature (the square root of the sum of squares):

> floor = √\[(0.5 × 5.20)² + (0.5 × 0.23)²\] = √(2.60² + 0.12²) ≈ 2.60 bp (3)

This is the noise floor: about **2.6 bp**. No fit can ever do better, because that much scatter is built into the published numbers. If the project\'s residuals come out at about 2.6 bp, the conclusion is not "the model is mediocre"; it is "the Bank is following its rule as exactly as this data can reveal".

**The floor is not a constant.** The rounding grid is fixed at 0.0001, but the SDR value has fallen over twenty years as the Pula crawled down, so the same grid now represents a larger relative error. In January 2001 the SDR value was 0.1427 and the floor was **1.02 bp**; today it is **2.61 bp**. Treating the floor as fixed would make the change detector too suspicious about recent years and too forgiving about early ones.

## 2.5 Is the Fit Good? {#is-the-fit-good}

The goodness of fit can be checked using a chi-squared test:

- If the chi-squared value per degree of freedom is about 1, the residuals are exactly the size that the rounding predicts, and the model describes everything else.
- If it is well above 1, there is something in the data that the model does not describe.

# 3. Recovering the Rule {#recovering-the-rule}

## 3.1 Why the Project Works in Logarithms {#why-the-project-works-in-logarithms}

The crawl is a **percentage** change per year. Percentage changes are multiplied rather than added, and multiplication quickly becomes awkward. To solve this, logarithms are used to convert multiplication into addition. This is essentially converting a graph into linear form: if a quantity *P* falls by a fixed percentage each day, then ln(*P*) falls by a **fixed amount** each day, so a graph of ln(*P*) against time is a straight line.

*Note:* a change of 0.01% corresponds to ln(1.0001) ≈ 0.0001.

## 3.2 Basis Points {#basis-points}

A **basis point** (bp) is a unit of measure equal to **one-hundredth of one percent (0.01%, or 0.0001)**. It is used mainly to express changes or differences in [*interest rates*](https://www.investopedia.com/terms/b/basispoint.asp) and bond yields. Basis points are used throughout this write-up.

## 3.3 The Basket: Which Average to Use? {#the-basket-which-average-to-use}

The basket is half rand and half SDR, but there are two ways to average the two: the arithmetic average and the geometric average. A good test is to see whether the answer survives when the rates are read the other way round (Pula per unit instead of units per Pula).

*Table 3: Basket rates per Pula, and the same rates flipped*

| **Currency** | **Per Pula** | **Flipped: Pula per unit** |
|--------------|--------------|----------------------------|
| Rand         | 1.2300       | 0.8130                     |
| SDR          | 0.0555       | 18.0180                    |

### 3.3.1 Arithmetic Average {#arithmetic-average}

The arithmetic average is the sum of a collection of numbers divided by how many numbers there are; here it adds the two rates and halves the result. For example:

> (1.2300 + 0.0555) × 0.5 = 0.6428 (4)

The problem with this method is that if the rates are read the other way round, the two routes give different answers, so the result depends on an arbitrary choice of direction:

*Table 4: Arithmetic average, calculated both ways*

| **Route** | **Working** | **Answer** |
|----|----|----|
| Average, then flip | (1.2300 + 0.0555) ÷ 2 = 0.6428, then 1 ÷ 0.6428 | **1.5557** |
| Flip, then average | (0.8130 + 18.0180) ÷ 2 | **9.4155** |

### 3.3.2 Geometric Average {#geometric-average}

The geometric average (or geometric mean) multiplies the numbers together and takes the *n*-th root, which makes it ideal for tracking compound growth rates and investment returns. For two numbers it takes the square root of their product:

> √(1.2300 × 0.0555) = 0.2613 (5)

This method has the added advantage that it still holds true when the rates are flipped, so the direction in which the currencies are read has no effect on the result:

*Table 5: Geometric average, calculated both ways*

| **Route**          | **Working**                                  | **Answer** |
|--------------------|----------------------------------------------|------------|
| Average, then flip | √(1.2300 × 0.0555) = 0.2613, then 1 ÷ 0.2613 | **3.8274** |
| Flip, then average | √(0.8130 × 18.0180)                          | **3.8274** |

For this reason, the project uses the geometric average.

## 3.4 The Rule {#the-rule}

Let:

- *Z* = rand per Pula on a given day
- *S* = SDR per Pula on a given day
- *lz* = ln(*Z*)
- *ls* = ln(*S*)
- *w* = the weight on the rand
- *t* = time in days

The rule says that the log of the basket falls in a straight line:

> w·lz + (1 − w)·ls = a + b·t (6)
>
> w·lz + ls − w·ls = a + b·t
>
> ls + w(lz − ls) = a + b·t

Now let *u* = *lz* − *ls*:

> ls + w·u = a + b·t
>
> ls = a + b·t − w·u (7)

The best line is then found by least squares, that is, by choosing the values that minimise the sum of squared residuals.

*Note on u:* *Z* is rand per Pula and *S* is SDR per Pula. Dividing them cancels the Pula completely, leaving rand per SDR: a pure world-market quantity with no Botswana in it. Whatever the Bank of Botswana does, it cannot move *u*. The input is therefore genuinely independent of the thing being explained. In statistics this property is called **exogeneity**, and it is why the estimates here can be trusted more than those in a typical economics regression.

The output is three numbers, and they are read as follows:

- the constant is *a*, the starting level, which is of no interest;
- the coefficient on *t* is *b*, the crawl per day; multiply it by 365 and by 100 to get percent per year;
- the coefficient on *u* is −*w*, so the **weight is minus that coefficient**.

## 3.5 Results {#results}

On recent data this gives *w* = 0.496 and a crawl of −2.76 percent a year, against the Bank\'s announced 0.50 and −2.76. That agreement, obtained from the published numbers alone, is the first real result of the project.

# 4. Detecting Changes in the Rule {#detecting-changes-in-the-rule}

## 4.1 The Problem {#the-problem}

The economy and markets are always in motion, so it would be unreasonable to keep the crawl the same throughout; inflation can grow too quickly. The crawl therefore has to be changed from time to time, which means the history is made up of **segments**: stretches where the rule was constant, separated by **changepoints**. The Bank does not say when the segments begin and end, but the changepoints can be found, and the segments determined, by fitting the model to the dated data.

The difficulty is that a model with more segments always achieves a lower score, but the lowest score is not always the best answer. To fix this, a penalty is introduced:

> total cost = cost of all segments + penalty × (number of changepoints) (8)

Now a changepoint has to earn its place by improving the fit by more than the penalty costs. The project uses a standard choice called BIC, 3 × ln(*n*), where 3 is the number of parameters that a new segment adds, and then multiplies it by 10. The multiplier was settled by checking the detected dates against the Bank\'s announcements.

## 4.2 The Exact Method {#the-exact-method}

Using dynamic programming, the problem is solved for the shortest stretch of data first, and each longer solution is then built from the solutions already found:

> F\[τ\] = min over every possible s of { F\[s\] + cost(segment from s to τ) + penalty } (9)

This method is slow, however, because it tries every earlier start. That is roughly *n*²/2 segment costs; for the 5,254 days of Botswana data, that is about 12.9 million.

Two things make it manageable. Each segment cost is computed in a small, fixed number of steps using **running totals** (cumulative sums are kept, so any segment\'s totals are a simple subtraction). There is also a **minimum segment length** of 90 days, which also stops the detector from inventing tiny segments that are too short to pin down a crawl anyway.

Still, *n*²/2 grows badly: ten times more data means a hundred times more work.

## 4.3 PELT: Discarding Starts That Cannot Win {#pelt-discarding-starts-that-cannot-win}

**PELT** (Pruned Exact Linear Time) is the fast version, and the important word is **exact**: it gives the identical answer, not an approximation.

The key observation is that the slow method keeps reconsidering starts that have already lost badly. PELT keeps a list of the starts still in contention and permanently drops the rest.

The rule for dropping a start *s*, checked after computing F\[τ\], is:

> if F\[s\] + cost(s to τ) \> F\[τ\], then s can never win again (10)

This relies on one property of the cost, which is worth stating because it is the hinge of the whole argument:

> cost(s to t) + cost(t to T) ≤ cost(s to T) (11)

In words: splitting a segment in two never fits worse than leaving it whole. That is clearly true for least squares, because the two-piece fit can always copy the one-piece fit by using the same parameters in both halves, so it can only do better or the same.

**The subtlety.** The argument above says that *s* is beaten by a rival, namely τ. But with a minimum segment length, τ is not always available as a start. So the project does not drop a dead start immediately; it keeps it for another 90 steps first.

Is that caution necessary? The project tested it on 16,987 randomly generated hard cases. The delayed version was exact every time. The textbook version, which drops starts immediately, **got the wrong answer 546 times** (3.2 percent). It never failed on the actual Botswana data, which is precisely why a proof matters more than testing.

# 5. Findings {#findings}

## 5.1 What the Detector Found {#what-the-detector-found}

Four results came out of this, and they are the substance of this write-up.

1.  **It finds every announced change.** The Bank published a list of fifteen crawl rates since 2005; fourteen changes fall inside the window studied. The detector, which never sees that list, finds all fourteen. Eight land within a week of the start of the announced month, and thirteen within a month.
2.  **A basket reweighting nobody announced.** The weight changed from 0.45 to 0.50 on 6 January 2025 while the crawl stayed put, six months before the next announced change. (This is only "unannounced" once the Bank\'s press releases have been checked, which is still on the to-do list.)
3.  **The Bank steps the peg on trading days, not calendar days.** This was found by looking at what was left over after fitting: the leftovers formed a sawtooth that reset every week, with exactly the size and shape predicted if the Bank moves the rate once per published day while the model assumes every calendar day.
4.  **The first year ran fast.** From July 2006 onwards, the recovered crawl matches the announced one to within a few hundredths of a percentage point, for twenty years. In the framework\'s first thirteen months, it runs about half a percentage point faster than announced.

# 6. Tracking and Predicting {#tracking-and-predicting}

## 6.1 The Kalman Filter, Without the Matrices {#the-kalman-filter-without-the-matrices}

The detector assumes that the weight and crawl are constant and then jump. A **Kalman filter** assumes instead that they drift slightly every day, and it updates its estimate each time a new number arrives.

The idea in one paragraph: you have a guess about something you cannot see directly, and a measurement of it that is noisy. Each day you do two things. **Predict:** move your guess forward using what you know about how it changes, and become slightly less certain, because it may have drifted. **Update:** compare your prediction with the new measurement, and move your guess part of the way towards the measurement. How far you move depends on which you trust more. If your prediction is solid and the measurement noisy, you barely move; if the reverse is true, you move a lot. That trust ratio is the **Kalman gain**.

A satnav does exactly this: it combines where it expects you to be (from your last position and speed) with a noisy GPS fix, and trusts each according to its uncertainty.

Two details are worth noting. First, the project\'s filter tracks the rule\'s implied rate rather than the basket index, because when the weight changes, the index jumps by about 15 percent for purely definitional reasons while the Pula itself does not move. Second, the filter\'s own settings were fitted on 2005 to 2010 only and then frozen, so its performance after 2010 is an honest test rather than a fit to the same data.

Its main result: the filter says the same thing as the detector about the weight, including the January 2025 change, and it also settled the trading-day question over the whole history, favouring trading-day stepping by a decisive margin.

## 6.2 What "Predicting" Means Here {#what-predicting-means-here}

Once the rule is known, tomorrow\'s Pula rate depends on tomorrow\'s rand and SDR, which are world-market numbers that nobody can forecast. So predicting the actual rate would just be predicting the currency markets, which is not the point.

Instead, the project predicts the **rule**: tomorrow\'s basket index will equal today\'s plus one step of crawl. No market forecast is needed, and the claim can be checked against the Bank\'s own published numbers the next day. It is a claim about the Bank\'s behaviour, not about the market.

Every prediction is written to a file at 09:00, before the Bank publishes, and pushed to a public repository in the same minute. The file is never edited, only added to.

## 6.3 The Backtest and the Trap It Avoids {#the-backtest-and-the-trap-it-avoids}

A **backtest** replays history: pretend you were making a prediction every day since 2005, using only what was known at the time, and see how you would have done. It gives thousands of test predictions instead of the handful in the live log.

The trap is **lookahead**: accidentally using information from the future. It is the equivalent of marking your own exam after seeing the answers, and it is easy to do without noticing. If you estimate the weight using all the data up to 2026 and then test a prediction made in 2009, you have used the future.

The project\'s backtest proves that it avoids this, rather than merely claiming it. At several dates it reruns everything twice: once with all later data **deleted**, and once with all later data **replaced by random numbers**. Every prediction made before that date must come out bit-for-bit identical. If anything had peeked, the numbers would change. All five models pass.

Results over 3,725 predictions from 2011 onwards: the model currently writing the live log is out by 3.56 bp on average, and the Kalman filter by 2.68 bp. Both are close to the rounding floor, which means the rule is being recovered about as well as the published data allows.

*Table 6: Next-day prediction error with no lookahead (out-of-sample RMS error, 2011--2026)*

| **Model**                       | **RMS error (bp)** |
|---------------------------------|--------------------|
| Kalman filter                   | 2.68               |
| Trading-day crawl               | 3.54               |
| Live model (calendar-day crawl) | 3.56               |
| No crawl (baseline)             | 3.61               |
| Trading-day, fitted anchor      | 3.75               |

# 7. Mistakes and the Limits of What Has Been Shown {#mistakes-and-the-limits-of-what-has-been-shown}

## 7.1 Six Things That Did Not Work {#six-things-that-did-not-work}

A paper that reports only successes is less useful, and less believable, than one that reports how the successes were arrived at. Here are six from this project.

1.  **A constant noise floor.** The first version of the change detector treated the 2.6 bp floor as a fixed number. It is not: it was 1.02 bp in 2001, because the SDR value was larger then. A detector using a constant floor judges recent data against too tight a standard, so it invents changes in recent years and misses them in early ones. This was found by checking a derived quantity against the data rather than trusting the derivation.
2.  **A changepoint that was not one.** The detector places a break on 16 November 2005, where the weight moves from 0.662 to 0.651 and the crawl from −5.23 to −5.20: no real change. The surrounding fit is poor (chi-squared per degree of freedom of 26.9), and the detector is spending a changepoint patching a misfit. A detector will always do this when the model is wrong in a region, and reporting it is part of reporting the method honestly.
3.  **The textbook pruning rule.** Dropping a beaten start immediately, as the standard description of PELT does, is wrong when there is a minimum segment length. It lost the true answer in 3.2 percent of hard cases. It never failed on the real data, so only the proof or a deliberate search for hard cases could have caught it.
4.  **Anchoring predictions on the fitted value.** Starting each prediction from the regression\'s own fitted value, rather than from the published number, sounds better because it avoids one day\'s rounding. It is worse: the 90-day regression lags every policy change, so near a break its error is 5.20 bp, against 3.64 bp elsewhere. The Kalman filter gets the benefit without the lag, because it updates daily.
5.  **Reading the penalty off a curve.** The textbook way to choose the penalty is to plot the number of changes found against the penalty and look for a plateau. There was no clean plateau, because real changes vary enormously in size, so each increase in the penalty removes the next-weakest real change. The penalty was settled against the Bank\'s announcements instead, and the paper says so.
6.  **An unchecked count.** An earlier draft said that eleven of the fourteen detected dates landed within about a week. Nobody had counted. It was eight. Every number now has a command that reproduces it.

## 7.2 What the Project Has Not Shown {#what-the-project-has-not-shown}

The live prediction has not been running for long; it becomes evidence only over months. Until then, the backtest carries the accuracy claims.

- The backtest tests a claim about the Bank, not a forecast of the currency market. It would be worthless as a trading result and is not offered as one.
- The Kalman filter\'s structure was chosen after looking at the data. Only its fitted values are out of sample, not the decision to use that model.
- Two findings rest on documents nobody has checked yet: the January 2025 reweighting and the first year\'s crawl rate.
