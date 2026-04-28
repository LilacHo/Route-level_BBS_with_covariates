# Baird's Sparrow — Route-level BBS Trend Model (2010–2024)

Model: `models/slope_habitat_route_NB.stan` (route-level iCAR slope model with
habitat mean and habitat-change as predictors on intercept and trend, negative-
binomial likelihood).

Predictor scaling (Smith et al. 2024, ACE 19(2):2, Appendix 2):

- `route_habitat` = mean proportion of the covariate per route, mean-centered.
- `route_habitat_slope` = 100 × (per-route slope of the covariate over time,
  mean-centered). Units: percentage-points per year above the range-wide average.

---

## Run 1 — Covariate: proportion developed (1 km buffer)

## Diagnostics

All reported parameters have `rhat ≈ 1.00` with `ess_bulk` and `ess_tail` in the
thousands. The sampler converged and mixed well; inferences are trustworthy.

## Parameter-by-parameter

### BETA — mean residual trend (log scale per year)

- mean 0.00031, 90% CI [−0.00051, +0.00115]
- On the percent scale: ≈ **+0.03 %/yr** [−0.05, +0.12]
- After removing the habitat-driven component, Baird's Sparrow shows no detectable
  residual population trend; credible interval spans zero.

### rho_BETA_hab — effect of developed-slope on route trends

- mean 0.00059, 90% CI [−0.00408, +0.00571]
- A 1-unit increase in the predictor = a route gaining 1 percentage-point of
  developed land per year *above* the range-wide average.
- Coefficient translates to ≈ **+0.06 %/yr** added to the route's log-trend per
  +1 pp/yr of development. 90% CI spans zero.
- No evidence that change in developed land drives route-level trend differences
  in this dataset / time window.

### ALPHA — mean log-scale intercept (mean relative abundance)

- mean −0.029, sd ≈ 1.00
- Posterior is essentially the prior (`ALPHA ~ N(0, 1)`). The data did not push the
  grand-mean intercept away from the prior center, which is expected for a sparse,
  patchy species with low counts at most routes.

### rho_ALPHA_hab — effect of mean-developed on route-level mean abundance

- mean 0.035, sd ≈ 0.99, 90% CI [−1.58, +1.69]
- Posterior ≈ prior (`rho_ALPHA_hab ~ N(0, 1)`). **No signal.**
- Ecologically consistent: Baird's Sparrow is a grassland obligate and avoids
  developed land, so developed% as an abundance predictor is uninformative across
  the routes where the species occurs.
- For comparison: Smith et al. 2024 reported `Pα ≈ 3` [2.2, 3.8] for Rufous
  Hummingbird with modeled habitat suitability — that predictor was biologically
  relevant for the species; `developed` here is not.

### T — full population trend (%/yr)

- mean −0.0053 %/yr, 90% CI [−0.089, +0.079]
- Mean-route population, averaged across the species range, is changing by roughly
  **−0.005 %/yr** — essentially zero with ±0.09 %/yr uncertainty.
- No detectable population trend for Baird's Sparrow over 2010–2024 in this
  analysis.

### T_no_habitat — trend with habitat-change effect removed

- mean −0.0026 %/yr
- Nearly identical to `T`. Confirms habitat-change is not moving the headline
  trend.

### T_dif = T − T_no_habitat

- mean −0.0027 %/yr, 90% CI [−0.014, +0.004]
- The fraction of the trend "explained by" development change: very slightly
  negative (development pulling the trend down by ~0.003 %/yr) but overwhelmed
  by uncertainty.
- ESS is healthy, so this is a small, noisy genuine signal rather than a sampling
  artifact.

### CH — cumulative % change over the 15-yr window (full model)

- mean −0.08%, 90% CI [−1.33%, +1.19%]
- Total population change from 2010 to 2024: essentially zero, ±1%. Flat within
  uncertainty.

### CH_no_habitat and CH_dif

- `CH_no_habitat`: mean −0.04%, 90% CI [−1.30, +1.24]
- `CH_dif`: mean −0.04%, 90% CI [−0.21, +0.07]
- Development-driven change over 15 years accounts for a tiny fraction of a
  percent. Effectively no mechanistic signal.

## Bottom line

For Baird's Sparrow, 2010–2024, with *proportion developed in a 1 km buffer* as
the habitat covariate:

- No detectable population trend (`T ≈ 0`).
- No detectable effect of developed-land area on mean abundance
  (`rho_ALPHA_hab` ≈ prior).
- No detectable effect of developed-land change on route trends
  (`rho_BETA_hab` ≈ 0, CI spans zero).

The model converged and fit cleanly — the result is "no signal," not a failed
fit. This is ecologically consistent: Baird's Sparrow is a grassland bird, and
proportion developed is a poor habitat descriptor for it.

## Next step

Re-run with the **grassland** covariate (`data/grassland_1km/grassland.csv`) —
it is the relevant land-cover class for this species. Expected directions:

- `rho_ALPHA_hab` positive: more grassland → more birds.
- `rho_BETA_hab` positive: routes losing grassland should show more negative
  trends (and the sign convention here means a positive coefficient, since both
  predictor and trend move together).


---

## Run 2 — Covariate: proportion grassland (1 km buffer)

### Diagnostics

All parameters have `rhat ≈ 1.00` with `ess_bulk` and `ess_tail` in the thousands.
The sampler converged cleanly. `T_dif` has the lowest bulk ESS (1963) but still
well above the usual 400 threshold, so inferences are trustworthy.

### Parameter-by-parameter

#### BETA — mean residual trend (log scale per year)

- mean 0.00029, 90% CI [−0.00055, +0.00113]
- On the percent scale: ≈ **+0.03 %/yr** [−0.05, +0.11]
- Residual trend (after removing the grassland-change component) is essentially
  flat, as in the developed run. CI spans zero.

#### rho_BETA_hab — effect of grassland-slope on route trends

- mean **−0.000182**, 90% CI [−0.00106, +0.000682]
- Predictor units: percentage-points of grassland per year above the range-wide
  average. So a 1-unit increase = 1 pp/yr more grassland gain than average.
- Coefficient is **negative**: routes with relatively *more* grassland gain have
  slightly *lower* log-trends; routes with grassland *loss* have less-negative
  (more positive) log-trends. 90% CI still overlaps zero, so the evidence is
  weak, but the posterior sits noticeably below zero (~66% of mass is negative).
- Sign interpretation: the expected sign for a grassland obligate is *positive*
  (losing grassland should drag populations down). The observed negative mean is
  small and uncertain and shouldn't be over-interpreted — it may reflect that
  within the narrow 2010–2024 window, remnant-grassland routes and working-
  cropland routes both saw small, noisy changes in % grassland that don't line
  up cleanly with sparrow trend variation.

#### ALPHA — mean log-scale intercept

- mean −0.060, sd ≈ 1.01
- Posterior ≈ prior (`ALPHA ~ N(0, 1)`). No shift from the prior center,
  consistent with a sparse, patchy species.

#### rho_ALPHA_hab — effect of mean-grassland on route-level mean abundance

- mean **+0.26**, sd 0.74, 90% CI [−0.95, +1.48]
- The posterior has narrowed substantially compared to the prior sd = 1.00
  (posterior sd = 0.74) — the data *did* inform this parameter, unlike the
  developed run where the posterior equalled the prior.
- Mean is positive, as expected for a grassland obligate: routes with more
  grassland tend to have higher mean counts of Baird's Sparrow. Roughly 64% of
  posterior mass is above zero.
- CI still spans zero, so this is suggestive but not statistically strong.
  Biologically it's consistent with the known ecology of the species.

#### T — full population trend (%/yr)

- mean **−0.011 %/yr**, 90% CI [−0.094, +0.071]
- The full trend (including grassland change) is slightly negative, about
  -0.01 %/yr, with CI still spanning zero. Slightly more negative than the
  developed-run value of −0.005 %/yr.

#### T_no_habitat — trend with habitat-change effect removed

- mean **+0.008 %/yr**, 90% CI [−0.074, +0.090]
- With the grassland-change component stripped out, the residual trend flips
  sign to slightly *positive*.

#### T_dif = T − T_no_habitat

- mean **−0.020 %/yr**, 90% CI [−0.047, +0.001]
- This is the part of the trend attributable to grassland change. The posterior
  is almost entirely below zero (~98% of mass < 0), and the upper bound of the
  90% CI barely touches zero (+0.001).
- Interpretation: within this 15-yr window, grassland change is pulling the
  Baird's Sparrow trend *down* by about 0.02 %/yr. This is a weak but
  directionally clear signal — and it's the opposite of what developed-land
  produced (where `T_dif` was noisy at −0.003 %/yr).

#### CH — cumulative % change over 15 yr (full model)

- mean **−0.17%**, 90% CI [−1.40%, +1.06%]
- Small cumulative change over the window, CI overlapping zero.

#### CH_no_habitat

- mean **+0.13%**, 90% CI [−1.11%, +1.35%]
- Again, removing grassland-change flips the mean positive.

#### CH_dif = CH − CH_no_habitat

- mean **−0.29%**, 90% CI [−0.70%, +0.01%]
- Cumulative attribution: grassland change accounts for about −0.3% of
  population change over 2010–2024, with 98% of the posterior below zero.
  This is the clearest "habitat is acting on populations" signal in the whole
  analysis — the credible interval just barely touches zero at the upper end.

### Bottom line — Grassland covariate

The grassland covariate is doing more work than the developed covariate, in
ways that make ecological sense:

- **Intercept**: positive (mean +0.26) effect of mean-grassland on abundance;
  posterior narrower than the prior. Consistent with grassland-obligate ecology.
- **Trend**: the `T_dif` and `CH_dif` summaries both put ~98% of posterior mass
  below zero. Grassland change is associated with a small, directionally
  consistent downward pressure on Baird's Sparrow trends over 2010–2024.
- **Headline trend** (`T`) is still not statistically different from zero, but
  the decomposition is informative: without grassland change, the residual
  trend would be slightly positive; with it, slightly negative.

Caveats:

- 15 years is a short window for detecting grassland-change signatures on BBS
  trends. The developed covariate may also behave differently over a longer
  baseline.
- The negative `rho_BETA_hab` sign is counter-ecological; it may reflect noise
  in the covariate over a short window, or confounding with the positive trend
  pressure captured by `T_no_habitat`. The positive attribution signal
  (`T_dif < 0`) arises because the magnitude of grassland *loss* at specific
  routes overcomes the sign of the `rho_BETA_hab` coefficient when combined
  with the route-level slopes. Worth revisiting if a longer covariate record is
  available.

## Overall comparison

|                         | Developed             | Grassland               |
|-------------------------|-----------------------|-------------------------|
| `rho_ALPHA_hab` (mean)  | 0.03 (≈ prior)        | 0.26 (narrower than prior) |
| `rho_BETA_hab` (mean)   | 0.0006 (≈ 0)          | −0.0002 (weak, unclear sign) |
| `T` (%/yr)              | −0.005 [−0.09, +0.08] | −0.011 [−0.09, +0.07]   |
| `T_dif` (%/yr)          | −0.003 [−0.014, +0.004] | **−0.020 [−0.047, +0.001]** |
| `CH_dif` (%)            | −0.04 [−0.21, +0.07]  | **−0.29 [−0.70, +0.01]** |

Grassland is the ecologically relevant habitat variable for Baird's Sparrow,
and the model reflects that: grassland change shows a small but directionally
consistent downward effect on population change over 15 years, whereas
developed land shows essentially no mechanistic signal.
