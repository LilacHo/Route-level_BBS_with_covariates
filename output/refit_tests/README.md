# Route-Level BBS Refit Tests

This note documents the model-refit tests used to improve convergence and
posterior predictive diagnostics for the route-level BBS habitat models.

## Purpose

The original negative-binomial route-level model had two separate issues:

1. Some large species did not mix well, especially variance hyperparameters
   such as `sdalpha`, `sdbeta`, and habitat-effect scale terms.
2. Posterior predictive checks showed too few zero counts in replicated data,
   especially for abundant but patchily detected species.

The tests here evaluate whether stricter data filtering, standardized
covariates, stronger initialization, and zero-inflated likelihoods improve
diagnostics.

## Test Species

Four species were selected to span the diagnostic gradient while keeping run
time manageable:

| Code | Species | Reason included |
|---|---|---|
| BOSP | Botteri's Sparrow | Clean control species |
| CCLO | Chestnut-collared Longspur | Marginal original convergence |
| UPSA | Upland Sandpiper | Poor original convergence |
| BOBO | Bobolink | Failed original convergence and strong zero problem |

## Shared Abundance Model

For observation `i` on route `s = route[i]`, the expected count is modeled on
the log scale as

```text
E[i] = alpha[s]
       + beta[s] * (year[i] - fixedyear)
       + obs[observer[i]]
       + eta * firstyr[i]
```

The abundance likelihood in the baseline test model is

```text
count[i] ~ NegBinomial2Log(E[i], phi)
```

where `phi` is the negative-binomial dispersion parameter.

The route intercept and slope are decomposed as

```text
alpha[s] = alpha_resid[s] + alpha_hab[s]
beta[s]  = beta_resid[s]  + beta_hab[s]
```

with

```text
alpha_resid[s] = ALPHA + sdalpha * alpha_raw[s]
beta_resid[s]  = BETA  + sdbeta  * beta_raw[s]
```

and habitat components

```text
alpha_hab[s] =
  (rho_ALPHA_hab + sdrho_alpha_hab * rho_alpha_raw_hab[s])
  * route_habitat[s]

beta_hab[s] =
  (rho_BETA_hab + sdrho_beta_hab * rho_beta_raw_hab[s])
  * route_habitat_slope[s]
```

The route-level spatial effects use the intrinsic CAR prior

```text
p(b) proportional to exp(-0.5 * sum_edges (b[node1] - b[node2])^2)
```

with a soft sum-to-zero constraint:

```text
sum(b) ~ Normal(0, 0.001 * nroutes)
```

## Refit Settings

The test runner uses:

```text
min_max_route_years = 3
standardized route_habitat and route_habitat_slope
adapt_delta = 0.99
max_treedepth = 15
stable non-zero initial values for scale parameters
tighter scale priors for intercept and habitat-effect SDs
```

These changes mainly address mixing problems caused by weakly informed routes
and broad scale priors.

## Tested Likelihoods

### 1. NB_test

This is the improved negative-binomial model:

```text
count[i] ~ NegBinomial2Log(E[i], phi)
```

It tests whether convergence can be improved without changing the observation
model.

### 2. ZINB_test

This model adds one species-wide structural-zero probability:

```text
Pr(count[i] = 0) =
  zi + (1 - zi) * NegBinomial2LogPMF(0 | E[i], phi)

Pr(count[i] = y > 0) =
  (1 - zi) * NegBinomial2LogPMF(y | E[i], phi)
```

with

```text
zi = inv_logit(zi_logit)
zi_logit ~ Normal(-2, 1.5)
```

This model helped slightly, but a single zero-inflation probability was too
blunt for several species.

### 3. ZINB_route_test

This model allows structural-zero probability to vary by route:

```text
logit(zi_route[s]) =
  zi_intercept + sd_zi_route * zi_route_raw[s]
```

and observation `i` uses `zi_route[route[i]]`.

This helped CCLO, but UPSA and BOBO still underpredicted zeros.

### 4. ZINB_route_mu_test

This model allows zero inflation to vary by route and expected log abundance:

```text
logit(zi_obs[i]) =
  zi_intercept
  + sd_zi_route * zi_route_raw[route[i]]
  + zi_log_mu * E[i]
```

The prior is

```text
zi_log_mu ~ Normal(-0.5, 0.75)
```

A negative `zi_log_mu` means low expected abundance is allowed to have higher
structural-zero probability. This directly targets the remaining zero-count
misfit for UPSA and BOBO.

## Diagnostics Used

Convergence targets:

```text
R-hat < 1.01
bulk ESS > 400
tail ESS > 400
```

Sampler-health targets:

```text
0 divergent transitions
0 or very few max-treedepth hits
E-BFMI > 0.2
```

Posterior predictive checks compare observed and replicated:

```text
mean(count)
sd(count)
proportion(count == 0)
max(count)
```

Bayesian p-values near 0 or 1 indicate poor fit.

LOO is reported when available, but high Pareto-k counts mean PSIS-LOO should
be treated as supportive rather than definitive.

## Findings So Far

### Sampler Health

All recent ZINB route and route-plus-abundance fits had:

```text
0 divergences
0 max-treedepth hits
healthy E-BFMI values
```

So the main remaining limitations are mixing and model fit, not HMC geometry.

### Convergence

The route-plus-abundance model (`ZINB_route_mu_test`) produced:

| Code | max R-hat | min bulk ESS | Interpretation |
|---|---:|---:|---|
| BOSP | 1.003 | 1147 | Pass |
| CCLO | 1.011 | 439 | Very close, slight R-hat issue |
| UPSA | 1.026 | 211 | Does not fully pass |
| BOBO | 1.011 | 656 | Very close, slight R-hat issue |

The abundance-dependent zero model improves fit but adds parameters, so UPSA in
particular may need longer sampling or a simpler final zero model.

### Posterior Predictive Checks

The major improvement is in the zero-count PPC:

| Code | Route-ZINB p_zero | Route+mu ZINB p_zero | Interpretation |
|---|---:|---:|---|
| BOSP | 0.868 | 0.796 | Both acceptable; simpler model is enough |
| CCLO | 0.104 | 0.384 | Route+mu improves zeros |
| UPSA | 0.000 | 0.416 | Route+mu fixes the zero failure |
| BOBO | 0.000 | 0.436 | Route+mu fixes the zero failure |

This strongly supports the idea that the remaining zero inflation depends on
expected abundance, not only route identity.

The spread check also improved or remained acceptable:

```text
BOSP p_sd = 0.592
CCLO p_sd = 0.648
UPSA p_sd = 0.912
BOBO p_sd = 0.784
```

### LOO

Compared with route-only ZINB, route-plus-abundance ZINB improved expected
predictive accuracy for CCLO, UPSA, and BOBO:

| Code | Route-ZINB elpd | Route+mu ZINB elpd | Change |
|---|---:|---:|---:|
| BOSP | -130.1 | -130.6 | -0.5 |
| CCLO | -1715.9 | -1694.2 | +21.6 |
| UPSA | -6865.2 | -6852.6 | +12.7 |
| BOBO | -16077.2 | -16020.6 | +56.5 |

However, Pareto-k warnings remain high for CCLO, UPSA, and BOBO, so LOO should
not be the only basis for model selection.

## Current Recommendation

Use the simpler improved NB or route-ZINB model for species that already pass
PPC zero checks, such as BOSP.

For species with strong excess-zero problems, especially UPSA and BOBO, the
best candidate so far is `ZINB_route_mu_test`, because it fixes the zero PPC
and improves LOO. Before using it as final, refit difficult species with longer
sampling:

```text
iter_warmup = 3000
iter_sampling = 3000
chains = 4 or 6
adapt_delta = 0.99
```

If UPSA remains poorly mixed, consider simplifying the zero model by keeping
the abundance-dependent term but reducing route-level zero effects, or fitting
only route-level abundance effects for species where the route random zero
effect is weak.

