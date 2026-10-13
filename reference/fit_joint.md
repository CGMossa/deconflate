# Fit the maximum-entropy distribution of disease combinations

The model needs the probability of every combination of diseases. The
association rows do not determine it (for three or more diseases there
are 2^n - 1 free probabilities), so the package uses the maximum-entropy
distribution that matches the disease probabilities and every
association row: the distribution that assumes as little as possible
beyond the inputs. Pairs that no row mentions (unknown pairs) are not
constrained: their association is whatever the fit implies given the
others. A pair given an odds ratio of 1 is constrained to be
independent.

## Usage

``` r
fit_joint(
  model,
  tol = 1e-10,
  max_iter = 10000L,
  max_diseases = 20L,
  backend = c("exact", "sampled"),
  n_samples = 50000L,
  n_chains = 1000L,
  burn_in = 50L,
  fit_iter = 300L,
  calibrate = TRUE,
  seed = NULL,
  start = NULL,
  pool = c("random", "fixed")
)
```

## Arguments

- model:

  A
  [`cm_population()`](https://rasmussenphilip.github.io/deconflate/reference/cm_population.md)
  or
  [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md).

- tol:

  Convergence tolerance: maximum absolute deviation of any event or
  disease probability from its target.

- max_iter:

  Maximum number of full sweeps.

- max_diseases:

  Safety limit on the number of diseases for the exact backend (the
  table has 2^n cells).

- backend:

  `"exact"` (default) or `"sampled"`.

- n_samples:

  Sampled backend: number of combinations drawn after calibration.

- n_chains:

  Sampled backend: number of parallel Gibbs chains.

- burn_in:

  Sampled backend: Gibbs sweeps discarded before calibration and before
  the final draws.

- fit_iter:

  Sampled backend: maximum number of moment-matching iterations.

- calibrate:

  Sampled backend: adjust the sample weights to match the targets
  exactly?

- seed:

  Optional random seed (sampled backend).

- start:

  Optional earlier exact fit to start from (e.g. the fit at the central
  inputs, when refitting for a draw of the inputs). It is used only when
  it has the same diseases and association rows (only their values may
  differ); the result is then the same maximum-entropy distribution,
  found in fewer sweeps.

- pool:

  `"random"` (default) or `"fixed"`: how several association rows about
  the same quantity are pooled (see
  [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)).

## Value

A `cm_joint` object: `cells` (0/1 matrix of combinations), `prob`
(probability of each combination), `converged`, `iterations`,
`max_residual`, `backend`, the `targets` it was fitted to (used by
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
to check that a supplied joint matches the model), the pooling of
association rows (`pooling`) and, for the sampled backend,
`diagnostics`.

## Details

Every association row fixes the probability of one event: for a pair
row, the probability of both diseases; for a row within a subgroup, of
both diseases within it; for a combination, of all its diseases together
(see
[`cm_associations()`](https://rasmussenphilip.github.io/deconflate/reference/cm_associations.md)).
Several rows about the same event are pooled (see
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)).
For rows whose event probability depends on other probabilities (e.g. an
odds ratio within a subgroup, whose margins are themselves fitted), the
targets are recomputed at the current fit until they no longer move.

## Exact backend

`backend = "exact"` enumerates all 2^n combinations. Starting from
independence (the product of the disease probabilities; slide 12 of the
2023 seminar), iterative proportional fitting on events rescales the
combination probabilities until every event probability and disease
probability is matched. The result is the maximum-entropy distribution
with those event probabilities. Enumeration limits this backend to about
20 diseases.

## Sampled backend

`backend = "sampled"` fits the same model without enumerating
combinations, for larger numbers of diseases, when every association row
is a crude pair row:

1.  Calibration: the model's parameters (a main effect per disease and
    an interaction per constrained pair) are fitted by Monte Carlo
    moment matching. Gibbs sampling from many parallel chains estimates
    the probabilities and pairwise tables implied by the current
    parameters, and the parameters are updated (by the differences in
    logits and log odds ratios between the targets and the estimates)
    until the estimates match the targets within Monte Carlo error.

2.  Sampling: `n_samples` combinations are drawn from the calibrated
    model.

3.  Optionally (`calibrate = TRUE`, the default), the weights of the
    sampled combinations are adjusted by IPF so that the disease
    probabilities and pairwise tables match the targets exactly on the
    sample (raking). The remaining Monte Carlo error then affects only
    higher-order structure.

The result has the same form as the exact one (`cells` are the distinct
sampled combinations, `prob` their weights). Its `diagnostics` give the
constraint residuals of the raw sample and after raking, R-hat and
effective sample sizes across chains, and the Monte Carlo errors.
Results have Monte Carlo error: compare runs with different seeds, and
with the exact backend where n allows.

## Convergence

If the association rows cannot all hold together, the fit does not
converge. Non-convergence after `max_iter` sweeps suggests, but does not
prove, infeasibility (see
[`check_feasibility()`](https://rasmussenphilip.github.io/deconflate/reference/check_feasibility.md)).
For the sampled backend, a calibration that does not converge is
reported as unresolved, not as evidence that the inputs are infeasible.

## Examples

``` r
j <- fit_joint(example_supplement())
j
#> <cm_joint> 3 diseases, 8 combinations (exact backend)
#>   Converged: TRUE after 28 sweeps (max residual 8.28e-12)
#>   Association rows: 3
# \donttest{
js <- fit_joint(example_supplement(), backend = "sampled", n_samples = 20000, seed = 1)
js$diagnostics$summary
#>   n_samples n_chains n_unique fit_iterations fit_converged residual_sample
#> 1     20000     1000        8             17          TRUE     0.001220607
#>   residual_calibrated max_rhat  min_ess    max_mcse
#> 1        5.075046e-11 1.004167 17279.93 0.003048884
# }
```
