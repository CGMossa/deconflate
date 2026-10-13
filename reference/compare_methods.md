# Compare the model with the published approximation

Runs the package's model and the earlier approximations on the same
inputs and tabulates the results side by side. For additive impacts:

- `"deconflate"`: the model of
  [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md);

- `"published"`: the proportional approximation used in Rasmussen et al.
  (2022, 2024) (eq. 16 of the 2022 paper). It is not part of the model;
  it is kept here (and in the `reproduce_*()` functions) to compare
  earlier results with the exact solution. It needs one crude row per
  disease and a crude pair row for every pair, can mask incompatible
  inputs, and is undefined when its denominator is zero.

## Usage

``` r
compare_methods(
  model,
  methods = NULL,
  event_model = FALSE,
  overall_risk = NULL,
  n_draws = 0,
  seed = NULL,
  ...
)
```

## Arguments

- model:

  A
  [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md).

- methods:

  Methods to compare (default: all for the kind of impacts).

- event_model, overall_risk:

  As in
  [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md).

- n_draws:

  Number of draws (default 0: point estimates only).

- seed:

  Optional random seed for the draws.

- ...:

  Passed to
  [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)'s
  model (e.g. `joint`, `pool`, or arguments of
  [`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md)).

## Value

A `cm_comparison` object with `impacts` (raw and adjusted values, one
column per method), `change` (relative change from raw), `long` (long
format with sign-change flags), `totals` (naive and adjusted aggregate
per method; for event impacts, the attributable risk), `diagnostics`,
`failed` (methods that could not be run, or gave an undefined result,
with reasons), `undefined`, `methods`, `results` and, with draws,
`draws` (summary of each quantity by method, rejections).

## Details

For event impacts (`event_model = TRUE`): `"deconflate"` (the snapshot
model of
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)),
`"first_order"` (its log-linear approximation, from the pairwise tables)
and `"published"` (the approach of Rasmussen et al. 2024: `HR - 1`
adjusted with eq. 16), each with the attributable risk at
`overall_risk`. The first-order and published methods need one hazard
(or rate) ratio per disease and every pair.

Methods that cannot be run on the inputs are reported with the reason,
and results with non-finite values are kept apart in `undefined`.

With `n_draws`, the methods are also applied to the same draws of the
uncertain inputs (a draw is rejected if any method fails on it), and the
means and 95% intervals are reported in `draws`.

## Examples

``` r
compare_methods(example_supplement())
#> <cm_comparison> methods: published, deconflate
#> Impacts: yield
#> Units: %
#> 
#> Adjusted values:
#>  disease raw published deconflate
#>       d1 2.5      2.07       2.14
#>       d2 5.0      3.70       3.39
#>       d3 7.5      6.75       6.93
#> 
#> Totals:
#>      method raw_sum adjusted_total
#>   published     2.5          2.111
#>  deconflate     2.5          2.109
```
