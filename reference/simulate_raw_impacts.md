# Simulate the raw impacts that studies would report

For validation: given a population whose disease combinations follow the
maximum-entropy distribution of `model` (see
[`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md)),
and known true additive impacts (and optional interactions), computes
what a study would estimate for each impact row: the crude difference
between the two groups the row compares, or the regression coefficient
of an adjusted row (see
[`cm_impacts()`](https://rasmussenphilip.github.io/deconflate/reference/cm_impacts.md)
for what a row means). With `n = NULL` the exact population values are
returned; with `n` set, `n` animals are sampled, which adds sampling
error.

## Usage

``` r
simulate_raw_impacts(
  model,
  true_impacts,
  interactions = NULL,
  rows = NULL,
  adjusted_for = NA_character_,
  n = NULL,
  joint = NULL,
  units = NULL,
  ...
)
```

## Arguments

- model:

  A
  [`cm_population()`](https://rasmussenphilip.github.io/deconflate/reference/cm_population.md)
  or
  [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md)
  (its impacts are ignored).

- true_impacts:

  Named numeric vector of true impacts, one per disease (any units).

- interactions:

  Optional named numeric vector of true interactions, named by
  combination (e.g. `c("d1;d2" = 0.5)`), in the same units.

- rows:

  Optional data frame of the rows to simulate, with column `disease` and
  optional `among` and `adjusted_for`. Default: one row per disease
  (with `adjusted_for`) and one crude row per interaction.

- adjusted_for:

  For the default rows: adjustment sets (ids separated by `";"`, or
  `"all"`), one per disease or recycled.

- n:

  Optional number of animals to sample.

- joint:

  Optional pre-computed
  [`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md)
  result. It must match the model and have converged; a joint fitted
  here that does not converge stops with class
  `deconflate_nonconvergence`.

- units:

  Optional units label for the returned impacts.

- ...:

  Not used (arguments of earlier versions give a message).

## Value

A
[`cm_impacts()`](https://rasmussenphilip.github.io/deconflate/reference/cm_impacts.md)
object with the simulated raw impacts.

## Details

Adjusting the returned raw impacts with
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
recovers `true_impacts` and `interactions` exactly (for `n = NULL`),
whenever the rows determine them.

## Examples

``` r
raw <- simulate_raw_impacts(example_supplement(), c(d1 = 2, d2 = 4, d3 = 6),
                            interactions = c("d1;d2" = 1))
raw
#>   row disease among adjusted_for measure    value study               source
#> 1   1      d1  <NA>         <NA>    <NA> 2.666100  <NA> simulated (expected)
#> 2   2      d2  <NA>         <NA>    <NA> 5.569634  <NA> simulated (expected)
#> 3   3      d3  <NA>         <NA>    <NA> 6.691225  <NA> simulated (expected)
#> 4   4   d1;d2  <NA>         <NA>    <NA> 8.082420  <NA> simulated (expected)
#>   note
#> 1 <NA>
#> 2 <NA>
#> 3 <NA>
#> 4 <NA>
deconflate(cm_model(example_supplement(), raw))$adjusted
#>   disease        kind      raw adjusted     change rows
#> 1      d1     disease 2.666100        2 -0.2498407    1
#> 2      d2     disease 5.569634        4 -0.2818199    1
#> 3      d3     disease 6.691225        6 -0.1033032    1
#> 4   d1;d2 interaction 8.082420        1         NA    1
```
