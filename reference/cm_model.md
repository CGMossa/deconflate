# Combine a population with an impact table

Combine a population with an impact table

## Usage

``` r
cm_model(
  population,
  impacts = NULL,
  associations = NULL,
  distributions = NULL,
  ...
)
```

## Arguments

- population:

  A
  [`cm_population()`](https://rasmussenphilip.github.io/deconflate/reference/cm_population.md).
  For convenience, a
  [`cm_diseases()`](https://rasmussenphilip.github.io/deconflate/reference/cm_diseases.md)
  object can be given instead, together with `associations`.

- impacts:

  A
  [`cm_impacts()`](https://rasmussenphilip.github.io/deconflate/reference/cm_impacts.md)
  object.

- associations:

  Used only when `population` is a
  [`cm_diseases()`](https://rasmussenphilip.github.io/deconflate/reference/cm_diseases.md)
  object; see
  [`cm_population()`](https://rasmussenphilip.github.io/deconflate/reference/cm_population.md).

- distributions:

  Optional named list of distributions
  ([distributions](https://rasmussenphilip.github.io/deconflate/reference/distributions.md))
  for uncertain inputs, used when
  [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
  runs draws. Names are keys: `"prob:<disease>"` (the disease's `value`,
  on the scale it was entered), `"associations:<row>"` and
  `"impacts:<row>"` (a row of those tables, by its number in `row`). For
  convenience, `"assoc:<d1>:<d2>"` names the crude row of a pair and
  `"impact:<id>"` the row of a disease or combination, when there is
  exactly one.
  [`cm_read_inputs()`](https://rasmussenphilip.github.io/deconflate/reference/cm_read_inputs.md)
  fills this in from the `dist` columns of the tables. When `population`
  is a model and `distributions` is not given, its distributions are
  kept, except those of impacts that are replaced.

- ...:

  Not used (arguments of earlier versions give a message).

## Value

A `cm_model` object (which is also a `cm_population`).

## Examples

``` r
m <- cm_model(
  cm_diseases(c("d1", "d2", "d3"), c(0.10, 0.15, 0.20)),
  cm_impacts(c("d1", "d2", "d3"), c(2.5, 5, 7.5), units = "%"),
  associations = cm_associations(c("d1", "d2"), c("d2", "d3"), c(2, 3)),
  distributions = list("impact:d1" = dist_normal(2.5, 0.5))
)
m
#> <cm_model>
#>   Diseases: 3 (d1, d2, d3)
#>   Association rows: 2 (OR: 2)
#>   Disease pairs: 3 (2 with a crude pair row, 1 unknown)
#>   Impacts: (unlabelled) [%] (additive)
#>   Impact rows: 3 (0 adjusted, 0 within subgroups)
#>   Uncertain inputs (with a distribution): 1
```
