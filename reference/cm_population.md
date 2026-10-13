# Describe the population: diseases and their associations

A population holds the disease probabilities and the association rows.
Combine it with one impact table in
[`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md).
Every row is checked against the diseases (see
[`cm_associations()`](https://rasmussenphilip.github.io/deconflate/reference/cm_associations.md)
for what a row means).

## Usage

``` r
cm_population(diseases, associations = NULL, ...)
```

## Arguments

- diseases:

  A
  [`cm_diseases()`](https://rasmussenphilip.github.io/deconflate/reference/cm_diseases.md)
  object.

- associations:

  Optional
  [`cm_associations()`](https://rasmussenphilip.github.io/deconflate/reference/cm_associations.md)
  object.
  [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
  needs at least one association; the sensitivity tools (e.g.
  [`screen_associations()`](https://rasmussenphilip.github.io/deconflate/reference/screen_associations.md))
  also run without any.

- ...:

  Not used (arguments of earlier versions give a message).

## Value

A `cm_population` object.

## Examples

``` r
pop <- cm_population(
  cm_diseases(c("d1", "d2", "d3"), c(0.10, 0.15, 0.20)),
  cm_associations(c("d1", "d2"), c("d2", "d3"), c(2, 3))
)
pop
#> <cm_population>
#>   Diseases: 3 (d1, d2, d3)
#>   Association rows: 2 (OR: 2)
#>   Disease pairs: 3 (2 with a crude pair row, 1 unknown)
```
