# Set or replace one three-way association

Sets how much `disease3` changes the association of `disease1` and
`disease2`: an `OR_ratio` row (the odds ratio of disease1 and disease2
among animals with disease3, divided by that among animals without it;
see
[`cm_associations()`](https://rasmussenphilip.github.io/deconflate/reference/cm_associations.md)).
Replaces such a row if present.

## Usage

``` r
set_three_way(model, disease1, disease2, disease3, ratio)
```

## Arguments

- model:

  A
  [`cm_population()`](https://rasmussenphilip.github.io/deconflate/reference/cm_population.md)
  or
  [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md).

- disease1, disease2, disease3:

  Disease ids.

- ratio:

  Positive ratio of the two odds ratios.

## Value

The modified object. Distributions of replaced rows are dropped.
