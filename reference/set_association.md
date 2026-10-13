# Set or replace one association

Replaces the crude pair rows of a pair (one disease on each side, no
subgroup, not adjusted) by one row with the given value, or adds that
row. Other rows about the pair (e.g. within subgroups) are kept.

## Usage

``` r
set_association(model, disease1, disease2, value, measure = "OR")
```

## Arguments

- model:

  A
  [`cm_population()`](https://rasmussenphilip.github.io/deconflate/reference/cm_population.md)
  or
  [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md).

- disease1, disease2:

  Disease ids.

- value:

  Association value.

- measure:

  Measure (see
  [`cm_associations()`](https://rasmussenphilip.github.io/deconflate/reference/cm_associations.md)).

## Value

The modified object. Distributions of replaced rows are dropped (the
scenario value is fixed).
