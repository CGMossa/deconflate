# Set or replace one pairwise interaction

Sets the interaction of two diseases to `value` (in the units of the
impacts): a row for the combination adjusted for `all`, which is the
interaction itself (see
[`cm_impacts()`](https://rasmussenphilip.github.io/deconflate/reference/cm_impacts.md)).
Replaces such a row if present.

## Usage

``` r
set_interaction(model, disease1, disease2, value)
```

## Arguments

- model:

  A
  [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md).

- disease1, disease2:

  Disease ids.

- value:

  Interaction value (same units as the impacts).

## Value

The modified model. Distributions of replaced rows are dropped.
