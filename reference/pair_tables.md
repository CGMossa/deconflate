# Pairwise 2x2 tables for every disease pair

The tables that the crude pair rows of the associations give (one
disease in `disease`, one in `associated_with`, no `among`, not adjusted
for diseases of the model). Several rows of one pair are pooled (see
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)).

## Usage

``` r
pair_tables(model, pool = c("random", "fixed"))
```

## Arguments

- model:

  A
  [`cm_population()`](https://rasmussenphilip.github.io/deconflate/reference/cm_population.md)
  or
  [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md).

- pool:

  `"random"` (default) or `"fixed"`: how several rows of one pair are
  pooled.

## Value

A data frame with one row per unordered pair: the measure (or
`"pooled"`), its `status` (`"specified"`: a crude pair row;
`"other rows"`: only other kinds of rows mention the pair, e.g. within a
subgroup; `"unknown"`: no row mentions it), the association rows used,
the joint probability `p11` (`NA` unless specified), the implied odds
ratio, and the excess probabilities
`ep_2_given_1 = P(d2 | d1) - P(d2 | not d1)` and `ep_1_given_2`.
