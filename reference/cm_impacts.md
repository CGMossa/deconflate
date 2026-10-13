# Describe raw impact estimates

One impact table holds the estimates of one outcome, one row per
estimate. A disease can have several rows (e.g. from several studies, or
crude and adjusted), and a row can be about a combination of diseases or
a subgroup. There are two kinds of impact table:

## Usage

``` r
cm_impacts(
  disease,
  value,
  ...,
  among = NA_character_,
  adjusted_for = NA_character_,
  measure = NULL,
  study = NA_character_,
  source = NA_character_,
  note = NA_character_,
  label = NULL,
  units = NULL
)
```

## Arguments

- disease:

  Character vector of disease ids or combinations (`"a;b"`).

- value:

  Numeric raw impacts.

- ...:

  Not used (arguments of earlier versions give a message). The arguments
  after it must be named.

- among, adjusted_for:

  Optional subgroups and adjustment sets (see What a row means).

- measure:

  `NULL` (default) for additive impacts; for event impacts, one measure
  per row or recycled.

- study, source, note:

  Optional study label, citation and note.

- label, units:

  Optional table-level metadata (e.g. `label = "milk yield loss"`,
  `units = "% of yield"`), carried into the results.

## Value

A `cm_impacts` data frame (attributes `label`, `units` and `kind`,
`"additive"` or `"event"`).

## Details

- **Additive impacts** (no `measure`): amounts in the outcome's own
  units, the same for every row (e.g. kg of milk, percent of yield,
  days, a welfare score). The package does not convert units; results
  come back in the units supplied. Adjust them with
  [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md).

- **Event impacts** (`measure` given): comparisons of the risk of an
  event (e.g. death or culling), as hazard ratios (`"HR"`), rate ratios
  (`"rate_ratio"`), risk ratios (`"RR"`), odds ratios (`"OR"`) or risk
  differences (`"RD"`), which can be mixed. Adjust them with
  `deconflate(..., event_model = TRUE)` (snapshot hazard model).

## What a row means

The model has an impact for each disease and an interaction for each
combination that appears in `disease`. A row compares animals with all
the diseases in `disease` against animals with none of them (for one
disease: with against without), within the subgroup `among` (as in
[`cm_associations()`](https://rasmussenphilip.github.io/deconflate/reference/cm_associations.md),
including `"no others"`, e.g. "only this disease against healthy
animals").

- Crude (empty `adjusted_for`): the difference in the outcome between
  the two groups (event impacts: the ratio or difference of their
  hazards or risks).

- Adjusted (`adjusted_for` names diseases, separated by `";"`): additive
  impacts are coefficients of a linear regression of the outcome on the
  disease (for a combination: its members, their smaller combinations
  and the combination itself) and the diseases in `adjusted_for`; event
  impacts are comparisons within strata of those diseases (for a
  combination, the interaction contrast), combined with
  Mantel-Haenszel-type weights. `"all"` means every other disease and
  combination term of the model: the row is then the term's own value (a
  combination adjusted for `all` is its interaction; for event impacts
  this holds for hazard and rate ratios). Within a subgroup `among`
  naming present diseases, it is the sum of the term and the terms that
  combine it with those diseases (e.g. `a` among `c`, adjusted for all,
  is the impact of `a` plus the `a;c` interaction, when the model has
  one). Names that are not diseases of the model are factors outside the
  model, ignored with a note.

Every disease needs at least one row of its own (for no effect: 0 for
additive impacts and risk differences, 1 for ratios). The rows must
determine every impact and interaction;
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
says which are not determined otherwise. Several rows about the same
impacts are fitted together by weighted least squares (see
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)).

A Cox hazard ratio estimated over follow-up is not exactly a snapshot
comparison (see
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md),
"Event impacts"); entering one is an approximation and the user's
assumption. Risk ratios, odds ratios and risk differences refer to the
period of the overall risk given to
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md).

## Examples

``` r
cm_impacts(c("d1", "d2", "d3"), c(2.5, 5, 7.5), label = "yield", units = "%")
#>   row disease among adjusted_for measure value study source note
#> 1   1      d1  <NA>         <NA>    <NA>   2.5  <NA>   <NA> <NA>
#> 2   2      d2  <NA>         <NA>    <NA>   5.0  <NA>   <NA> <NA>
#> 3   3      d3  <NA>         <NA>    <NA>   7.5  <NA>   <NA> <NA>
# Two estimates for d1, and the interaction of d1 and d2
cm_impacts(c("d1", "d1", "d2", "d3", "d1;d2"), c(2.5, 2.1, 5, 7.5, 1),
           adjusted_for = c(NA, "d2", NA, NA, "all"), units = "%")
#>   row disease among adjusted_for measure value study source note
#> 1   1      d1  <NA>         <NA>    <NA>   2.5  <NA>   <NA> <NA>
#> 2   2      d1  <NA>           d2    <NA>   2.1  <NA>   <NA> <NA>
#> 3   3      d2  <NA>         <NA>    <NA>   5.0  <NA>   <NA> <NA>
#> 4   4      d3  <NA>         <NA>    <NA>   7.5  <NA>   <NA> <NA>
#> 5   5   d1;d2  <NA>          all    <NA>   1.0  <NA>   <NA> <NA>
# Event impacts: culling, with mixed measures
cm_impacts(c("d1", "d2", "d3"), c(1.5, 1.3, 0.04), measure = c("HR", "RR", "RD"),
           label = "culling")
#>   row disease among adjusted_for measure value study source note
#> 1   1      d1  <NA>         <NA>      HR  1.50  <NA>   <NA> <NA>
#> 2   2      d2  <NA>         <NA>      RR  1.30  <NA>   <NA> <NA>
#> 3   3      d3  <NA>         <NA>      RD  0.04  <NA>   <NA> <NA>
```
