# Describe associations between diseases

Each row is one estimate from a study: an association between two
diseases, or between combinations of diseases, possibly within a
subgroup or adjusted for other diseases. Several rows can describe the
same association (they are pooled). Pairs that no row mentions are
unknown: the model fills in their association from the others (see
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)).
To state that two diseases are unrelated, give the pair an odds ratio of
1.

## Usage

``` r
cm_associations(
  disease,
  associated_with = NA_character_,
  value,
  measure = "OR",
  ...,
  among = NA_character_,
  adjusted_for = NA_character_,
  study = NA_character_,
  source = NA_character_,
  note = NA_character_
)
```

## Arguments

- disease, associated_with:

  Disease ids or combinations (`"a;b"`).

- value:

  Numeric association values (required).

- measure:

  Character, one per row or recycled (see What a row means).

- ...:

  Not used (arguments of earlier versions give a message). The arguments
  after it must be named.

- among:

  Optional subgroups (see What a row means).

- adjusted_for:

  Optional adjustment sets (see What a row means).

- study, source, note:

  Optional study label (used in diagnostics), citation and free-text
  note.

## Value

A `cm_associations` data frame (with the row number of each row in
column `row`).

## What a row means

`disease` (G) is the outcome and `associated_with` (H) the comparison:
animals with all of H against animals with none of H (for one disease,
with against without). Each can be one disease or a combination, with
ids separated by `";"` (e.g. `"MAS;SCK"`). `among` narrows the
population to a subgroup: ids that are present (`"SCK"`), absent
(`"no SCK"`), both (`"SCK;no MET"`), or `"no others"` (none of the
diseases not named in the row).

- `"OR"` (default), `"RR"`, `"RD"`, `"phi"`: odds ratio, risk ratio
  `P(G | H) / P(G | no H)`, risk difference and phi correlation of the
  2x2 table of G by H, within `among`.

- `"prevalence"`: `P(G | H)` within `among`; without `associated_with`,
  the prevalence of G (a combination, or a disease within a subgroup).
  The prevalence of one disease in the whole population belongs in
  [`cm_diseases()`](https://rasmussenphilip.github.io/deconflate/reference/cm_diseases.md).

- `"OR_ratio"`: the odds ratio of G and H among `among` divided by their
  odds ratio among the complement of `among` (one disease, present or
  absent): how much a third disease changes an association, as a study
  that reports stratum-specific odds ratios would give it.

`adjusted_for`: an odds ratio adjusted for other diseases of the model
(ids separated by `";"`, or `"all"` for every other disease), from a
logistic regression of G on H and those diseases, is an exact
conditional association and is used as such (G and H must be single
diseases). Names that are not diseases of the model (e.g. `"parity"`)
are factors outside the model: they are ignored, with a note, and a row
adjusted only for them is used as a crude association.

Measures are applied to the modelled population's own probabilities
("transported"): this assumes the measure is the same in the source
study and in the modelled population.

## Examples

``` r
cm_associations(c("d1", "d2"), c("d2", "d3"), c(2, 3))
#>   row disease associated_with among adjusted_for measure value study source
#> 1   1      d1              d2  <NA>         <NA>      OR     2  <NA>   <NA>
#> 2   2      d2              d3  <NA>         <NA>      OR     3  <NA>   <NA>
#>   note
#> 1 <NA>
#> 2 <NA>
# An odds ratio within a subgroup, and how much d3 changes it
cm_associations(c("d1", "d1"), c("d2", "d2"), c(2.5, 1.8),
                measure = c("OR", "OR_ratio"), among = c("d3", "d3"))
#>   row disease associated_with among adjusted_for  measure value study source
#> 1   1      d1              d2    d3         <NA>       OR   2.5  <NA>   <NA>
#> 2   2      d1              d2    d3         <NA> OR_ratio   1.8  <NA>   <NA>
#>   note
#> 1 <NA>
#> 2 <NA>
```
