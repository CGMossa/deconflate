# Describe the diseases in a system

One row per disease: how common it is in the whole population. Rows
about two or more diseases (associations, combinations, prevalences
within subgroups) belong in
[`cm_associations()`](https://rasmussenphilip.github.io/deconflate/reference/cm_associations.md).

## Usage

``` r
cm_diseases(
  disease,
  value,
  measure = "prevalence",
  label = NA_character_,
  time_horizon = NA_character_,
  reference_population = NA_character_,
  source = NA_character_,
  note = NA_character_,
  ...
)
```

## Arguments

- disease:

  Character vector of unique disease ids (no `|`, `;` or `:`; `"all"`
  and `"others"` are reserved).

- value:

  Numeric vector of disease occurrence values, interpreted according to
  `measure`.

- measure:

  Character, one per disease (or recycled):

  - `"prevalence"`: proportion of animals affected at a representative
    point in time (used as the probability directly).

  - `"probability"`: probability of at least one case within
    `time_horizon` (used directly).

  - `"incidence_rate"`: expected number of cases per animal within
    `time_horizon` (e.g. `0.48` for 48 cases per 100 lactations).
    Converted to a probability assuming Poisson-distributed events,
    `1 - exp(-value)`, as in Rasmussen et al. (2024), eq. 1.

- label:

  Optional full names of the diseases (e.g. `"lameness"`).

- time_horizon:

  Character label of the period the probabilities refer to (e.g.
  `"lactation"`, `"year"`). Probabilities, associations and impacts
  should all refer to the same period.

- reference_population:

  Optional description of the population the values come from.

- source, note:

  Optional citation and free-text note.

- ...:

  Not used (arguments of earlier versions give a message).

## Value

A `cm_diseases` data frame with the converted probability in column
`prob` and the inputs retained.

## Examples

``` r
cm_diseases(c("SCK", "PTB"), c(0.4789, 0.1001),
            measure = c("incidence_rate", "prevalence"),
            label = c("subclinical ketosis", "paratuberculosis"),
            time_horizon = "lactation")
#>   disease      prob  value        measure               label time_horizon
#> 1     SCK 0.3805356 0.4789 incidence_rate subclinical ketosis    lactation
#> 2     PTB 0.1001000 0.1001     prevalence    paratuberculosis    lactation
#>   reference_population source note
#> 1                 <NA>   <NA> <NA>
#> 2                 <NA>   <NA> <NA>
```
