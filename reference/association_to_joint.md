# Joint probability of a disease pair from any supported measure

Joint probability of a disease pair from any supported measure

## Usage

``` r
association_to_joint(measure, value, p1, p2)
```

## Arguments

- measure:

  `"OR"`, `"RR"`, `"RD"`, `"phi"` or `"prevalence"` (see
  [`cm_associations()`](https://rasmussenphilip.github.io/deconflate/reference/cm_associations.md);
  the directional measures treat disease 1 as the outcome and disease 2
  as the comparison).

- value:

  Association value.

- p1, p2:

  Marginal probabilities of disease 1 and disease 2.

## Value

`P(d1 and d2)`. Errors with class `deconflate_infeasible` if the measure
is incompatible with the marginals.
