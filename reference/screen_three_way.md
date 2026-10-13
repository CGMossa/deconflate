# Screen three-way association scenarios

Pairwise associations do not identify how often three diseases occur
together; the model assumes as little as possible (maximum entropy).
This screen sets, for each triple in turn, how much the third disease
changes the association of the first two (an `OR_ratio` row, see
[`set_three_way()`](https://rasmussenphilip.github.io/deconflate/reference/set_three_way.md)),
refits the joint distribution (the other rows are still matched) and
reports the change in the aggregate and rankings.

## Usage

``` r
screen_three_way(
  model,
  ratios = c(0.5, 2),
  triples = NULL,
  event_model = FALSE,
  overall_risk = NULL
)
```

## Arguments

- model:

  A
  [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md).

- ratios:

  Odds ratio ratios to try.

- triples:

  Optional list of character vectors of three disease ids (the
  association of the first two, among the third); default all triples
  (which can be many).

- event_model, overall_risk:

  As in
  [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md).

## Value

A `cm_screen` data frame as in
[`screen_associations()`](https://rasmussenphilip.github.io/deconflate/reference/screen_associations.md).

## Details

When the pairs of a triple all have crude pair rows, such a row keeps
their tables fixed, so it changes additive results only through
interactions or rows about combinations or subgroups. When a pair is
unknown, its fitted table changes, and additive results can change even
without interactions. Three-way associations also matter for event
impacts (`event_model = TRUE`).
