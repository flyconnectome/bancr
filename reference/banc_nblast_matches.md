# Read NBLAST match results from CAVE

Read NBLAST match results from CAVE

## Usage

``` r
banc_nblast_matches(
  dataset = c("malecns", "fafb", "hemibrain", "manc", "fanc"),
  ...
)
```

## Arguments

- dataset:

  Character, which cross-species NBLAST comparison to query.

- ...:

  Additional arguments passed to
  [`banc_cave_query`](https://flyconnectome.github.io/bancr/reference/banc_cave_query.md)

## Value

A data.frame with NBLAST match results (cell_match schema)

## Examples

``` r
if (FALSE) { # \dontrun{
# Get all maleCNS NBLAST matches
matches <- banc_nblast_matches("malecns")

# Get validated matches only
validated <- banc_nblast_matches("fafb") %>%
  dplyr::filter(validation == TRUE)
} # }
```
