# Fetch change log information for one or more neurons

Fetch change log information for one or more neurons

## Usage

``` r
banc_change_log(x, tz = "UTC", filtered = TRUE, OmitFailures = TRUE, ...)
```

## Arguments

- x:

  One or more banc ids in any format understandable by
  [`ngl_segments`](https://rdrr.io/pkg/fafbseg/man/ngl_segments.html)

- tz:

  Time zone for edit timestamps. Defaults to "UTC" i.e. Universal Time,
  Coordinated. Set to "" for your current timezone. See
  [`as.POSIXct`](https://rdrr.io/r/base/as.POSIXlt.html) for more
  details.

- filtered:

  Whether to filter out edits unlikely to relate to the current state of
  the neuron (default `TRUE`, see details).

- OmitFailures:

  Whether to omit neurons for which there is an API timeout or error.
  The default value (`TRUE`) will skip over errors, while `NA`) will
  result in a hard stop on error. See
  [`nlapply`](https://rdrr.io/pkg/nat/man/nlapply.html) for more
  details.

- ...:

  Additional arguments passed to
  [`flywire_fetch`](https://rdrr.io/pkg/fafbseg/man/flywire_fetch.html)

## Value

a `data.frame` See
`fabseg::`[`flywire_change_log`](https://rdrr.io/pkg/fafbseg/man/flywire_change_log.html)
for details

## Details

As of August 2021 this is a simple wrapper of
`fafbseg::`[`flywire_change_log`](https://rdrr.io/pkg/fafbseg/man/flywire_change_log.html).
For now the old (and less convenient format) available from the zetta
API can be obtained with the private `bancr:::banc_change_log_zetta`
function.

## Examples

``` r
if (FALSE) { # \dontrun{
banc_change_log("720575941477428704")
} # }
```
