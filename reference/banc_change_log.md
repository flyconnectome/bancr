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
# \donttest{
banc_change_log("720575941477428704")
#> Warning: running command ''/home/runner/.cache/R/reticulate/uv/cache/archive-v0/eKU3DhU7NupYnOukbdetF/bin/python' -m pip freeze' had status 1
#> Warning: using default setting
#>    operation_id           timestamp user_id
#> 0         10515 2023-11-28 22:20:43    2660
#> 1        127798 2023-12-20 06:15:55    2767
#> 2        144911 2023-12-26 00:28:30    4162
#> 3        147922 2023-12-26 06:04:13    2758
#> 4        147926 2023-12-26 06:05:40    2758
#> 5        147945 2023-12-26 06:06:56    2758
#> 6        147980 2023-12-26 06:17:43    2758
#> 7        147981 2023-12-26 06:17:44    2758
#> 8        147982 2023-12-26 06:17:47    2758
#> 9        147983 2023-12-26 06:17:47    2758
#> 10       147984 2023-12-26 06:17:50    2758
#> 11       147985 2023-12-26 06:17:51    2758
#> 12       147987 2023-12-26 06:17:56    2758
#> 13       147989 2023-12-26 06:18:02    2758
#> 14       147990 2023-12-26 06:18:27    2758
#> 15       147993 2023-12-26 06:18:32    2758
#> 16       148000 2023-12-26 06:20:22    2758
#>                          before_root_ids     after_root_ids is_merge
#> 0                     720575941480769421 720575941556888079    FALSE
#> 1                     720575941556888079 720575941474194206    FALSE
#> 2                     720575941474194206 720575941435522834    FALSE
#> 3                     720575941435522834 720575941488179406    FALSE
#> 4                     720575941488179406 720575941539449692    FALSE
#> 5                     720575941539449692 720575941557292526    FALSE
#> 6  720575941551446830 720575941557292526 720575941503347109     TRUE
#> 7  720575940543311189 720575941489220228 720575941483891974     TRUE
#> 8  720575941483891974 720575941503347109 720575941477813667     TRUE
#> 9  720575941125943612 720575941555221066 720575941474604865     TRUE
#> 10 720575941039420897 720575941474604865 720575941447470414     TRUE
#> 11 720575941477813667 720575941558385512 720575941515845036     TRUE
#> 12 720575941447470414 720575941515845036 720575941547060285     TRUE
#> 13 720575941547060285 720575941569217576 720575941504051269     TRUE
#> 14 720575941454232760 720575941504051269 720575941636113781     TRUE
#> 15 720575941513505412 720575941636113781 720575941636114037     TRUE
#> 16 720575941352474137 720575941636114037 720575941477428704     TRUE
#>                   user_name                     user_affiliation
#> 0             Jasper Phelps                    Wei-Chung Lee lab
#> 1           Allien Mae Gogo Mala Murthy Lab, Sebastian Seung Lab
#> 2  David Benjamin Conquilla Mala Murthy Lab, Sebastian Seung Lab
#> 3           Jacquilyn Laude Mala Murthy Lab, Sebastian Seung Lab
#> 4           Jacquilyn Laude Mala Murthy Lab, Sebastian Seung Lab
#> 5           Jacquilyn Laude Mala Murthy Lab, Sebastian Seung Lab
#> 6           Jacquilyn Laude Mala Murthy Lab, Sebastian Seung Lab
#> 7           Jacquilyn Laude Mala Murthy Lab, Sebastian Seung Lab
#> 8           Jacquilyn Laude Mala Murthy Lab, Sebastian Seung Lab
#> 9           Jacquilyn Laude Mala Murthy Lab, Sebastian Seung Lab
#> 10          Jacquilyn Laude Mala Murthy Lab, Sebastian Seung Lab
#> 11          Jacquilyn Laude Mala Murthy Lab, Sebastian Seung Lab
#> 12          Jacquilyn Laude Mala Murthy Lab, Sebastian Seung Lab
#> 13          Jacquilyn Laude Mala Murthy Lab, Sebastian Seung Lab
#> 14          Jacquilyn Laude Mala Murthy Lab, Sebastian Seung Lab
#> 15          Jacquilyn Laude Mala Murthy Lab, Sebastian Seung Lab
#> 16          Jacquilyn Laude Mala Murthy Lab, Sebastian Seung Lab
# }
```
