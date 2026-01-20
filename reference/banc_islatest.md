# Check if a banc root id is up to date

Check if a banc root id is up to date

## Usage

``` r
banc_islatest(x, timestamp = NULL, ...)
```

## Arguments

- x:

  FlyWire rootids in any format understandable to
  [`ngl_segments`](https://rdrr.io/pkg/fafbseg/man/ngl_segments.html)
  including as `integer64`

- timestamp:

  (optional) argument to set an endpoint - edits after this time will be
  ignored (see details).

- ...:

  Additional arguments passed to
  [`flywire_islatest`](https://rdrr.io/pkg/fafbseg/man/flywire_islatest.html)

## See also

Other banc-ids:
[`banc_cellid_from_segid()`](https://flyconnectome.github.io/bancr/reference/banc_cellid_from_segid.md),
[`banc_ids()`](https://flyconnectome.github.io/bancr/reference/banc_ids.md),
[`banc_latestid()`](https://flyconnectome.github.io/bancr/reference/banc_latestid.md),
[`banc_leaves()`](https://flyconnectome.github.io/bancr/reference/banc_leaves.md),
[`banc_rootid()`](https://flyconnectome.github.io/bancr/reference/banc_rootid.md),
[`banc_xyz2id()`](https://flyconnectome.github.io/bancr/reference/banc_xyz2id.md)

## Examples

``` r
banc_islatest("720575941520182775")
#> Warning: running command ''/home/runner/.cache/R/reticulate/uv/cache/archive-v0/rNCHVmTPN6nWgKMm8Xu-Y/bin/python' -m pip freeze' had status 1
#> Warning: using default setting
#> [1] FALSE
```
