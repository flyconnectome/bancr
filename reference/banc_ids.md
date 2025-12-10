# Return a vector of banc root ids from diverse inputs

Return a vector of banc root ids from diverse inputs

## Usage

``` r
banc_ids(x, integer64 = NA)
```

## Arguments

- x:

  A data.frame, URL or vector of ids

- integer64:

  Whether to return ids as
  [`bit64::integer64`](https://rdrr.io/pkg/bit64/man/bit64-package.html)
  or character vectors. Default value of NA leaves the ids unmodified.

## Value

A vector of ids

## See also

Other banc-ids:
[`banc_cellid_from_segid()`](https://flyconnectome.github.io/bancr/reference/banc_cellid_from_segid.md),
[`banc_islatest()`](https://flyconnectome.github.io/bancr/reference/banc_islatest.md),
[`banc_latestid()`](https://flyconnectome.github.io/bancr/reference/banc_latestid.md),
[`banc_leaves()`](https://flyconnectome.github.io/bancr/reference/banc_leaves.md),
[`banc_rootid()`](https://flyconnectome.github.io/bancr/reference/banc_rootid.md),
[`banc_xyz2id()`](https://flyconnectome.github.io/bancr/reference/banc_xyz2id.md)

## Examples

``` r
banc_ids(data.frame(rootid="648518346474360770"))
#> [1] "648518346474360770"
```
