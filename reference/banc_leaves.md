# Find the supervoxel identifiers of a banc neuron

Find the supervoxel identifiers of a banc neuron

## Usage

``` r
banc_leaves(x, integer64 = TRUE, ...)
```

## Arguments

- x:

  One or more FlyWire segment ids

- integer64:

  Whether to return ids as integer64 type (the default, more compact but
  a little fragile) rather than character (when `FALSE`).

- ...:

  additional arguments passed to
  [`flywire_leaves`](https://rdrr.io/pkg/fafbseg/man/flywire_leaves.html)

## Value

A vector of supervoxel ids

## See also

[`flywire_leaves`](https://rdrr.io/pkg/fafbseg/man/flywire_leaves.html)

Other banc-ids:
[`banc_cellid_from_segid()`](https://flyconnectome.github.io/bancr/reference/banc_cellid_from_segid.md),
[`banc_ids()`](https://flyconnectome.github.io/bancr/reference/banc_ids.md),
[`banc_islatest()`](https://flyconnectome.github.io/bancr/reference/banc_islatest.md),
[`banc_latestid()`](https://flyconnectome.github.io/bancr/reference/banc_latestid.md),
[`banc_rootid()`](https://flyconnectome.github.io/bancr/reference/banc_rootid.md),
[`banc_xyz2id()`](https://flyconnectome.github.io/bancr/reference/banc_xyz2id.md)

## Examples

``` r
if (FALSE) { # \dontrun{
svids=banc_leaves("720575941478275714")
head(svids)
} # }
```
