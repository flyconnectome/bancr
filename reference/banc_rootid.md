# Find the root identifier of a banc neuron

Find the root identifier of a banc neuron

## Usage

``` r
banc_rootid(x, integer64 = FALSE, version = NULL, timestamp = NULL, ...)
```

## Arguments

- x:

  One or more FlyWire segment ids

- integer64:

  Whether to return ids as integer64 type (more compact but a little
  fragile) rather than character (default `FALSE`).

- version:

  An optional CAVE materialisation version number. See details and
  examples.

- timestamp:

  An optional timestamp as a string or POSIXct, interpreted as UTC when
  no timezone is specified.

- ...:

  Additional arguments passed to the underlying functions and eventually
  to Python `cv$CloudVolume` object.

## Value

A vector of root ids (by default character)

## See also

[`flywire_rootid`](https://rdrr.io/pkg/fafbseg/man/flywire_rootid.html)

Other banc-ids:
[`banc_cellid_from_segid()`](https://flyconnectome.github.io/bancr/reference/banc_cellid_from_segid.md),
[`banc_ids()`](https://flyconnectome.github.io/bancr/reference/banc_ids.md),
[`banc_islatest()`](https://flyconnectome.github.io/bancr/reference/banc_islatest.md),
[`banc_latestid()`](https://flyconnectome.github.io/bancr/reference/banc_latestid.md),
[`banc_leaves()`](https://flyconnectome.github.io/bancr/reference/banc_leaves.md),
[`banc_xyz2id()`](https://flyconnectome.github.io/bancr/reference/banc_xyz2id.md)

## Examples

``` r
if (FALSE) { # \dontrun{
banc_rootid("73186243730767724")
} # }
```
