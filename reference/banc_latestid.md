# Find the latest id for a banc root id

Find the latest id for a banc root id

## Usage

``` r
banc_latestid(
  rootid,
  sample = 1000L,
  cloudvolume.url = NULL,
  Verbose = FALSE,
  ...
)

banc_updateids(
  x,
  root.column = "root_id",
  supervoxel.column = "supervoxel_id",
  position.column = "position",
  use.cave = TRUE,
  serial = FALSE,
  ...
)
```

## Arguments

- rootid:

  One ore more FlyWire rootids defining a segment (in any form
  interpretable by
  [`ngl_segments`](https://rdrr.io/pkg/fafbseg/man/ngl_segments.html))

- sample:

  An absolute or fractional number of supervoxel ids to map to rootids
  or `FALSE` (see details).

- cloudvolume.url:

  URL for CloudVolume to fetch segmentation image data. The default
  value of NULL chooses the flywire production segmentation dataset.

- Verbose:

  When set to `TRUE` prints information about what fraction of

- ...:

  Additional arguments passed to
  [`flywire_latestid`](https://rdrr.io/pkg/fafbseg/man/flywire_latestid.html)

- x:

  a `data.frame` with at least one of: `root_id`, `pt_root_id`,
  `supervoxel_id` and/or `pt_supervoxel_id`. Supervoxels will be
  preferentially used to update the `root_id` column. Else a vector of
  `BANC` root IDs.

- root.column:

  when `x` is a `data.frame`, the `root_id` column you wish to update

- supervoxel.column:

  when `x` is a `data.frame`, the `supervoxel_id` column you wish to use
  to update `root.column`

- position.column:

  when `x` is a `data.frame`, the `position` column with xyz values you
  wish to use to update `supervoxel.column`

- use.cave:

  read from the best established CAVE tables and join by
  `pt_supervoxel_id` to update `root_id`

- serial:

  if TRUE and x is a vector, calls `banc_updateids` on each ID in
  sequence to bufffer against connection failures. Slower.

## See also

[`banc_islatest`](https://flyconnectome.github.io/bancr/reference/banc_islatest.md)

Other banc-ids:
[`banc_cellid_from_segid()`](https://flyconnectome.github.io/bancr/reference/banc_cellid_from_segid.md),
[`banc_ids()`](https://flyconnectome.github.io/bancr/reference/banc_ids.md),
[`banc_islatest()`](https://flyconnectome.github.io/bancr/reference/banc_islatest.md),
[`banc_leaves()`](https://flyconnectome.github.io/bancr/reference/banc_leaves.md),
[`banc_rootid()`](https://flyconnectome.github.io/bancr/reference/banc_rootid.md),
[`banc_xyz2id()`](https://flyconnectome.github.io/bancr/reference/banc_xyz2id.md)

## Examples

``` r
if (FALSE) { # \dontrun{
banc_latestid("720575941520182775")
} # }
```
