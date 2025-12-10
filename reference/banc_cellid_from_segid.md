# Convert between BANC cell ids and root ids

Converts between BANC cell ids (should survive most edits) and root ids
(guaranteed to match just one edit state). See details.

## Usage

``` r
banc_cellid_from_segid(
  rootids = NULL,
  timestamp = NULL,
  version = NULL,
  cellid_table = NULL,
  rval = c("ids", "data.frame")
)

banc_segid_from_cellid(
  cellids = NULL,
  timestamp = NULL,
  version = NULL,
  rval = c("ids", "data.frame"),
  integer64 = FALSE,
  cellid_table = NULL
)
```

## Arguments

- rootids:

  banc root ids in any form understood by
  [`banc_ids`](https://flyconnectome.github.io/bancr/reference/banc_ids.md).
  The default value of NULL will return all cell ids.

- timestamp:

  An optional time stamp. You should give only one of `version` or
  `timestamp`. When both are missing, ids should match the live
  materialisation version including up to the second edits.

- version:

  An optional integer CAVE materialisation version. You should give only
  one of `version` or `timestamp`. When both are missing, ids should
  match the live materialisation version including up to the second
  edits.

- cellid_table:

  Optional name of cell id table (the default value of `NULL` should
  find the correct table).

- rval:

  Whether to return the cell ids or the whole of the CAVE table with
  additional columns.

- cellids:

  Integer cell ids between between 1 and around 20000 that *should*
  uniquely identify each cell in the dataset.

- integer64:

  Whether to return ids as
  [`bit64::integer64`](https://rdrr.io/pkg/bit64/man/bit64-package.html)
  or character vectors. Default value of NA leaves the ids unmodified.

## Value

Either a vector of ids or a data.frame depending on `rval`. For cell ids
the vector will be an integer for root ids (segment ids), a character
vector or an
[`bit64::integer64`](https://rdrr.io/pkg/bit64/man/bit64-package.html)
vector depending on the `integer64` argument.

## Details

CAVE/PyChunkedGraph assigns a 64 bit integer root id to all bodies in
the segmentation. These root ids are persistent in a computer science
sense, which is often the exact opposite of what neuroscientists might
imagine. Specifically, a given root id is matched to a single edit state
of a neuron. If the neuron is edited, then root id changes. In contrast,
cell ids do not change even in the face of edits. However, it is
important to understand that they correspond to a specific point on a
neuron, commonly the nucleus. If the nucleus is edited away from a the
rest of a neuron to which it previously belonged, then the cell id and
any associated edits will effectively with move it.

For further details see [banc
slack](https://banc-reconstruction.slack.com/archives/CLDH21J4U/p1690755500802509)
and [banc
wiki](https://github.com/htem/banc_auto_recon/wiki/Neuron-annotations#neuron_information).

## See also

Other banc-ids:
[`banc_ids()`](https://flyconnectome.github.io/bancr/reference/banc_ids.md),
[`banc_islatest()`](https://flyconnectome.github.io/bancr/reference/banc_islatest.md),
[`banc_latestid()`](https://flyconnectome.github.io/bancr/reference/banc_latestid.md),
[`banc_leaves()`](https://flyconnectome.github.io/bancr/reference/banc_leaves.md),
[`banc_rootid()`](https://flyconnectome.github.io/bancr/reference/banc_rootid.md),
[`banc_xyz2id()`](https://flyconnectome.github.io/bancr/reference/banc_xyz2id.md)

## Examples

``` r
if (FALSE) { # \dontrun{
 banc_cellid_from_segid(banc_latestid("720575941626035769"))
} # }
if (FALSE) { # \dontrun{
banc_cellid_from_segid(banc_latestid("720575941480769421"))
} # }
```
