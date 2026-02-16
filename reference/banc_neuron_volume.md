# Calculate neuron volume from CAVE L2 cache

Computes the total volume (in nm³) for one or more BANC neurons by
summing pre-computed `size_nm3` values from the CAVE L2 cache. This is
much faster than downloading full meshes.

## Usage

``` r
banc_neuron_volume(ids, OmitFailures = TRUE, ...)
```

## Arguments

- ids:

  A vector of one or more neuron root IDs.

- OmitFailures:

  Logical; if `TRUE`, neurons that fail are returned as `NA` rather than
  raising an error (default `TRUE`).

- ...:

  Additional arguments (currently unused).

## Value

A `data.frame` with columns:

- root_id:

  Character root ID

- volume_nm3:

  Total volume in cubic nanometers

## Details

For each neuron, `banc_neuron_volume` retrieves all Level 2 chunk IDs
via the chunkedgraph, then queries the L2 cache for the `size_nm3`
attribute of each chunk. The total neuron volume is the sum across all
chunks.

The L2 cache stores pre-computed statistics that are updated after every
CAVE edit, so results reflect the current segmentation state.

## See also

[`banc_read_l2skel`](https://flyconnectome.github.io/bancr/reference/banc_read_l2dp.md)
for L2 skeleton data

## Examples

``` r
if (FALSE) { # \dontrun{
# Single neuron
banc_neuron_volume("720575941478275714")

# Multiple neurons
ids <- c("720575941478275714", "720575941512951014")
banc_neuron_volume(ids)
} # }
```
