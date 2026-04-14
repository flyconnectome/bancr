# Summarise the connectivity of BANC neurons

Returns synaptically connected partners for specified neurons.
Understanding synaptic partnerships is crucial for analyzing neural
circuits in the Brain And Nerve Cord (BANC) connectome, revealing how
distributed control architecture coordinates behaviour across brain and
ventral nerve cord regions.

`banc_partners` returns details of each unitary synaptic connection
(including its xyz location).

## Usage

``` r
banc_partner_summary(
  rootids,
  partners = c("outputs", "inputs"),
  synapse_table = c("synapses_v2", "synapses_v3", "synapses_v1"),
  threshold = 0,
  remove_autapses = TRUE,
  cleft.threshold = 0,
  datastack_name = NULL,
  ...
)

banc_partners(
  rootids,
  partners = c("input", "output"),
  synapse_table = c("synapses_v2", "synapses_v3", "synapses_v1"),
  details = FALSE,
  ...
)
```

## Arguments

- rootids:

  Character vector specifying one or more BANC rootids. As a convenience
  this argument is passed to
  [`banc_ids`](https://flyconnectome.github.io/bancr/reference/banc_ids.md)
  allowing you to pass in data.frames, BANC URLs or simple ids.

- partners:

  Character vector, either "outputs" or "inputs" to specify the
  direction of synaptic connections to retrieve.

- synapse_table:

  Character, the name of the synapse CAVE table you wish to use.
  Defaults to the latest.

- threshold:

  Integer threshold for minimum number of synapses (default 0).

- remove_autapses:

  Logical, whether to remove self-connections (default TRUE).

- cleft.threshold:

  Numeric threshold for cleft filtering (default 0).

- datastack_name:

  An optional CAVE `datastack_name`. If unset a sensible default is
  chosen.

- ...:

  Additional arguments passed to
  [`flywire_partner_summary`](https://rdrr.io/pkg/fafbseg/man/flywire_partners.html)

- details:

  Logical. If `TRUE` and `synapse_table="synapses_v3"`, additionally
  fetch `mean_score` and `median_score` from the reference tables
  `synapses_v3_mean_score` and `synapses_v3_median_score` and merge them
  onto the returned data.frame by synapse id. Default `FALSE`. Note: the
  median-score join is substantially slower than the mean-score join (on
  the order of 10x), so only request `details=TRUE` when you need the
  per-synapse scores. Silently ignored for v2/v1.

## Value

a data.frame

## Details

note that the rootids you pass in must be up to date. See example.

## See also

[`flywire_partner_summary`](https://rdrr.io/pkg/fafbseg/man/flywire_partners.html),
[`banc_latestid`](https://flyconnectome.github.io/bancr/reference/banc_latestid.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic connectivity analysis
sample_id=banc_latestid("720575941478275714")
head(banc_partner_summary(sample_id))
head(banc_partner_summary(sample_id, partners='inputs'))

# Research application: Analyze descending neuron control circuits
library(dplyr)

# Get DNa02 descending neurons that control walking behavior
dna02_annotations <- banc_codex_annotations() %>%
  filter(cell_type == "DNa02")
dna02_id <- dna02_annotations$pt_root_id[1]

# Find their downstream targets in the VNC
dna02_outputs <- banc_partner_summary(dna02_id, partners='outputs') %>%
  slice_max(weight, n = 10)

# Visualize the circuit in neuroglancer
banc_partner_summary(sample_id, partners='inputs') %>%
  slice_max(weight, n = 20) %>%
  banc_scene(open=TRUE)
} # }
if (FALSE) { # \dontrun{
# plot input and output synapses of a neuron
nclear3d()
fpi=banc_partners(banc_latestid("720575941478275714"), partners='in')
points3d(banc_raw2nm(fpi$post_pt_position), col='cyan')
fpo=banc_partners(banc_latestid("720575941478275714"), partners='out')
points3d(banc_raw2nm(fpo$pre_pt_position), col='red')

# Compare results between the v2 (default) and v3 synapse tables
id <- banc_latestid("720575941478275714")
fpi_v2 <- banc_partners(id, partners='input', synapse_table="synapses_v2")
fpi_v3 <- banc_partners(id, partners='input', synapse_table="synapses_v3")
nrow(fpi_v2); nrow(fpi_v3)
# partner overlap
length(intersect(fpi_v2$pre_pt_root_id, fpi_v3$pre_pt_root_id))

# Pull v3 synapses with mean and median scores attached (slower)
fpi_v3d <- banc_partners(id, partners='input', synapse_table="synapses_v3",
                         details=TRUE)
head(fpi_v3d[, c("id", "mean_score", "median_score")])
} # }
```
