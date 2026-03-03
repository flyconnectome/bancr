# Read NBLAST match results from CAVE

Query cross-species NBLAST match results stored in CAVE tables. Each
table contains pairwise morphological similarity scores between BANC
neurons and neurons from another connectome dataset, computed using
NBLAST (Costa et al., 2016). Matches are identified by running all BANC
neurons against a target dataset and retaining hits above a score
threshold.

## Usage

``` r
banc_nblast_matches(
  dataset = c("malecns", "fafb", "hemibrain", "manc", "fanc"),
  ...
)
```

## Arguments

- dataset:

  Character, which cross-species NBLAST comparison to query. One of:

  `"malecns"`

  :   CAVE table: `banc_malecns_nblast_v2`. Matches to the male CNS
      (Takemura et al., 2024) v0.9 dataset, covering the complete male
      central nervous system (~75K neurons).

  `"fafb"`

  :   CAVE table: `banc_fafb_nblast_v2`. Matches to FAFB (Zheng et al.,
      2018; Dorkenwald et al., 2024) FlyWire v783 dataset, a complete
      female *Drosophila* brain connectome.

  `"hemibrain"`

  :   CAVE table: `banc_hemibrain_nblast_v2`. Matches to the hemibrain
      (Scheffer et al., 2020) v1.2.1 dataset, a dense reconstruction of
      half the *Drosophila* brain.

  `"manc"`

  :   CAVE table: `banc_manc_nblast_v2`. Matches to the male adult nerve
      cord (Takemura et al., 2024) MANC v1.2.1 dataset.

  `"fanc"`

  :   CAVE table: `banc_fanc_nblast_v2`. Matches to FANC (Azevedo et
      al., 2024) v1116, a female adult nerve cord dataset.

- ...:

  Additional arguments passed to
  [`banc_cave_query`](https://flyconnectome.github.io/bancr/reference/banc_cave_query.md),
  including `live` (default `TRUE`; set to `2` for real-time results or
  `FALSE` for the latest materialised version).

## Value

A `data.frame` following the CAVE `cell_match` schema with columns:

- `id`:

  CAVE annotation ID (integer).

- `pt_root_id`:

  Current BANC root ID at the time of query (automatically updated by
  CAVE when neurons are edited).

- `pt_supervoxel_id`:

  Supervoxel ID anchoring the annotation to the segmentation. Stable
  across root ID changes.

- `pt_position`:

  3D position in voxel coordinates (resolution 4 x 4 x 45 nm)
  identifying the BANC neuron.

- `query_root_id`:

  BANC root ID at the time the NBLAST was run. May differ from
  `pt_root_id` if the neuron has since been edited.

- `match_id`:

  Identifier of the matched neuron in the target dataset. Format varies:
  hemibrain/maleCNS use `bodyid` (integer as string), FAFB uses
  `root_783` (FlyWire root ID), MANC uses `bodyid`, FANC uses `cell_id`.
  Mirrored matches are prefixed with `"m"` (e.g. `"m12345"`).

- `score`:

  NBLAST similarity score (0-1). Higher is more similar. Typical
  thresholds: 0.3-0.4 for strong matches.

- `validation`:

  Logical. `TRUE` if the match has been manually validated by a human
  annotator, `FALSE` otherwise.

## Details

These tables are populated by the bancpipeline NBLAST workflow
(`banc-nblast-compile.R` and `banc-nblast-cave.R`). The compile step
runs NBLAST morphological comparisons and writes results to feather
files; the CAVE sync step uploads new results and removes stale entries
(where the BANC neuron's root ID has changed since the NBLAST was run).

Because CAVE tracks root ID changes via the `pt_supervoxel_id` anchor
point, the `pt_root_id` column always reflects the current segmentation.
Compare `pt_root_id` with `query_root_id` to identify entries that may
need re-running (the neuron was edited after the NBLAST).

## See also

[`banc_cave_query`](https://flyconnectome.github.io/bancr/reference/banc_cave_query.md)
for the underlying CAVE query function,
[`banc_cave_tables`](https://flyconnectome.github.io/bancr/reference/banc_cave_tables.md)
for listing all available CAVE tables.

## Examples

``` r
if (FALSE) { # \dontrun{
# Get all maleCNS NBLAST matches
matches <- banc_nblast_matches("malecns")

# Get validated matches only
validated <- banc_nblast_matches("fafb") %>%
  dplyr::filter(validation == TRUE)

# Find matches for specific BANC neurons
my_matches <- banc_nblast_matches("hemibrain") %>%
  dplyr::filter(pt_root_id %in% my_root_ids)

# Identify stale entries (neuron edited since NBLAST)
stale <- banc_nblast_matches("manc") %>%
  dplyr::filter(pt_root_id != query_root_id)
} # }
```
