# bancr 0.3.0

* new `banc_lm_scene()` for overlaying precomputed light-microscopy image
  layers on top of a public BANC Neuroglancer scene, with optional
  shortened-URL POST to `nglstate/api/v1/post`.
* new `banc_influence()` family (`banc_influence_arrow`,
  `banc_influence_duckdb`, `banc_influence_path`) for connectome influence
  scoring against parquet snapshots in GCS.
* read from BANC `synapses_v3`; with `details = TRUE`, `banc_synapses()` now
  returns neurotransmitter prediction scores.
* `bancsee()` gains a `clean_segments` option and surfaces NBLAST scores
  alongside synapse data; assorted scene-building fixes.
* read NBLAST match CAVE tables via `banc_nblast_matches()`.
* fix a Linux PATH_MAX crash in the Neuroglancer URL encoder when the
  scene JSON is long.
* internal: rework CI to use system Python via `actions/setup-python` (the
  reticulate miniconda env no longer ships pip), and clean up
  `R CMD check` Rd cross-reference / `\usage` warnings.

**Full Changelog**: https://github.com/flyconnectome/bancr/compare/v0.2.1...v0.3.0

# bancr 0.2.1

* fixes for `bancr::register_banc_coconat()` so that coconatfly can use new 
  banc dataset.

**Full Changelog**: https://github.com/flyconnectome/bancr/compare/v0.2.0...v0.2.1

# bancr 0.2.0

* first tagged release after BANC preprint

**Full Changelog**: https://github.com/flyconnectome/bancr/commits/v0.2.0
