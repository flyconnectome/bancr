# Pre-minted Spelunker viewer URLs for the BANC light-microscopy receptor library

A four-column flat record of pre-minted Spelunker
(`spelunker.cave-explorer.org`) viewer URLs for every light-microscopy
receptor / neuropeptide / orphan-GPCR layer currently mirrored under
`gs://lee-lab_brain-and-nerve-cord-fly-connectome/light_level/`. Each
layer is a Kondo et al. 2020 endogenously-tagged receptor GFP stack
bridged via CMTK + saalfeldlab template-building + Elastix into the BANC
voxel grid (see `neuronbridger::lm_to_banc_layer`). Each `ngl_link` is a
state ID resolvable forever by the BANC state server, so this table is a
stable shareable index.

For ad-hoc / interactive use, prefer `bancr:::banc_lm_volumes()` (the
live registry browser) +
[`banc_lm_scene`](https://flyconnectome.github.io/bancr/reference/banc_lm_scene.md)
(mints a fresh URL with custom render settings); `banc_lm_links` is the
canned snapshot for sharing or scripting.

## Usage

``` r
banc_lm_links
```

## Format

A `tibble` with one row per LM layer and columns:

- `source`:

  Upstream dataset citation (e.g. `"Kondo et al. 2020"`).

- `gene`:

  Receptor / peptide / gene short name (e.g. `"GluRIIA"`, `"5-HT2A"`,
  `"CapaR"`).

- `sample`:

  Sample identifier within the gene (Kondo's `no1` / `no2` prepared
  specimens).

- `ngl_link`:

  Spelunker URL of the form
  `https://spelunker.cave-explorer.org/#!middleauth+...nglstate/api/v1/<id>`,
  opening a scene that overlays the LM layer on the canonical public
  BANC scene.

## Source

Built by `data-raw/make_banc_lm_links.R` from the master registry at
`gs://lee-lab_brain-and-nerve-cord-fly-connectome/light_level/registry.json`.
Re-run the script (`Rscript data-raw/make_banc_lm_links.R`) to pick up
new layers or re-mint URLs with different defaults.

## Examples

``` r
if (FALSE) { # \dontrun{
# Browse all glutamate-receptor layers
subset(banc_lm_links, grepl("GluR|Nmdar|GluCl|mGluR|KaiR", gene))

# Open one in the browser
utils::browseURL(banc_lm_links$ngl_link[
  banc_lm_links$gene == "GluRIIA" & banc_lm_links$sample == "no1"])
} # }
```
