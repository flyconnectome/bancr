# Build a BANC Neuroglancer scene with a custom LM image layer

Construct a Neuroglancer scene that overlays a light-microscopy
precomputed image layer (e.g. the output of
`neuronbridger::nrrd_to_precomputed` applied to a registered confocal
stack) on top of the canonical public BANC scene from
[`banc_scene`](https://flyconnectome.github.io/bancr/reference/banc_scene.md)
(BANC EM, segmentation proofreading, region outlines, JRC2018F atlas,
FAFB / hemibrain / MANC imported meshes, synapse cloud, nuclei) and
either return a long fragment-encoded URL or, when `shorten = TRUE`,
POST the state to the BANC state server and return a shortened
`spelunker.cave-explorer.org/#!middleauth+...nglstate/api/v1/<id>` URL —
the same URL form
[`bancsee`](https://flyconnectome.github.io/bancr/reference/bancsee.md)
produces.

## Usage

``` r
banc_lm_scene(
  lm_url,
  layer_name = "LM data",
  shader = NULL,
  range = c(1, 30),
  opacity = 0.55,
  blend = c("additive", "default"),
  volume_rendering = c("on", "max", "off"),
  volume_rendering_depth_samples = 788L,
  volume_rendering_gain = NULL,
  ids = NULL,
  url = NULL,
  shorten = TRUE,
  open = FALSE
)
```

## Arguments

- lm_url:

  URL of the precomputed LM layer (`gs://...`, `https://...`, or
  `file:///...` for local viewing). `precomputed://` is added
  automatically if missing.

- layer_name:

  display name of the LM layer in Neuroglancer. Default `"LM data"`.

- shader:

  optional Neuroglancer shader string. If `NULL` (default) the layer
  uses Neuroglancer's built-in `emitGrayscale` shader, controlled by
  `shaderControls.normalized.range` (set via `range`).

- range:

  numeric length-2 vector `c(low, high)` setting
  `shaderControls.normalized.range` — the source-data intensity window
  mapped to 0..1 brightness. Default `c(1, 30)`, which suits LM volumes
  that have been Elastix-warped + clipped to BANC voxel space (B-spline
  ringing leaves most of the signal in the bottom tenth of the uint8
  range). Increase the upper bound for brighter source data; tighten it
  for sparse stains.

- opacity:

  layer opacity in `[0, 1]`; default `0.55` (matching the public BANC
  `JRC2018F atlas imported` layer).

- blend:

  layer blend mode (`"default"`, `"additive"`). Default `"additive"` so
  LM signal lights up where it overlaps the EM rather than occluding it.

- volume_rendering:

  one of `"on"` (default), `"max"` or `"off"`. Required for the layer to
  be visible in 3-D Neuroglancer views. Cross-section / orthogonal slice
  views ignore this setting.

- volume_rendering_depth_samples:

  integer; how many depth samples Neuroglancer uses when ray-tracing the
  volume in 3-D. Default `788` (the value the public BANC atlas uses).

- volume_rendering_gain:

  numeric; 3-D volume rendering gain. `NULL` (default) leaves it unset,
  which lets Neuroglancer's UI control it interactively.

- ids:

  optional vector of BANC root IDs to add to the
  `"segmentation proofreading"` layer that ships in the base scene.

- url:

  Spelunker base URL with a starting state to extend. Defaults to the
  public BANC scene used by
  [`banc_scene`](https://flyconnectome.github.io/bancr/reference/banc_scene.md).

- shorten:

  logical; if `TRUE` (default), POST the state via the internal
  `banc_shorturl()` and return a shortened URL. If `FALSE`, return the
  long fragment-encoded URL with the state JSON inlined.

- open:

  logical; if `TRUE`, open the result in the system browser via
  [`browseURL`](https://rdrr.io/r/utils/browseURL.html).

## Value

A character of length 1: a shortened
`spelunker.cave-explorer.org/#!middleauth+...nglstate/api/v1/<id>` URL
by default, or a long fragment URL when `shorten = FALSE`.

## Viewers (`ng.banc.community/view/` / `ng.banc.community/`)

The BANC team also serves a public viewer at
[ng.banc.community/view/](https://ng.banc.community/view/) and a private
one at [ng.banc.community/](https://ng.banc.community/). **These viewers
do not POST states dynamically** — they load named states from
`ngstate.banc.community/view/<state-name>`, which redirects to static
JSON files committed under
[the-BANC-fly-connectome/neuroglancer_states/view/](https://github.com/jasper-tms/the-BANC-fly-connectome/tree/main/neuroglancer_states/view).
To publish a state to `ng.banc.community/view/`, save the state JSON
returned by this function (set `shorten = FALSE` and decode the
fragment, or use
[`fafbseg::ngl_decode_scene`](https://rdrr.io/pkg/fafbseg/man/ngl_decode_scene.html))
and PR it into the BANC repo. For ad-hoc / dev sharing, the
`shorten = TRUE` URL is the practical pattern.

## Bucket policy

Lee-lab maintains a curated mirror of registered LM volumes at
`gs://lee-lab_brain-and-nerve-cord-fly-connectome/light_level/`, but
that bucket is **not public-write**. To produce a sharable URL of your
own you'll need either (a) write access from the lee-lab maintainers, or
(b) your own public-read GCS / S3 / static HTTP host. Once the
precomputed directory is reachable, pass its URL as `lm_url`.

## See also

[`banc_scene`](https://flyconnectome.github.io/bancr/reference/banc_scene.md),
[`bancsee`](https://flyconnectome.github.io/bancr/reference/bancsee.md),
`neuronbridger::nrrd_to_precomputed`

## Examples

``` r
if (FALSE) { # \dontrun{
# Convert + serve a registered LM volume (the upstream
# vignettes do this in detail):
neuronbridger::nrrd_to_precomputed(
  "CapaR_in_JRC2018U_HR.nrrd",
  output     = "/tmp/CapaR_pc",
  resolution = c(519, 519, 1000)   # JRC2018U_HR voxel size in nm
)
system("gsutil -m cp -r /tmp/CapaR_pc gs://your-bucket/lm/CapaR/")

u <- banc_lm_scene(
  "gs://your-bucket/lm/CapaR",
  layer_name = "Kondo 2020 - CapaR",
  open       = TRUE
)
} # }
```
