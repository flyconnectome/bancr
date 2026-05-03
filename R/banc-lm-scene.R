#' @title Build a BANC Neuroglancer scene with a custom LM image layer
#'
#' @description Construct a Neuroglancer scene that overlays a
#' light-microscopy precomputed image layer (e.g. the output of
#' \code{neuronbridger::nrrd_to_precomputed} applied to a registered
#' confocal stack) on top of the canonical public BANC scene from
#' \code{\link{banc_scene}} (BANC EM, segmentation proofreading, region
#' outlines, JRC2018F atlas, FAFB / hemibrain / MANC imported meshes,
#' synapse cloud, nuclei) and either return a long fragment-encoded URL
#' or, when \code{shorten = TRUE}, POST the state to the BANC state
#' server and return a shortened
#' \code{spelunker.cave-explorer.org/#!middleauth+...nglstate/api/v1/<id>}
#' URL — the same URL form
#' \code{\link{bancsee}} produces.
#'
#' @param lm_url URL of the precomputed LM layer (\code{gs://...},
#' \code{https://...}, or \code{file:///...} for local viewing).
#' \code{precomputed://} is added automatically if missing.
#' @param layer_name display name of the LM layer in Neuroglancer.
#' Default \code{"LM data"}.
#' @param shader optional Neuroglancer shader string. If \code{NULL}
#' (default) the layer uses Neuroglancer's built-in \code{emitGrayscale}
#' shader, controlled by \code{shaderControls.normalized.range} (set
#' via \code{range}).
#' @param range numeric length-2 vector \code{c(low, high)} setting
#' \code{shaderControls.normalized.range} — the source-data intensity
#' window mapped to 0..1 brightness. Default \code{c(1, 30)}, which
#' suits LM volumes that have been Elastix-warped + clipped to BANC
#' voxel space (B-spline ringing leaves most of the signal in the bottom
#' tenth of the uint8 range). Increase the upper bound for brighter
#' source data; tighten it for sparse stains.
#' @param opacity layer opacity in \code{[0, 1]}; default \code{0.55}
#' (matching the public BANC \code{JRC2018F atlas imported} layer).
#' @param blend layer blend mode (\code{"default"}, \code{"additive"}).
#' Default \code{"additive"} so LM signal lights up where it overlaps
#' the EM rather than occluding it.
#' @param volume_rendering one of \code{"on"} (default), \code{"max"} or
#' \code{"off"}. Required for the layer to be visible in 3-D
#' Neuroglancer views. Cross-section / orthogonal slice views ignore
#' this setting.
#' @param volume_rendering_depth_samples integer; how many depth
#' samples Neuroglancer uses when ray-tracing the volume in 3-D.
#' Default \code{788} (the value the public BANC atlas uses).
#' @param volume_rendering_gain numeric; 3-D volume rendering gain.
#' \code{NULL} (default) leaves it unset, which lets Neuroglancer's UI
#' control it interactively.
#' @param ids optional vector of BANC root IDs to add to the
#' \code{"segmentation proofreading"} layer that ships in the base
#' scene.
#' @param url Spelunker base URL with a starting state to extend.
#' Defaults to the public BANC scene used by \code{\link{banc_scene}}.
#' @param shorten logical; if \code{TRUE} (default), POST the state via
#' the internal \code{banc_shorturl()} and return a shortened URL. If
#' \code{FALSE}, return the long fragment-encoded URL with the state
#' JSON inlined.
#' @param open logical; if \code{TRUE}, open the result in the system
#' browser via \code{\link[utils]{browseURL}}.
#'
#' @return A character of length 1: a shortened
#' \code{spelunker.cave-explorer.org/#!middleauth+...nglstate/api/v1/<id>}
#' URL by default, or a long fragment URL when \code{shorten = FALSE}.
#'
#' @section Viewers (\code{ng.banc.community/view/} / \code{ng.banc.community/}):
#' The BANC team also serves a public viewer at
#' \href{https://ng.banc.community/view/}{ng.banc.community/view/} and a
#' private one at
#' \href{https://ng.banc.community/}{ng.banc.community/}. \strong{These
#' viewers do not POST states dynamically} — they load named states from
#' \code{ngstate.banc.community/view/<state-name>}, which redirects to
#' static JSON files committed under
#' \href{https://github.com/jasper-tms/the-BANC-fly-connectome/tree/main/neuroglancer_states/view}{the-BANC-fly-connectome/neuroglancer_states/view/}.
#' To publish a state to \code{ng.banc.community/view/}, save the
#' state JSON returned by this function (set \code{shorten = FALSE} and
#' decode the fragment, or use \code{fafbseg::ngl_decode_scene}) and PR
#' it into the BANC repo. For ad-hoc / dev sharing, the
#' \code{shorten = TRUE} URL is the practical pattern.
#'
#' @section Bucket policy:
#' Lee-lab maintains a curated mirror of registered LM volumes at
#' \code{gs://lee-lab_brain-and-nerve-cord-fly-connectome/light_level/},
#' but that bucket is \strong{not public-write}. To produce a sharable
#' URL of your own you'll need either (a) write access from the
#' lee-lab maintainers, or (b) your own public-read GCS / S3 / static
#' HTTP host. Once the precomputed directory is reachable, pass its
#' URL as \code{lm_url}.
#'
#' @examples
#' \dontrun{
#' # Convert + serve a registered LM volume (the upstream
#' # vignettes do this in detail):
#' neuronbridger::nrrd_to_precomputed(
#'   "CapaR_in_JRC2018U_HR.nrrd",
#'   output     = "/tmp/CapaR_pc",
#'   resolution = c(519, 519, 1000)   # JRC2018U_HR voxel size in nm
#' )
#' system("gsutil -m cp -r /tmp/CapaR_pc gs://your-bucket/lm/CapaR/")
#'
#' u <- banc_lm_scene(
#'   "gs://your-bucket/lm/CapaR",
#'   layer_name = "Kondo 2020 - CapaR",
#'   open       = TRUE
#' )
#' }
#' @seealso \code{\link{banc_scene}}, \code{\link{bancsee}},
#'   \code{neuronbridger::nrrd_to_precomputed}
#' @export
banc_lm_scene <- function(lm_url,
                          layer_name                       = "LM data",
                          shader                           = NULL,
                          range                            = c(1, 30),
                          opacity                          = 0.55,
                          blend                            = c("additive", "default"),
                          volume_rendering                 = c("on", "max", "off"),
                          volume_rendering_depth_samples   = 788L,
                          volume_rendering_gain            = NULL,
                          ids                              = NULL,
                          url                              = NULL,
                          shorten                          = TRUE,
                          open                             = FALSE) {
  blend            <- match.arg(blend)
  volume_rendering <- match.arg(volume_rendering)
  if (length(range) != 2L || !is.numeric(range))
    stop("`range` must be a length-2 numeric vector c(low, high).")
  if (!grepl("^precomputed://", lm_url))
    lm_url <- paste0("precomputed://", lm_url)

  # Start from the canonical public BANC scene as an ngscene object,
  # then append our LM image layer using fafbseg's helpers — the same
  # pattern bancsee() uses.
  scene <- fafbseg::ngl_decode_scene(if (is.null(url)) banc_scene() else banc_scene(url = url))

  lm_layer <- list(type    = "image",
                   source  = lm_url,
                   name    = layer_name,
                   tab     = "rendering",
                   opacity = opacity,
                   blend   = blend,
                   shaderControls = list(
                     normalized = list(range = list(range[[1]], range[[2]]))
                   ),
                   volumeRendering              = volume_rendering,
                   volumeRenderingDepthSamples  = volume_rendering_depth_samples)
  if (!is.null(shader))   lm_layer$shader              <- shader
  if (!is.null(volume_rendering_gain))
                          lm_layer$volumeRenderingGain <- volume_rendering_gain
  scene[["layers"]] <- c(scene[["layers"]], list(lm_layer))

  # Optionally light up some BANC seg IDs as visible segments
  if (!is.null(ids) && length(ids)) {
    layer_names <- names(fafbseg::ngl_layers(scene))
    target <- "segmentation proofreading"
    if (target %in% layer_names) {
      idx <- which(layer_names == target)
      cur <- scene[["layers"]][[idx]][["segments"]]
      scene[["layers"]][[idx]][["segments"]] <-
        unique(c(cur, banc_ids(ids)))
    }
  }

  u <- as.character(scene)
  if (shorten) u <- banc_shorturl(u)
  if (open) {
    utils::browseURL(u)
    invisible(u)
  } else {
    u
  }
}
