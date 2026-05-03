#' @title Build a Spelunker / ng.banc.community scene with a custom LM image layer
#'
#' @description Construct a Neuroglancer JSON state that overlays a
#' light-microscopy precomputed image layer (e.g. the output of
#' \code{neuronbridger::nrrd_to_precomputed} applied to a registered
#' confocal stack) on top of a public BANC scene that already contains the
#' BANC EM image, segmentation proofreading, region outlines and JRC2018F
#' atlas layers, and either return a long fragment-encoded URL, or — if a
#' Spelunker / CAVE token is available — POST the state to the
#' \code{nglstate/api/v1/post} endpoint and return the shortened
#' \code{...nglstate/api/v1/<id>} URL.
#'
#' @param lm_url URL of the precomputed LM layer (\code{gs://...},
#' \code{https://...}, or \code{file:///...} for local viewing).
#' \code{precomputed://} is added automatically if missing.
#' @param layer_name display name of the LM layer in Neuroglancer.
#' Default \code{"LM data"}.
#' @param shader optional Neuroglancer shader string. If \code{NULL} (the
#' default), a single-channel grayscale shader with a 4× contrast boost is
#' used (Neuroglancer's UI lets viewers tune it further).
#' @param opacity layer opacity in \code{[0, 1]}; default \code{0.6}.
#' @param blend layer blend mode (\code{"default"}, \code{"additive"}).
#' Default \code{"additive"} so the LM signal lights up where it overlaps
#' the EM rather than occluding it.
#' @param ids optional vector of BANC root IDs to add to the
#' \code{"segmentation proofreading"} layer that ships in the base state.
#' @param base_url Spelunker base URL with a starting state to extend.
#' Defaults to the public BANC scene used by \code{\link{banc_scene}}.
#' Pass an alternative scene URL to start from your own layered state.
#' @param shorten logical; if \code{TRUE}, POST the state to the
#' \code{nglstate/api/v1/post} endpoint and return a shortened URL of the
#' form \code{<base>/#!middleauth+https://global.daf-apis.com/nglstate/api/v1/<id>}.
#' Requires a valid token from \code{\link{banc_set_token}}.
#' @param viewer one of \code{"banc.ng.community/view"} (default; the
#' public BANC viewer at \code{https://banc.ng.community/view/}),
#' \code{"banc.ng.community"} (the private CAVE-authenticated BANC
#' viewer at \code{https://banc.ng.community/}, used by the BANC
#' team), \code{"spelunker"} (the upstream
#' \code{https://spelunker.cave-explorer.org/} viewer), or any other
#' base URL string (must end with \code{/}). Only affects the prefix
#' of the returned URL — the state JSON itself is identical and any
#' Spelunker-compatible viewer can render it via
#' \code{#!middleauth+https://...}.
#' @param open logical; if \code{TRUE}, open the result in the system
#' browser.
#'
#' @return A character of length 1: a long fragment URL by default, or a
#' shortened \code{nglstate/api/v1/<id>} URL when \code{shorten = TRUE}.
#'
#' @details The state is assembled by:
#' \enumerate{
#'   \item fetching the JSON for \code{base_url} (so the result inherits
#'     all standard public layers — BANC EM, segmentation proofreading,
#'     region outlines, JRC2018F atlas, FAFB / hemibrain / MANC imported
#'     meshes, synapse cloud, nuclei);
#'   \item appending an \code{image} layer pointing at \code{lm_url} with
#'     the supplied \code{shader}, \code{opacity} and \code{blend}; and
#'   \item optionally adding \code{ids} as visible segments in the
#'     proofreading layer.
#' }
#' If the base-state fetch fails (e.g. no token / offline), a minimal stub
#' state with just BANC EM + segmentation + LM is used instead.
#'
#' The companion converter \code{neuronbridger::nrrd_to_precomputed()}
#' takes any 3-D \code{.nrrd} (or in-memory array) and writes the
#' Neuroglancer "precomputed" format that this layer URL expects. See
#' \code{vignette("lm_layer_neuroglancer", package = "neuronbridger")}
#' for the full IS2 → JRC2018F → BANC pipeline.
#'
#' @section Bucket policy:
#' Lee-lab maintains a curated mirror of registered LM volumes at
#' \code{gs://lee-lab_brain-and-nerve-cord-fly-connectome/light_level/},
#' but that bucket is \strong{not public-write}. To produce a sharable
#' URL of your own you'll need either (a) write access to that bucket
#' (ask the lee-lab folks), or (b) your own public-read GCS / S3 / static
#' HTTP host. Once the precomputed directory is reachable over HTTPS,
#' pass its URL as \code{lm_url}.
#'
#' @examples
#' \dontrun{
#' # Convert + serve a registered LM volume:
#' neuronbridger::nrrd_to_precomputed(
#'   "CapaR_in_BANC.nrrd",
#'   output     = "/tmp/CapaR_pc",
#'   resolution = c(380, 380, 380)
#' )
#' system("gsutil -m cp -r /tmp/CapaR_pc gs://your-bucket/lm/CapaR/")
#'
#' u <- banc_lm_scene(
#'   "gs://your-bucket/lm/CapaR",
#'   layer_name = "Kondo 2020 - CapaR",
#'   viewer     = "ng.banc.community",
#'   shorten    = TRUE,
#'   open       = TRUE
#' )
#' }
#' @seealso \code{\link{banc_scene}},
#'   \code{neuronbridger::nrrd_to_precomputed}
#' @export
banc_lm_scene <- function(lm_url,
                          layer_name = "LM data",
                          shader     = NULL,
                          opacity    = 0.6,
                          blend      = c("additive", "default"),
                          ids        = NULL,
                          base_url   = paste0(
                            "https://spelunker.cave-explorer.org/#!middleauth+",
                            "https://global.daf-apis.com/nglstate/api/v1/",
                            "6431332029693952"),
                          shorten    = FALSE,
                          viewer     = "banc.ng.community/view",
                          open       = FALSE) {
  blend  <- match.arg(blend)
  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("Install the 'jsonlite' package.")
  if (!grepl("^precomputed://", lm_url))
    lm_url <- paste0("precomputed://", lm_url)
  if (is.null(shader))
    shader <- "void main() { emitGrayscale(toNormalized(getDataValue()) * 4.0); }"

  state <- tryCatch(banc_lm_fetch_state(base_url),
                    error = function(e) {
                      message("Could not fetch base state (", conditionMessage(e),
                              "); using minimal stub.")
                      banc_lm_minimal_state()
                    })

  state$layers <- c(state$layers,
                    list(banc_lm_image_layer(lm_url, layer_name, shader,
                                             opacity, blend)))

  if (!is.null(ids) && length(ids)) {
    seg_idx <- which(vapply(state$layers,
                            function(L) identical(L$name, "segmentation proofreading"),
                            logical(1)))
    if (length(seg_idx) == 1L) {
      cur <- state$layers[[seg_idx]]$segments
      state$layers[[seg_idx]]$segments <- unique(c(cur, banc_ids(ids)))
    } else {
      state$layers <- c(state$layers,
                        list(banc_lm_seg_layer(banc_ids(ids))))
    }
  }

  out <- if (shorten) banc_lm_shorten(state, base_url, viewer)
         else        banc_lm_long_url(state, base_url, viewer)
  if (open) utils::browseURL(out)
  invisible(out)
}

# ---------- internals ---------------------------------------------------

banc_lm_image_layer <- function(url, name, shader, opacity, blend) {
  list(type    = "image",
       source  = url,
       name    = name,
       shader  = shader,
       opacity = opacity,
       blend   = blend)
}

banc_lm_seg_layer <- function(ids) {
  list(type     = "segmentation",
       source   = "precomputed://gs://lee-lab_brain-and-nerve-cord-fly-connectome/segmentation",
       name     = "BANC seg",
       segments = as.character(ids))
}

banc_lm_minimal_state <- function() {
  list(
    dimensions = list(x = list(4e-9, "m"),
                      y = list(4e-9, "m"),
                      z = list(45e-9, "m")),
    layers = list(
      list(type   = "image",
           source = "precomputed://gs://zetta_lee_fly_cns_001_kisuk/final/v2/image",
           name   = "BANC EM"),
      list(type   = "segmentation",
           source = paste0("precomputed://gs://lee-lab_brain-and-nerve-cord-fly-connectome/",
                           "segmentation"),
           name   = "segmentation proofreading")
    ),
    selectedLayer = list(layer = "BANC EM"),
    layout = "xy-3d"
  )
}

banc_lm_fetch_state <- function(base_url) {
  if (!requireNamespace("fafbseg", quietly = TRUE))
    stop("Install 'fafbseg'.")
  parts <- strsplit(sub("#!middleauth+", "?", base_url, fixed = TRUE),
                    "?", fixed = TRUE)[[1]]
  json <- fafbseg::flywire_fetch(parts[2], token = banc_token(),
                                 return = "text", cache = TRUE)
  jsonlite::fromJSON(json, simplifyVector = FALSE)
}

banc_lm_view_base <- function(base_url, viewer) {
  switch(viewer,
    `banc.ng.community/view` = "https://banc.ng.community/view/",
    `banc.ng.community`      = "https://banc.ng.community/",
    spelunker                = "https://spelunker.cave-explorer.org/",
    {
      # Either a custom base URL, or fall back to whatever was in base_url
      if (grepl("^https?://", viewer)) {
        if (!endsWith(viewer, "/")) viewer <- paste0(viewer, "/")
        return(viewer)
      }
      base <- sub("#!.*$", "", base_url)
      if (!nzchar(base)) "https://banc.ng.community/view/" else base
    })
}

banc_lm_long_url <- function(state, base_url, viewer) {
  base <- banc_lm_view_base(base_url, viewer)
  json <- jsonlite::toJSON(state, auto_unbox = TRUE, null = "null")
  paste0(base, "#!", utils::URLencode(as.character(json), reserved = TRUE))
}

banc_lm_shorten <- function(state, base_url, viewer) {
  if (!requireNamespace("httr", quietly = TRUE))
    stop("Install 'httr'.")
  base <- banc_lm_view_base(base_url, viewer)
  json <- jsonlite::toJSON(state, auto_unbox = TRUE, null = "null")
  res <- httr::POST(
    "https://global.daf-apis.com/nglstate/api/v1/post",
    body   = json, encode = "raw",
    httr::add_headers(Authorization = paste("Bearer", banc_token()),
                      `Content-Type` = "application/json"))
  httr::stop_for_status(res)
  url <- httr::content(res, as = "text", encoding = "UTF-8")
  url <- gsub('"', "", url, fixed = TRUE)
  paste0(base, "#!middleauth+", url)
}
