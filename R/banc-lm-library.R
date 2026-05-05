#' @title Browse the curated BANC LM volume library
#'
#' @description Hidden helper that fetches and filters the lee-lab
#' light-microscopy volume registry. The registry is a single JSON
#' manifest at
#' \code{gs://lee-lab_brain-and-nerve-cord-fly-connectome/light_level/registry.json}
#' listing every LM volume that has been Elastix-warped into BANC voxel
#' coordinates and uploaded to the public-read bucket. Each entry
#' carries enough metadata (\code{dataset}, \code{gene},
#' \code{sample}, \code{channel}, \code{gs_url},
#' \code{source_filename}) to identify a volume and feed its URL
#' straight into \code{\link{banc_lm_scene}}.
#'
#' @param dataset Optional dataset filter (case-insensitive substring
#'   match against the \code{dataset} column). Current master folders:
#'   \code{kondo_et_al_2020} (Kondo et al. 2020 endogenous
#'   neurotransmitter-receptor tags), \code{deng_et_al_2019} (Deng et
#'   al. 2019 dense-core-vesicle / neuropeptide drivers).
#' @param gene Optional gene / receptor / peptide filter
#'   (case-insensitive substring match against \code{gene}).
#' @param channel Optional channel filter (e.g. \code{"no2"} for the
#'   GFP/neuron channel from Kondo). \code{NULL} returns all channels.
#' @param refresh Logical; if \code{TRUE} bypass the per-session cache
#'   and re-fetch the manifest from GCS.
#'
#' @return A \code{tibble} with one row per matching volume and columns:
#' \describe{
#'   \item{\code{dataset}}{Master folder (e.g. \code{kondo_et_al_2020}).}
#'   \item{\code{gene}}{Receptor / peptide / gene short name.}
#'   \item{\code{sample}}{Sample identifier within the gene.}
#'   \item{\code{channel}}{Source-stack channel
#'     (\code{no1} = NC82 reference, \code{no2} = GFP / neuron).}
#'   \item{\code{name}}{Display name on GCS (\code{{gene}_{sample}}).}
#'   \item{\code{gs_url}}{\code{gs://...} precomputed URL — feed
#'     directly to \code{banc_lm_scene(lm_url = ...)}.}
#'   \item{\code{alignment_space}}{Always \code{"BANC"} for entries in
#'     this registry (volumes are stored on the BANC voxel grid).}
#'   \item{\code{voxdims_nm}}{Length-3 numeric, output voxel size in nm
#'     (typically \code{c(400, 400, 400)}).}
#'   \item{\code{source_filename}}{Original filename on the upstream
#'     repository (e.g. G-Node DOI host).}
#'   \item{\code{source_path}}{Original storage path on the upstream
#'     repository, for provenance.}
#'   \item{\code{uploaded}}{ISO date the volume entered the registry.}
#' }
#'
#' @section Registry schema:
#' The registry JSON is a single object with two top-level fields:
#' \code{schema_version} (integer, currently \code{1}) and
#' \code{volumes} (array of entries matching the columns documented
#' above). New uploads append to \code{volumes}; readers tolerate
#' unknown fields. Manifest path:
#' \code{gs://lee-lab_brain-and-nerve-cord-fly-connectome/light_level/registry.json}
#' (public-read; the parallel HTTPS URL
#' \code{https://storage.googleapis.com/lee-lab_brain-and-nerve-cord-fly-connectome/light_level/registry.json}
#' is what this function actually fetches).
#'
#' @section Status:
#' Hidden / internal — not exported. Intended as an exploratory helper
#' while the registry is being populated. Promote to an exported API
#' (and rename) once the schema and library scope stabilise.
#'
#' @examples
#' \dontrun{
#' # All currently-registered volumes:
#' bancr:::banc_lm_volumes()
#'
#' # Just the Kondo glutamate-receptor neuron channels:
#' bancr:::banc_lm_volumes(dataset = "kondo", channel = "no2")
#'
#' # Pick one and open it in a BANC scene:
#' v <- bancr:::banc_lm_volumes(gene = "GluRIA")
#' banc_lm_scene(v$gs_url[1], layer_name = v$name[1], open = TRUE)
#' }
#' @keywords internal
#' @noRd
banc_lm_volumes <- function(dataset = NULL,
                            gene    = NULL,
                            channel = NULL,
                            refresh = FALSE) {
  manifest <- .banc_lm_registry(refresh = refresh)
  if (!nrow(manifest)) return(manifest)

  if (!is.null(dataset))
    manifest <- manifest[grepl(dataset, manifest$dataset, ignore.case = TRUE), , drop = FALSE]
  if (!is.null(gene))
    manifest <- manifest[grepl(gene, manifest$gene, ignore.case = TRUE), , drop = FALSE]
  if (!is.null(channel))
    manifest <- manifest[manifest$channel %in% channel, , drop = FALSE]

  manifest
}

# Per-session cache for the registry JSON.
.banc_lm_cache <- new.env(parent = emptyenv())

.banc_lm_registry_url <- function() {
  "https://storage.googleapis.com/lee-lab_brain-and-nerve-cord-fly-connectome/light_level/registry.json"
}

.banc_lm_registry <- function(refresh = FALSE) {
  if (!refresh && !is.null(.banc_lm_cache$df)) return(.banc_lm_cache$df)

  url <- .banc_lm_registry_url()
  res <- try(httr::GET(url), silent = TRUE)
  if (inherits(res, "try-error") || httr::status_code(res) == 404L) {
    df <- .empty_lm_registry()
    .banc_lm_cache$df <- df
    return(df)
  }
  if (httr::status_code(res) >= 400L)
    stop("Could not fetch BANC LM registry (", httr::status_code(res),
         "): ", url)

  json <- jsonlite::fromJSON(httr::content(res, as = "text", encoding = "UTF-8"),
                             simplifyVector = TRUE)
  vols <- json[["volumes"]]
  if (is.null(vols) || !length(vols)) {
    df <- .empty_lm_registry()
  } else {
    df <- tibble::as_tibble(vols)
    expected <- c("dataset","gene","sample","channel","name","gs_url",
                  "alignment_space","voxdims_nm","source_filename",
                  "source_path","uploaded")
    for (col in setdiff(expected, names(df))) df[[col]] <- NA
    df <- df[, expected, drop = FALSE]
  }
  .banc_lm_cache$df <- df
  df
}

.empty_lm_registry <- function() {
  tibble::tibble(
    dataset         = character(),
    gene            = character(),
    sample          = character(),
    channel         = character(),
    name            = character(),
    gs_url          = character(),
    alignment_space = character(),
    voxdims_nm      = list(),
    source_filename = character(),
    source_path     = character(),
    uploaded        = character()
  )
}
