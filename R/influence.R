#' Query BANC all-to-all influence scores
#'
#' @description Retrieve pre-computed influence scores between BANC neurons from
#'   partitioned parquet files stored on Google Cloud Storage. Influence scores
#'   quantify how much a "seed" (upstream) neuron's activity affects a
#'   "target" (downstream) neuron's steady-state response, based on the
#'   connectome's synaptic weight matrix (see Bates et al. 2020).
#'
#' @details The parquet files contain columns \code{upstream_id},
#'   \code{downstream_id}, and \code{raw_influence}. Adjusted influence is
#'   computed on-the-fly as \code{log(raw_influence) + const}, floored at 0.
#'
#'   Data is read from a local cache directory. On first use (or when
#'   \code{force_download=TRUE}), parquet files are downloaded from GCS using
#'   \code{gsutil}. Subsequent calls read directly from the cache.
#'
#'   The \code{"arrow"} method uses \code{\link[arrow]{open_dataset}} with
#'   predicate pushdown to scan only relevant chunks. The \code{"duckdb"}
#'   method registers a DuckDB view over the parquet files for fast SQL-based
#'   filtering. For small queries (few upstream/downstream IDs), both methods
#'   are fast; for large scans, \code{"duckdb"} may be faster.
#'
#' @param upstream_ids Character vector of upstream (seed) neuron root IDs.
#'   If \code{NULL}, all upstream neurons are included (use with caution —
#'   dataset is very large).
#' @param downstream_ids Character vector of downstream (target) neuron root
#'   IDs. If \code{NULL}, all downstream neurons are included.
#' @param const Numeric constant for adjusted influence calculation.
#'   Adjusted influence = \code{max(0, log(raw_influence) + const)}.
#'   Default 24 corresponds to a minimum meaningful influence of
#'   \code{exp(-24)} (approx 3.78e-11).
#' @param min_score Minimum adjusted influence score to return. Pairs with
#'   adjusted influence below this threshold are filtered out. Default 0
#'   returns all pairs with non-zero adjusted influence.
#' @param method Character, either \code{"arrow"} or \code{"duckdb"}.
#'   Controls which backend is used to read the parquet files.
#' @param local_path Path to a local directory containing the parquet files.
#'   If \code{NULL} (default), uses a cache directory under
#'   \code{tools::R_user_dir("bancr", "cache")}.
#' @param force_download Logical. If \code{TRUE}, re-download parquet files
#'   from GCS even if a local cache exists. Default \code{FALSE}.
#'
#' @return A \code{data.frame} with columns:
#' \describe{
#'   \item{upstream_id}{Character. Root ID of the upstream (seed) neuron.}
#'   \item{downstream_id}{Character. Root ID of the downstream (target) neuron.}
#'   \item{raw_influence}{Numeric. Raw steady-state influence score.}
#'   \item{adjusted_influence}{Numeric. \code{max(0, log(raw_influence) + const)}.}
#' }
#'
#' @export
#' @seealso \code{\link[influencer]{adjust_influence}} for grouped adjusted
#'   influence calculations.
#'
#' @references Bates, A.S., Schlegel, P., Roberts, R.J.V. et al. Complete
#'   Connectomic Reconstruction of Olfactory Projection Neurons in the Fly
#'   Brain. \emph{Curr Biol} 30, 3183-3199.e6 (2020).
#'
#' @examples
#' \dontrun{
#' # Get influence of one neuron on all targets
#' inf <- banc_influence(upstream_ids = "720575941521131930")
#' head(inf)
#'
#' # Get influence between specific pairs
#' inf <- banc_influence(
#'   upstream_ids = c("720575941521131930", "720575941478275714"),
#'   downstream_ids = c("720575941555924992")
#' )
#'
#' # Use duckdb backend
#' inf <- banc_influence(
#'   upstream_ids = "720575941521131930",
#'   method = "duckdb"
#' )
#'
#' # Only return strong connections
#' inf <- banc_influence(
#'   upstream_ids = "720575941521131930",
#'   min_score = 5
#' )
#' }
banc_influence <- function(upstream_ids = NULL,
                           downstream_ids = NULL,
                           const = 24,
                           min_score = 0,
                           method = c("arrow", "duckdb"),
                           local_path = NULL,
                           force_download = FALSE) {

  method <- match.arg(method)

  if (is.null(upstream_ids) && is.null(downstream_ids))
    stop("At least one of upstream_ids or downstream_ids must be specified")

  if (!is.null(upstream_ids))
    upstream_ids <- as.character(upstream_ids)
  if (!is.null(downstream_ids))
    downstream_ids <- as.character(downstream_ids)

  # Resolve local parquet directory
  parquet_dir <- banc_influence_path(local_path = local_path,
                                     force_download = force_download)

  # Read and filter using selected method
  if (method == "arrow") {
    res <- banc_influence_arrow(parquet_dir, upstream_ids, downstream_ids)
  } else {
    res <- banc_influence_duckdb(parquet_dir, upstream_ids, downstream_ids)
  }

  # Compute adjusted influence
  res$adjusted_influence <- log(res$raw_influence) + const
  res$adjusted_influence[res$adjusted_influence < 0] <- 0

  # Filter by minimum score
  if (min_score > 0) {
    res <- res[res$adjusted_influence >= min_score, , drop = FALSE]
  }

  res
}


# --- Internal helpers ---

#' Resolve local path to influence parquet files, downloading from GCS if needed
#' @keywords internal
banc_influence_path <- function(local_path = NULL, force_download = FALSE) {
  gs_url <- "gs://brain-and-nerve-cord_exports/brain_and_nerve_cord/influence/all_to_all/"

  if (is.null(local_path)) {
    local_path <- file.path(tools::R_user_dir("bancr", "cache"),
                            "influence", "all_to_all")
  }

  # Check if we already have parquet files
  existing <- list.files(local_path, pattern = "\\.parquet$", full.names = TRUE)

  if (!length(existing) || force_download) {
    message("Downloading influence parquet files from GCS...")
    message("  Source: ", gs_url)
    message("  Destination: ", local_path)
    dir.create(local_path, recursive = TRUE, showWarnings = FALSE)
    cmd <- sprintf("gsutil -m rsync %s %s", gs_url, local_path)
    status <- system(cmd)
    if (status != 0)
      stop("gsutil download failed (exit code ", status,
           "). Ensure gsutil is installed and authenticated.")
    existing <- list.files(local_path, pattern = "\\.parquet$", full.names = TRUE)
    message(sprintf("  Downloaded %d parquet files", length(existing)))
  }

  if (!length(existing))
    stop("No parquet files found in: ", local_path)

  local_path
}


#' Read influence scores using Arrow dataset
#' @keywords internal
banc_influence_arrow <- function(parquet_dir, upstream_ids, downstream_ids) {
  check_package_available("arrow")
  ds <- arrow::open_dataset(parquet_dir, format = "parquet")

  # Build filter expression
  if (!is.null(upstream_ids) && !is.null(downstream_ids)) {
    res <- ds |>
      dplyr::filter(.data$upstream_id %in% upstream_ids,
                     .data$downstream_id %in% downstream_ids) |>
      dplyr::collect()
  } else if (!is.null(upstream_ids)) {
    res <- ds |>
      dplyr::filter(.data$upstream_id %in% upstream_ids) |>
      dplyr::collect()
  } else {
    res <- ds |>
      dplyr::filter(.data$downstream_id %in% downstream_ids) |>
      dplyr::collect()
  }

  as.data.frame(res)
}


#' Read influence scores using DuckDB
#' @keywords internal
banc_influence_duckdb <- function(parquet_dir, upstream_ids, downstream_ids) {
  check_package_available("duckdb")
  check_package_available("DBI")

  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  # Register parquet files as a view
  parquet_glob <- file.path(parquet_dir, "*.parquet")
  DBI::dbExecute(con, sprintf(
    "CREATE VIEW influence AS SELECT * FROM parquet_scan('%s')",
    parquet_glob
  ))

  # Build SQL query with parameterised filters
  clauses <- character()
  if (!is.null(upstream_ids)) {
    ids_sql <- paste0("'", upstream_ids, "'", collapse = ", ")
    clauses <- c(clauses, sprintf("upstream_id IN (%s)", ids_sql))
  }
  if (!is.null(downstream_ids)) {
    ids_sql <- paste0("'", downstream_ids, "'", collapse = ", ")
    clauses <- c(clauses, sprintf("downstream_id IN (%s)", ids_sql))
  }

  where <- paste(clauses, collapse = " AND ")
  sql <- sprintf("SELECT * FROM influence WHERE %s", where)

  DBI::dbGetQuery(con, sql)
}

