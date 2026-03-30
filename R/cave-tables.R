#' Read BANC CAVE-tables, good sources of metadata
#'
#' CAVE tables query functions that track neurons across segmentation changes
#' so that annotations and neuron entities can be stably tracked together.
#' The Brain And Nerve Cord (BANC) dataset represents the first complete
#' connectome including both brain and ventral nerve cord of a limbed animal,
#' comprising approximately 160,000 neurons across the entire central nervous system.
#'
#' @param rootids Character vector specifying one or more BANC rootids. As a
#'   convenience this argument is passed to \code{\link{banc_ids}} allowing you
#'   to pass in data.frames, BANC URLs or simple ids.
#' @param nucleus_ids Character vector specifying one or more BANC nucleus ids.
#'   The nucleus (\url{https://en.wikipedia.org/wiki/Cell_nucleus}) contains
#'   the cell body and provides a stable reference point for neuron identification.
#' @param rawcoords Logical, whether or not to convert from raw coordinates into nanometers. Default is `FALSE`.
#' @param select A regex term for the name of the table you want
#' @param datastack_name  Defaults to "brain_and_nerve_cord". See https://global.daf-apis.com/info/ for other options.
#' @param table Character, possible alternative tables for the sort of data frame the function returns. One must be chosen.
#' @param edgelist_view Character, name of prepared CAVE view that computes the proofread-neuron edgelist.
#' @param simplify logical, if \code{TRUE} then the proportion of presynaptic connections for each transmitter type is returned, for each query neuron.
#' @param ... Additional arguments passed to
#'   \code{fafbseg::\link{flywire_cave_query}} or \code{bancr:::get_cave_table_data}.
#'
#' @return A \code{data.frame} describing a CAVE-table related to the BANC project.
#' In the case of \code{banc_cave_tables}, a vector is returned containing the names of
#' all query-able cave tables.
#'
#' @details
#' CAVE tables store rich metadata supporting analysis of distributed neural
#' control across the entire central nervous system. For more information about
#' CAVE infrastructure, see \url{https://www.caveconnecto.me/CAVEclient/}.
#'
#' @seealso \code{fafbseg::\link{flywire_cave_query}}
#'
#' @export
#' @examples
#' \dontrun{
#' all_banc_soma_positions <- banc_nuclei()
#' points3d(nat::xyzmatrix(all_banc_soma_positions$pt_position))
#' }
#' @importFrom magrittr "%>%"
banc_cave_tables <- function(datastack_name = NULL,
                             select = NULL){
  if(is.null(datastack_name))
    datastack_name=banc_datastack_name()
  fac <- fafbseg::flywire_cave_client(datastack_name = datastack_name)
  dsinfo <- fac$info$get_datastack_info()
  tt <- fac$annotation$get_tables()
  if(!is.null(select)){
    chosen_tables <- grep(select, tt)
    if (length(chosen_tables) == 0)
      stop(sprintf("I cannot find a '%s' table for datastack: ", select),
           datastack_name, "
Please ask for help on #annotation_infrastructure https://flywire-forum.slack.com/archives/C01M4LP2Y2D")
    if (length(chosen_tables) == 1)
      return(tt[chosen_tables])
    chosen <- tt[rev(chosen_tables)[1]]
    warning(sprintf("Multiple candidate '%s' tables. Choosing: ", select),
            chosen)
    return(chosen)
  }else{
    return(tt)
  }
}

# Safe position conversion: skips empty/NA values instead of erroring
safe_raw2nm_position <- function(x) {
  if (length(x) == 0) return(x)
  valid <- !is.na(x) & nzchar(x)
  if (!any(valid)) return(x)
  x[valid] <- xyzmatrix2str(banc_raw2nm(x[valid]))
  x
}

#' @rdname banc_cave_tables
#' @export
banc_cave_views <- function(datastack_name = NULL,
                            select = NULL){
  if(is.null(datastack_name))
    datastack_name=banc_datastack_name()
  fac <- fafbseg::flywire_cave_client(datastack_name = datastack_name)
  dsinfo <- fac$info$get_datastack_info()
  tt <- unique(names(fac$materialize$get_views()))
  if(!is.null(select)){
    chosen_tables <- grep(select, tt)
    if (length(chosen_tables) == 0)
      stop(sprintf("I cannot find a '%s' view for datastack: ", select),
           datastack_name, "
Please ask for help on #annotation_infrastructure https://flywire-forum.slack.com/archives/C01M4LP2Y2D")
    if (length(chosen_tables) == 1)
      return(tt[chosen_tables])
    chosen <- tt[rev(chosen_tables)[1]]
    warning(sprintf("Multiple candidate '%s' views. Choosing: ", select),
            chosen)
    return(chosen)
  }else{
    return(tt)
  }
}

### edgelist ###

#' @rdname banc_cave_tables
#' @details
#' \code{banc_edgelist} returns a data frame of neuron-neuron connections where
#' the pre (presynaptic) neuron is upstream of the post (postsynaptic) neuron.
#' This edgelist contains synaptic connectivity data crucial for understanding
#' distributed neural control and behaviour-centric neural modules across the
#' brain-VNC boundary.
#' @export
banc_edgelist <- function(edgelist_view = c("synapses_v2_backbone_proofread_and_peripheral_nerves_counts",
                                            "synapses_v2_backbone_proofread_counts",
                                            "synapses_250226_backbone_proofread_and_peripheral_nerves_counts",
                                            "synapses_250226_backbone_proofread_counts",
                                            "synapses_v1_backbone_proofread_counts"),
                          ...){
  edgelist_view <- match.arg(edgelist_view)
  el <- with_banc(cave_view_query(edgelist_view, fetch_all_rows= TRUE, ...))
  el <- el %>%
    dplyr::arrange(dplyr::desc(.data$n))
  if(nrow(el)==500000|nrow(el)==1000000){
    warning("edgelist is exactly ", nrow(el), " rows, which is suspicious")
  }
  el
}

### mitochondria ###

#' @rdname banc_cave_tables
#' @export
banc_mitochondria <- function(rootids = NULL,
                              table = "mitochondria_v1",
                              rawcoords = FALSE, ...){
  cavec <- fafbseg:::check_cave()
  client <- try(cavec$CAVEclient(datastack_name=banc_datastack_name()))
  if(is.null(rootids)){
    res <- with_banc(get_cave_table_data(table, fetch_all_rows = TRUE, ...))
    if(nrow(res)==500000|nrow(res)==1000000){
      warning("dataframe is exactly ", nrow(res), " rows, which is suspicious")
    }
  }else{
    res <- client$materialize$tables[[table]](pt_root_id=rootids)$query()
  }
  if (isTRUE(rawcoords))
    res
  else {
    res %>% dplyr::mutate(dplyr::across(dplyr::ends_with("position"),
                          safe_raw2nm_position))
  }
}

### nuclei ###

#' @rdname banc_cave_tables
#' @export
banc_nuclei <- function(rootids = NULL,
                         nucleus_ids = NULL,
                         table = c("both","somas_v1a","somas_v1b"),
                         rawcoords = FALSE,
                         ...) {
  table <- match.arg(table)
  if(table=="both"){
    ba <- banc_nuclei(table="somas_v1a", nucleus_ids=nucleus_ids,rawcoords=rawcoords,...)
    bb <- banc_nuclei(table="somas_v1b", nucleus_ids=nucleus_ids,rawcoords=rawcoords,...)
    return(plyr::rbind.fill(bb,ba))
  }
  if (!is.null(rootids) & !is.null(nucleus_ids))
    stop("You must supply only one of rootids or nucleus_ids!")
  res <- if (is.null(rootids) && is.null(nucleus_ids))
    banc_cave_query(table = table, ...)
  else if (!is.null(rootids)) {
    rootids <- banc_ids(rootids)
    nuclei <- if (length(rootids) < 200)
      banc_cave_query(table =  table,
                      filter_in_dict = list(pt_root_id=rootids),
                      ...)
    else
      banc_cave_query(table =  table,
                      live = TRUE,
                      ...)
    if (nrow(nuclei) == 0)
      return(nuclei)
    nuclei <- nuclei %>%
      dplyr::right_join(data.frame(pt_root_id = as.integer64(rootids)),
                        by = "pt_root_id") %>%
      dplyr::select(colnames(nuclei))
    if (length(rootids) < 200) {
      nuclei
    }
    else {
      nuclei %>%
        dplyr::mutate(
          pt_root_id = with_banc(flywire_updateids(
            .data$pt_root_id,
            svids = .data$pt_supervoxel_id)))
    }
  } else {
    nuclei <- banc_cave_query(table = table,
                              filter_in_dict = list(id=nucleus_ids),
                              ...)
    nuclei %>%
      dplyr::right_join(data.frame(id = as.integer64(nucleus_ids)), by = "id") %>%
      dplyr::select(colnames(nuclei))
  }
  res$pt_position <- sapply(res$pt_position, paste, collapse=", ")
  # res$pt_position_ref <- sapply(res$pt_position_ref, paste, collapse=", ")
  res <- res %>%
    dplyr::rename(nucleus_id = .data$id,
                  nucleus_position = .data$pt_position,
                  root_id = .data$pt_root_id) %>%
    dplyr::filter(.data$valid=="t")
  if (isFALSE(rawcoords) && nrow(res) > 0) {
    valid_pos <- !is.na(res$nucleus_position) & nzchar(res$nucleus_position)
    res$nucleus_position_nm <- NA_character_
    if (any(valid_pos)) {
      res$nucleus_position_nm[valid_pos] <- gsub("\\(|\\)", "",
        apply(banc_raw2nm(res$nucleus_position[valid_pos]), 1, paste_coords))
    }
  }
  res
}

### neuron meta data ###

#' @rdname banc_cave_tables
#' @details
#' \code{banc_cell_info} accesses the cell_info CAVE table containing non-centralised
#' annotations from the research community for connectome neurones. These annotations
#' represent diverse contributions from researchers studying specific neural circuits
#' and cell types in the BANC dataset.
#' @export
#' @importFrom dplyr mutate ends_with across
#' @importFrom nat xyzmatrix2str
banc_cell_info <- function(rootids = NULL, rawcoords = FALSE, ...){
  table <- "cell_info"
  res <- with_banc(get_cave_table_data(table, ...))
  if (isTRUE(rawcoords))
    res
  else {
    res %>% mutate(across(ends_with("position"),
                          safe_raw2nm_position))
  }
}

#' @rdname banc_cave_tables
#' @export
banc_proofreading_notes <- function(rootids = NULL, rawcoords = FALSE, ...){
  table <- "proofreading_notes"
  res <- with_banc(get_cave_table_data(table, ...))
  if (isTRUE(rawcoords))
    res
  else {
    res %>% mutate(across(ends_with("position"),
                          safe_raw2nm_position))
  }
}

#' @rdname banc_cave_tables
#' @export
banc_cell_ids <- function(rootids = NULL,  ...){
  with_banc(get_cave_table_data('cell_ids', rootids, ...))
}

#' @rdname banc_cave_tables
#' @export
banc_neck_connective_neurons <- function(rootids = NULL,
                                         table = c("neck_connective_y92500", "neck_connective_y121000"),
                                         ...){
  table <- match.arg(table)
  with_banc(get_cave_table_data(table, rootids, ...))
}

#' @rdname banc_cave_tables
#' @export
banc_peripheral_nerves <- function(rootids = NULL, ...){
  with_banc(get_cave_table_data("peripheral_nerves", rootids, ...))
}

#' @rdname banc_cave_tables
#' @export
banc_backbone_proofread <- function(rootids = NULL, ...){
  with_banc(get_cave_table_data("backbone_proofread", rootids, ...))
}

#' Read NBLAST match results from CAVE
#'
#' Query cross-species NBLAST match results stored in CAVE tables. Each table
#' contains pairwise morphological similarity scores between BANC neurons and
#' neurons from another connectome dataset, computed using NBLAST
#' (Costa et al., 2016). Matches are identified by running all BANC neurons
#' against a target dataset and retaining hits above a score threshold.
#'
#' @param dataset Character, which cross-species NBLAST comparison to query.
#'   One of:
#'   \describe{
#'     \item{\code{"malecns"}}{CAVE table: \code{banc_malecns_nblast_v2}. Matches
#'       to the male CNS (Takemura et al., 2024) v0.9 dataset, covering the
#'       complete male central nervous system (~75K neurons).}
#'     \item{\code{"fafb"}}{CAVE table: \code{banc_fafb_nblast_v2}. Matches to
#'       FAFB (Zheng et al., 2018; Dorkenwald et al., 2024) FlyWire v783 dataset,
#'       a complete female \emph{Drosophila} brain connectome.}
#'     \item{\code{"hemibrain"}}{CAVE table: \code{banc_hemibrain_nblast_v2}. Matches
#'       to the hemibrain (Scheffer et al., 2020) v1.2.1 dataset, a dense
#'       reconstruction of half the \emph{Drosophila} brain.}
#'     \item{\code{"manc"}}{CAVE table: \code{banc_manc_nblast_v2}. Matches to the
#'       male adult nerve cord (Takemura et al., 2024) MANC v1.2.1 dataset.}
#'     \item{\code{"fanc"}}{CAVE table: \code{banc_fanc_nblast_v2}. Matches to
#'       FANC (Azevedo et al., 2024) v1116, a female adult nerve cord dataset.}
#'   }
#' @param ... Additional arguments passed to \code{\link{banc_cave_query}},
#'   including \code{live} (default \code{TRUE}; set to \code{2} for real-time
#'   results or \code{FALSE} for the latest materialised version).
#'
#' @return A \code{data.frame} following the CAVE \code{cell_match} schema with
#'   columns:
#'   \describe{
#'     \item{\code{id}}{CAVE annotation ID (integer).}
#'     \item{\code{pt_root_id}}{Current BANC root ID at the time of query
#'       (automatically updated by CAVE when neurons are edited).}
#'     \item{\code{pt_supervoxel_id}}{Supervoxel ID anchoring the annotation to
#'       the segmentation. Stable across root ID changes.}
#'     \item{\code{pt_position}}{3D position in voxel coordinates
#'       (resolution 4 x 4 x 45 nm) identifying the BANC neuron.}
#'     \item{\code{query_root_id}}{BANC root ID at the time the NBLAST was run.
#'       May differ from \code{pt_root_id} if the neuron has since been edited.}
#'     \item{\code{match_id}}{Identifier of the matched neuron in the target
#'       dataset. Format varies: hemibrain/maleCNS use \code{bodyid} (integer as
#'       string), FAFB uses \code{root_783} (FlyWire root ID), MANC uses
#'       \code{bodyid}, FANC uses \code{cell_id}. Mirrored matches are prefixed
#'       with \code{"m"} (e.g. \code{"m12345"}).}
#'     \item{\code{score}}{NBLAST similarity score (0-1). Higher is more
#'       similar. Typical thresholds: 0.3-0.4 for strong matches.}
#'     \item{\code{validation}}{Logical. \code{TRUE} if the match has been
#'       manually validated by a human annotator, \code{FALSE} otherwise.}
#'   }
#'
#' @details
#' These tables are populated by the bancpipeline NBLAST workflow
#' (\code{banc-nblast-compile.R} and \code{banc-nblast-cave.R}). The compile
#' step runs NBLAST morphological comparisons and writes results to feather
#' files; the CAVE sync step uploads new results and removes stale entries
#' (where the BANC neuron's root ID has changed since the NBLAST was run).
#'
#' Because CAVE tracks root ID changes via the \code{pt_supervoxel_id} anchor
#' point, the \code{pt_root_id} column always reflects the current segmentation.
#' Compare \code{pt_root_id} with \code{query_root_id} to identify entries that
#' may need re-running (the neuron was edited after the NBLAST).
#'
#' @export
#' @seealso \code{\link{banc_cave_query}} for the underlying CAVE query
#'   function, \code{\link{banc_cave_tables}} for listing all available CAVE
#'   tables.
#' @examples
#' \dontrun{
#' # Get all maleCNS NBLAST matches
#' matches <- banc_nblast_matches("malecns")
#'
#' # Get validated matches only
#' validated <- banc_nblast_matches("fafb") %>%
#'   dplyr::filter(validation == TRUE)
#'
#' # Find matches for specific BANC neurons
#' my_matches <- banc_nblast_matches("hemibrain") %>%
#'   dplyr::filter(pt_root_id %in% my_root_ids)
#'
#' # Identify stale entries (neuron edited since NBLAST)
#' stale <- banc_nblast_matches("manc") %>%
#'   dplyr::filter(pt_root_id != query_root_id)
#' }
banc_nblast_matches <- function(dataset = c("malecns", "fafb", "hemibrain", "manc", "fanc"),
                                ...) {
  dataset <- match.arg(dataset)
  table_name <- paste0("banc_", dataset, "_nblast_v2")
  banc_cave_query(table_name, ...)
}

# hidden
# Search for BANC neurons matching given target dataset IDs by NBLAST score.
# Uses CAVE NBLAST tables (banc_{dataset}_nblast_v2).
# Results are cached in a session-level environment to avoid repeated slow CAVE queries.
#
# @param match_ids Character vector of IDs in the target dataset to search for
# @param dataset Target dataset: "fafb", "manc", "hemibrain", "malecns", "fanc"
# @param top_n Number of top BANC hits to return (default 10)
# @param min_score Minimum NBLAST score threshold (default 0)
# @param cache If TRUE (default), cache the CAVE table for the session
# @param ... Additional arguments passed to banc_nblast_matches()
# @return Data frame of top BANC matches, one row per BANC neuron, with columns:
#   banc_id, max_score, best_match_id, n_matches
banc_nblast_search <- function(match_ids,
                               dataset = c("fafb", "manc", "hemibrain", "malecns", "fanc"),
                               top_n = 10L,
                               min_score = 0,
                               cache = TRUE,
                               ...) {
  dataset <- match.arg(dataset)
  match_ids <- as.character(match_ids)

  # Session-level cache for CAVE NBLAST tables (slow to query)
  cache_env <- ".banc_nblast_cache"
  if (!exists(cache_env, envir = .GlobalEnv)) {
    assign(cache_env, new.env(parent = emptyenv()), envir = .GlobalEnv)
  }
  cache_store <- get(cache_env, envir = .GlobalEnv)

  # Query CAVE for just the rows matching our target IDs (server-side filter)
  cache_key <- paste0("nblast_search_", dataset, "_", length(match_ids))
  if (cache && exists(cache_key, envir = cache_store)) {
    message(sprintf("Using cached %s NBLAST search results", dataset))
    hits <- get(cache_key, envir = cache_store)
  } else {
    table_name <- paste0("banc_", dataset, "_nblast_v2")
    message(sprintf("Querying CAVE %s for %d target IDs...", table_name, length(match_ids)))
    # Use Python CAVE client directly (banc_cave_query has serialization issues
    # with filter_in_dict on string columns)
    fcc <- banc_cave_client()
    ids_py <- reticulate::r_to_py(as.list(match_ids))
    py_env <- reticulate::py
    py_env$nblast_table <- table_name
    py_env$nblast_ids <- ids_py
    py_env$fcc_obj <- fcc
    reticulate::py_run_string("
_df = fcc_obj.materialize.query_table(nblast_table, filter_in_dict={'match_id': list(nblast_ids)})
for c in ['pt_root_id', 'query_root_id', 'pt_supervoxel_id']:
    if c in _df.columns:
        _df[c] = _df[c].astype(str)
_df['match_id'] = _df['match_id'].astype(str)
")
    nblast <- as.data.frame(reticulate::py$`_df`)
    hits <- nblast[nblast$score > min_score, ]
    if (cache) {
      assign(cache_key, hits, envir = cache_store)
      message(sprintf("  Cached %d hits for session", nrow(hits)))
    }
  }

  if (nrow(hits) == 0) {
    message(sprintf("No NBLAST matches found for the given IDs in %s.", dataset))
    return(data.frame(banc_id = character(), max_score = numeric(),
                      best_match_id = character(), n_matches = integer()))
  }

  message(sprintf("  Found %d NBLAST hits across %d BANC neurons",
                  nrow(hits), length(unique(hits$pt_root_id))))

  # For each BANC neuron, take max score across the requested match IDs
  top_hits <- hits %>%
    dplyr::mutate(banc_id = as.character(.data$pt_root_id)) %>%
    dplyr::group_by(.data$banc_id) %>%
    dplyr::summarise(
      max_score = max(.data$score, na.rm = TRUE),
      best_match_id = .data$match_id[which.max(.data$score)],
      n_matches = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(.data$max_score)) %>%
    dplyr::slice_head(n = top_n)

  top_hits
}

# hidden
banc_cave_cell_types <- function(cave_id = NULL, invert = FALSE, ...){
  banc.cell.info <- banc_cell_info(rawcoords = TRUE, ...)
  if(!is.null(cave_id)){
    if(invert){
      banc.cell.info <- banc.cell.info %>%
        dplyr::filter(!(.data$user_id %in% cave_id))
    }else{
      banc.cell.info <- banc.cell.info %>%
        dplyr::filter(.data$user_id %in% cave_id)
    }
  }
  banc.cell.info$pt_position <- sapply(banc.cell.info$pt_position, paste, collapse=", ")
  banc.cell.info.mod <- banc.cell.info %>%
    dplyr::filter(.data$valid == 't') %>%
    dplyr::rowwise() %>%
    dplyr::mutate(pt_position = paste0(.data$pt_position,collapse=",")) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(.data$pt_root_id) %>%
    dplyr::arrange(.data$pt_position, .data$tag2, .data$tag) %>%
    dplyr::mutate(side =  dplyr::case_when(
      grepl("^soma side",.data$tag2) ~ gsub("soma on |soma on ","",.data$tag),
      TRUE ~ NA
    )) %>%
    dplyr::mutate(cell_type = dplyr::case_when(
      grepl("neuron identity", .data$tag2) ~ .data$tag,
      !grepl(",",.data$tag) ~ .data$tag,
      TRUE ~ NA
    )) %>%
    dplyr::mutate(user_id = dplyr::case_when(
      !is.na(.data$cell_type) ~ .data$user_id,
      TRUE ~ NA
    )) %>%
    dplyr::mutate(cell_type = gsub("\\
.*|\\*.*","",.data$cell_type)) %>%
    dplyr::mutate(cell_class = dplyr::case_when(
      grepl("ascending|descending|descending|ascending", .data$tag) ~ .data$tag,
      grepl("sensory neuron|motor neuron|^trachea|^glia|^endocrine", .data$tag) ~ .data$tag,
      grepl("sensory neuron|motor neuron|^trachea|^glia|^endocrine", .data$tag) ~ .data$tag2,
      grepl("motor neuron", .data$tag) ~ "motor",
      grepl("endocrine", .data$tag) ~ "endocrine",
      grepl("central neuron", .data$tag2) ~ .data$tag,
      grepl("^innervates|^intersegmental", .data$tag) ~ .data$tag,
      TRUE ~ NA
    )) %>%
    dplyr::mutate(super_class = dplyr::case_when(
      grepl("ascend", .data$cell_class) ~ "ascending",
      grepl("descend", .data$cell_class) ~ "descending",
      grepl("sensory neuron", .data$cell_class) ~ "sensory",
      grepl("motor neuron", .data$cell_class) ~ "motor",
      grepl("endocrine", .data$cell_class) ~ "endocrine",
      grepl("efferent", .data$cell_class) ~ "efferent",
      grepl("optic", .data$cell_class) ~ "optic",
      grepl("optic", .data$tag) ~ "optic",
      grepl("optic", .data$tag2) ~ "optic",
      grepl("central", .data$cell_class) ~ "intrinsic",
      grepl("glia", .data$cell_class) ~ "glia",
      TRUE ~ NA
    )) %>%
    dplyr::mutate(notes = paste(unique(na.omit(sort(.data$tag))), collapse = ", "),
                  cell_class = paste(unique(na.omit(sort(.data$cell_class))), collapse = ", "),
                  super_class = paste(unique(na.omit(sort(.data$super_class))), collapse = ", "),
                  cell_type = paste(unique(na.omit(sort(.data$cell_type))), collapse = ", "),
                  side = paste(unique(na.omit(sort(.data$side))), collapse = ", "),
                  user_id = paste(unique(na.omit(sort(.data$user_id))), collapse = ", ")) %>%
    dplyr::ungroup() %>%
    dplyr::rename(cell_id = .data$id, root_id = .data$pt_root_id, supervoxel_id = .data$pt_supervoxel_id, position = .data$pt_position) %>%
    dplyr::distinct(.data$root_id, .data$supervoxel_id, .data$side, .data$super_class, .data$cell_class, .data$cell_type, .keep_all = TRUE) %>%
    dplyr::select(.data$cell_id, .data$root_id, .data$supervoxel_id, .data$position, .data$side, .data$super_class, .data$cell_class, .data$cell_type, .data$user_id, .data$notes) %>%
    dplyr::left_join(banc_users %>% dplyr::distinct(.data$pi_lab, .data$cave_id) %>% dplyr::mutate(cave_id=as.character(.data$cave_id)),
                     by=c("user_id"="cave_id")) %>%
    dplyr::rename(cell_type_source = .data$pi_lab)
  banc.cell.info.mod
}

# # # Updated cell_type_source column based on CAVE
# banc.cell.info.mod <- banc_cave_cell_types()
# banc.cell.info.mod <- subset(banc.cell.info.mod, ! user_id %in% c(355,52))
# bc.all <- banctable_query("SELECT _id, root_id, cell_type, other_names, super_class, cell_class, proofread, region, cell_type_source from banc_meta")
# bc.all$cell_type_source <- unlist(sapply(bc.all$cell_type_source ,function(x) paste(unlist(x),collapse=", ")))
# bc.ct <- bc.all %>%
#   dplyr::left_join(banc.cell.info.mod %>%
#                      dplyr::mutate(root_id=as.character(root_id)) %>%
#                      dplyr::distinct(root_id, cell_type, cell_type_source),
#                    by = "root_id") %>%
#   dplyr::mutate(
#     other_names = ifelse(is.na(other_names),'',other_names),
#     cell_type_source.y = gsub("Rachel Wilson Lab", "Wilson lab", cell_type_source.y),
#     cell_type_source.y = ifelse(is.na(cell_type_source.y),NA,tolower(cell_type_source.y)),
#     cell_type_source.x = ifelse(is.na(cell_type_source.x),NA,tolower(cell_type_source.x)),
#     cell_type_source.x = ifelse(grepl("NA|na|princeton|community|CAVE|Princeton",cell_type_source.x),NA,cell_type_source.x),
#     cell_type_source.x = ifelse(cell_type_source.x%in%c("","NA"),NA,cell_type_source.x),
#     cell_type_source.y = ifelse(cell_type_source.y%in%c("","NA"),NA,cell_type_source.y)) %>%
#   dplyr::mutate(cell_type = dplyr::case_when(
#     is.na(cell_type.x) ~ cell_type.y,
#     is.na(cell_type.y) ~ cell_type.x,
#     TRUE ~ cell_type.x),
#   ) %>%
#   dplyr::rowwise() %>%
#   dplyr::mutate(other_names = dplyr::case_when(
#     (!is.na(cell_type.x)&!is.na(cell_type.y)) &  (cell_type.y!= cell_type.x) ~ paste(sort(unique(c(unlist(strsplit(other_names,split=", ")),cell_type.y))),collapse=", "),
#     TRUE ~ other_names
#   )) %>%
#   dplyr::mutate(
#     cell_type_source.y = cell_type_source.y,
#     cell_type_source.x = cell_type_source.x,
#     cell_type_source = dplyr::case_when(
#     is.na(cell_type_source.x) ~ cell_type_source.y,
#     is.na(cell_type_source.y) ~ cell_type_source.x,
#     cell_type_source.x=="NA" ~ cell_type_source.y,
#     cell_type_source.y=="NA" ~ cell_type_source.x,
#     cell_type_source.x=="cave"&!is.na(cell_type_source.y) ~ cell_type_source.y,
#     cell_type_source.x=="community"&!is.na(cell_type_source.y) ~ cell_type_source.y,
#     cell_type_source.x==""&!is.na(cell_type_source.y) ~ cell_type_source.y,
#     !is.na(cell_type_source.x)&!is.na(cell_type_source.y) ~ paste(sort(unique(c(cell_type_source.x,cell_type_source.y)),
#                                                                        decreasing=TRUE),
#                                                                   collapse=","),
#     TRUE ~ cell_type_source.x
#   )) %>%
#   dplyr::filter(!is.na(cell_type_source), cell_type_source!="") %>%
#   dplyr::distinct(`_id`, root_id, .keep_all = TRUE) %>%
#   dplyr::select(`_id`, root_id, cell_type, other_names, cell_type_source,
#                   super_class, cell_class, proofread, region) %>%
#   dplyr::mutate(other_names = gsub("^,|^ ,|^ ","",other_names),
#                 cell_type_source = ifelse(cell_type_source=='151184',NA,cell_type_source))
#
# #Add cell type source labels
# bc.update <- as.data.frame(bc.ct)
# bc.update[is.na(bc.update)] <- ''
# banctable_update_rows(base='banc_meta',
#                       table = "banc_meta",
#                       df = bc.update[,c("_id","cell_type", "other_names", "cell_type_source")],
#                       append_allowed = FALSE,
#                       chunksize = 1000)


### neurotransmitters ###

#' @rdname banc_cave_tables
#' @export
#' @importFrom dplyr mutate ends_with across
#' @importFrom nat xyzmatrix2str
banc_nt_prediction <- function(rootids = NULL,
                               table = "synapses_250226_nt_prediction_5",
                               simplify = TRUE,
                               rawcoords = TRUE, ...){
  cavec <- fafbseg:::check_cave()
  client <- try(cavec$CAVEclient(datastack_name=banc_datastack_name()))
  if(is.null(rootids)){
    res <- with_banc(get_cave_table_data(table,
                                         rootids = rootids,
                                         fetch_all_rows = TRUE, ...))
    if(nrow(res)==500000|nrow(res)==1000000){
      warning("dataframe is exactly ", nrow(res), " rows, which is suspicious")
    }
  }else{
    res <- data.frame()
    for(rootid in rootids){
      res <- client$materialize$tables[[table]](pre_pt_root_id=rootid)$query()
      res$pre_pt_root_id <- rootid
      res <- plyr::rbind.fill(res, res)
      #ids <- fafbseg:::rids2pyint(rootid)
      # pyres <- if (method == "cave")
      #   reticulate::py_call(vol$chunkedgraph$get_roots, supervoxel_ids = ids,
      #                       ...)
      # else reticulate::py_call(vol$get_roots, ids, ...)
      #
      # res[[length(res) + 1]] = pyids2bit64(pyres, as_character = !integer64)
    }
  }
  if (!isTRUE(rawcoords)) {
    res <- res %>% dplyr::mutate(across(ends_with("position"),
                                        safe_raw2nm_position))
  }
  res <- res %>%
    dplyr::mutate(pre_pt_supervoxel_id = as.character(.data$pre_pt_supervoxel_id),
                  pre_pt_root_id = as.character(.data$pre_pt_root_id),
                  post_pt_supervoxel_id = as.character(.data$post_pt_supervoxel_id),
                  post_pt_root_id = as.character(.data$post_pt_root_id))
  if(simplify){
    nt.cols <- c("gaba", "serotonin", "acetylcholine", "dopamine", "octopamine", "glutamate", "histamine", "tyramine")
    res <- res %>%
      dplyr::filter(.data$valid_ref == 't', .data$valid == 't') %>%
      dplyr::arrange(dplyr::desc(.data$value)) %>%
      dplyr::distinct(.data$pre_pt_root_id, .data$id_ref, .keep_all = TRUE) %>%
      dplyr::group_by(.data$pre_pt_root_id, .data$tag) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop_last") %>%
      dplyr::mutate(count = sum(.data$n), prop = .data$n / .data$count) %>%
      dplyr::ungroup() %>%
      dplyr::select(.data$pre_pt_root_id, .data$count, .data$tag, .data$prop) %>%
      dplyr::mutate(prop = round(.data$prop, 4))

    if (!requireNamespace("tidyr", quietly = TRUE)) {
      stop("Package 'tidyr' is required for this function. Please install it with: install.packages('tidyr')")
    }

    res <- res %>%
      tidyr::pivot_wider(
        names_from = .data$tag,
        values_from = .data$prop,
        values_fill = 0
      )
    missing_nt_cols <- setdiff(nt.cols, names(res))
    if(length(missing_nt_cols) > 0) {
      res[missing_nt_cols] <- 0
    }
    res <- res %>%
      dplyr::select(.data$pre_pt_root_id, dplyr::all_of(nt.cols), dplyr::everything()) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        top_nt_p = max(dplyr::c_across(dplyr::all_of(nt.cols)), na.rm = TRUE),
        top_nt = nt.cols[which.max(dplyr::c_across(dplyr::all_of(nt.cols)))]
      ) %>%
      dplyr::ungroup()
  }
  res
}

### Make/edit cave tables ###\

# Validtate positions
banc_validate_positions <- function(positions,
                                    units = c("raw","nm")){
  # Input validation
  units <- match.arg(units)
  if(is.null(positions)) {
    stop("The 'positions' parameter cannot be NULL. Please provide 3D coordinates.")
  }

  # Validate positions format
  positions <- nat::xyzmatrix(positions)
  if(is.data.frame(positions)) {
    # For dataframes, check for X,Y,Z columns
    req_cols <- c("X", "Y", "Z")
    if(!all(tolower(colnames(positions)) %in% tolower(req_cols))) {
      stop("When providing a dataframe, it must contain columns named 'X', 'Y', 'Z' ")
    }
  } else if(is.vector(positions) && is.numeric(positions)) {
    # For vectors, check length
    if(length(positions) != 3) {
      stop("When providing a numeric vector, it must have exactly 3 elements (X,Y,Z)")
    }
  } else if(is.matrix(positions) && is.numeric(positions)) {
    # For matrices, check dimensions
    if(ncol(positions) != 3) {
      stop("When providing a matrix, it must have exactly 3 columns (X,Y,Z)")
    }
    positions <- as.data.frame(positions)
  } else {
    stop("'positions' must be either a dataframe with X,Y,Z columns, a numeric vector of length 3,
         or a matrix with 3 columns")
  }
  if(is.null(nrow(positions))){
    positions <- unlist(c(positions))
  }else if(nrow(positions)==1){
    positions <- unlist(c(positions))
  }

  # convert
  if(units=="nm"){
    positions <- banc_nm2raw(positions)
  }
  positions
}

#' Annotate positions as backbone proofread
#'
#' @description Mark specific positions as backbone proofread in the CAVE annotation system.
#' @param positions 3D coordinates in BANC space
#' @param user_id Integer user ID for the annotation
#' @param units Character, coordinate units - either "raw" or "nm"
#' @param proofread Logical, whether to mark as proofread (default TRUE)
#' @param datastack_name Optional datastack name
#' @return Data frame of added annotations
#' @export
#' @examples
#' \dontrun{
#' # Add an annotation to a point in raw voxel space
#' banc_annotate_backbone_proofread(c(117105, 240526, 5122), user_id = 355, units = "raw")
#'
#' # Add an annotation to a point in nm
#' banc_annotate_backbone_proofread(c(468420, 962104 ,230490), user_id = 355, units = "nm")
#'
#' # deannotate a point, only from points added with given user_id. Use user_id = NULL to remove from full pool
#' banc_deannotate(c(468420, 962104, 230490), user_id = 355, units = "nm", table = "backbone_proofread")
#' }

banc_annotate_backbone_proofread <- function(positions,
                                              user_id,
                                              units = c("raw", "nm"),
                                              proofread = TRUE,
                                              datastack_name = NULL)
{
  positions <- banc_validate_positions(positions = positions,
                                               units = units)
  cavec = fafbseg:::check_cave()
  np = reticulate::import("numpy")
  pd = reticulate::import("pandas")
  client = banc_service_account(datastack_name)
  annotations <- banc_backbone_proofread(live = 2) %>% dplyr::filter(.data$proofread ==
                                                                       eval(proofread))
  if (!nrow(annotations)) {
    stop("no annotations collected")
  }
  curr.positions <- do.call(rbind, annotations$pt_position)
  curr.positions <- as.data.frame(curr.positions)
  colnames(curr.positions) <- c("X", "Y", "Z")
  curr.positions$id <- annotations$id
  stage <- client$annotation$stage_annotations("backbone_proofread")
  if (is.data.frame(positions)) {
    positions.orig <- positions
    positions <- dplyr::anti_join(positions, as.data.frame(curr.positions),
                                  by = c("X", "Y", "Z"))
    cat("given positions already in backbone_proofread:",
        nrow(positions.orig) - nrow(positions), "
")
    if (!nrow(positions)) {
      stop("all positions already marked:", nrow(positions.orig))
    }
    valid_ids = banc_xyz2id(positions, rawcoords = TRUE)
    valid_ids_not_0 = valid_ids[valid_ids != "0"]
    positions = positions[valid_ids != "0", ]
    if (sum(valid_ids == "0")) {
      warning("given positions with invalid root_id: ",
              sum(valid_ids == "0"))
    }
    if (!nrow(positions)) {
      stop("no valid positions given")
    }

    result_ind <- integer(0)

    # Create a progress barn
    if (!requireNamespace("progress", quietly = TRUE)) {
      stop("Package 'progress' is required for this function. Please install it with: install.packages('progress')")
    }
    pb <- progress::progress_bar$new(
      format = "[:bar] :percent | ETA: :eta | :current/:total positions | Elapsed: :elapsedfull",
      total = nrow(positions),
      clear = FALSE,
      width = 80
    )

    for (i in 1:nrow(positions)) {
      # Update progress bar
      pb$tick()

      this_pos <- unlist(positions[i, ])
      # this_pos <- np$array(positions[i, ])
      this_id <- as.numeric(valid_ids_not_0[i])
      stage$add(valid = TRUE, pt_position = np$array(this_pos),
                user_id = as.integer(user_id), valid_id = this_id,
                proofread = proofread)
      this_result <- client$annotation$upload_staged_annotations(stage)
      result_ind <- c(result_ind, this_result)
      stage$clear_annotations()
    }
  }
  else {
    matching_rows <- which(apply(curr.positions[, 1:3], 1,
                                 function(row) all(row == positions)))
    point_exists <- length(matching_rows) > 0
    if (point_exists) {
      stop("given position already marked")
    }
    valid_id = banc_xyz2id(positions, rawcoords = TRUE)
    if (valid_id == "0") {
      stop("given position does not return a valid root_id")
    }
    stage$add(valid = TRUE, pt_position = np$array(positions),
              user_id = as.integer(user_id), valid_id = as.numeric(valid_id),
              proofread = proofread)
    result_ind <- client$annotation$upload_staged_annotations(stage)
  }

  if (is.data.frame(positions)) {
    pause_seconds <- nrow(positions) * 0.1
  }
  else {
    pause_seconds <- 0.1
  }
  Sys.sleep(pause_seconds)
  annotations <- banc_backbone_proofread(live = 2)
  annotations.new <- annotations %>% dplyr::filter(.data$id %in%
                                                     result_ind)
  cat("annotated", nrow(annotations.new), "entities with backbone proofread:",
      proofread, "")
  return(annotations.new)
}

# hidden
banc_deannotate_cave_table <- function(annotation_ids = NULL,
                                       positions = NULL,
                                       table = "proofreading_notes",
                                       user_id = NULL,
                                       units = c("raw","nm"),
                                       datastack_name = NULL,
                                       use_admin_creds = FALSE){
  # Validate positions
  if(!is.null(positions)){
    positions <- banc_validate_positions(positions=positions, units=units)
    # Read table and check annotations are added
    annotations <- with_banc(get_cave_table_data(table, live = 2))
    if(!is.null(user_id)){
      annotations <- annotations %>%
        dplyr::filter(.data$user_id %in% !!user_id)
    }
    curr.positions <- do.call(rbind,annotations$pt_position)
    if(is.data.frame(positions)){
      curr.positions <- as.data.frame(curr.positions)
      colnames(curr.positions) <- c("X","Y","Z")
      curr.positions$id <- annotations$id
      matches <- dplyr::left_join(positions, as.data.frame(curr.positions), by = c("X","Y","Z"))
      annotation_ids <- c(na.omit(matches$id))
    }else{
      matching_rows <- which(apply(curr.positions, 1, function(row) all(row == positions)))
      point_exists <- length(matching_rows) > 0
      annotation_ids <- annotations$id[matching_rows]
    }
    cat("pt_positions in ",  table, " match to", length(annotation_ids), "given points\n")
  }else if(is.null(annotation_ids)){
    stop("either annotation_ids or positions must be given")
  }
  if(length(annotation_ids)){
    # get table
    if (use_admin_creds) {
      client = banc_service_account(datastack_name)
    }
    else {
      client = fafbseg::flywire_cave_client(datastack_name = datastack_name)
    }

    # Delete specified IDs
    result <- client$annotation$delete_annotation(table, annotation_ids)

    # Read table and check annotations are added
    annotations.new <- with_banc(get_cave_table_data(table, live = 2)) %>%
      dplyr::filter(.data$id %in% annotation_ids)
    if(nrow(annotations.new)){
      warning('not all given positions removed from : missing annotation_ids')
    }
    cat("deannotated", length(result), "entities, valid set to FALSE\n")
    return(result)
  }else{
    invisible()
  }
}

# hidden

# # Get positions from the BANC seatable
# banc.meta <- banctable_query("SELECT root_id, position, status, l2_nodes from banc_meta") %>%
# dplyr::filter(grepl("DEBRIS",status),
#               !is.na(position),
#               l2_nodes<10)
# positions <- nat::xyzmatrix(banc.meta$position)
#
# # Or define them somehow yourself
# positions <- structure(list(X = c(167501L, 166837L, 99722L),
#                             Y = c(143868L,149546L, 239936L),
#                             Z = c(5456L, 2485L, 6792L)),
#                        row.names = c(NA, 3L),
#                        class = "data.frame")
#
# # Annotate points marked as debris
# banc_annotate_proofreading_notes(positions = positions,
#                                  user_id = 355,
#                                  label = "debris",
#                                  units = "raw",
#                                  use_admin_creds = TRUE)
#
# # See what is in the table
# banc.proofreading.notes <- banc_proofreading_notes(rawcoords = TRUE)
#
# # Flag for deletion
# positions.str <- nat::xyzmatrix2str(positions)
# noted.positions <- nat::xyzmatrix2str(banc.proofreading.notes$pt_position)
# flagged.positions <- noted.positions %in% positions.str
# banc.proofreading.notes.delete <- banc.proofreading.notes[flagged.positions,] %>%
#   dplyr::filter(tag == "debris")
#
# # Deannotate the points
# banc_deannotate_cave_table(annotation_ids = banc.proofreading.notes.delete$id,
#                           positions = nat::xyzmatrix(banc.proofreading.notes.delete$pt_position),
#                           user_id = 355,
#                           units = "raw",
#                           table = "proofreading_notes",
#                           use_admin_creds = TRUE)

# Hidden
banc_annotate_proofreading_notes <- function(positions,
                                              user_id,
                                              label,
                                              units = c("raw", "nm"),
                                              datastack_name = NULL,
                                              use_admin_creds = FALSE){
  positions <- banc_validate_positions(positions = positions,
                                       units = units)
  cavec <- fafbseg:::check_cave()
  np <- reticulate::import("numpy")
  pd <- reticulate::import("pandas")
  annotations <- banc_proofreading_notes(live = 2, rawcoords = TRUE) %>%
    dplyr::filter(.data$tag == eval(label))
  if(nrow(annotations)){
    curr.positions <- do.call(rbind, annotations$pt_position)
    curr.positions <- as.data.frame(curr.positions)
    colnames(curr.positions) <- c("X", "Y", "Z")
    curr.positions$id <- annotations$id
  }
  if (is.data.frame(positions)) {
      positions.orig <- positions
      positions <- dplyr::anti_join(positions, as.data.frame(curr.positions),
                                    by = c("X", "Y", "Z"))
      cat("given positions already present with the same tag:",
          nrow(positions.orig) - nrow(positions), "")
      if (!nrow(positions)) {
        stop("all positions already marked:", nrow(positions.orig))
      }
  }
  banc_annotate_bound_tag_user_cave_table(positions,
                                          user_id=user_id,
                                          column = "tag",
                                          tag = label,
                                          table = "proofreading_notes",
                                          use_admin_creds = use_admin_creds)
}

# hidden
banc_annotate_bound_tag_user_cave_table <- function(positions,
                                     column,
                                     tag,
                                     table,
                                     user_id,
                                     datastack_name = NULL,
                                     use_admin_creds = FALSE){
  if(is.null(datastack_name)){
    datastack_name <- banc_datastack_name()
  }
  if (use_admin_creds) {
    client <- banc_service_account(datastack_name)
  }else {
    client <- fafbseg::flywire_cave_client(datastack_name = datastack_name)
  }
  np <- reticulate::import("numpy")
  stage <- client$annotation$stage_annotations(table)
  cat("checking ", nrow(positions),"positions....\n")
  valid_ids <- banc_xyz2id(positions, rawcoords = TRUE)
  valid_ids_not_0 <- valid_ids[valid_ids != "0"]
  positions <- positions[valid_ids != "0", ]
    if (sum(valid_ids == "0")) {
      warning("given positions with invalid root_id: ",
              sum(valid_ids == "0"))
    }
    if (!nrow(positions)) {
      stop("no valid positions given")
    }
    result_ind <- integer(0)
    # Create a progress bar
    if (!requireNamespace("progress", quietly = TRUE)) {
      stop("Package 'progress' is required for this function. Please install it with: install.packages('progress')")
    }
    pb <- progress::progress_bar$new(
      format = "[:bar] :percent | ETA: :eta | :current/:total positions | Elapsed: :elapsedfull",
      total = nrow(positions),
      clear = FALSE,
      width = 80
    )
    cat("updating ", table,  " ...\n")
    for (i in 1:nrow(positions)) {
      # Update progress bar
      pb$tick()
      this_pos <- unlist(positions[i, ])
      # this_pos <- np$array(positions[i, ])
      this_id <- as.numeric(valid_ids_not_0[i])
      do.call(stage$add, c(
        list(
          valid = TRUE,
          pt_position = np$array(this_pos),
          #valid_id = this_id
          user_id = as.integer(user_id)
        ),
        setNames(list(tag), column)
      ))
      this_result <- client$annotation$upload_staged_annotations(stage)
      result_ind <- c(result_ind, this_result)
      stage$clear_annotations()
    }
  if (is.data.frame(positions)) {
    pause_seconds <- nrow(positions) * 0.1
  }else {
    pause_seconds <- 1
  }
  cat("checking result ...")
  Sys.sleep(5+pause_seconds)
  annotations <- try({with_banc(get_cave_table_data(table, live = 2))})
  if(is.null(annotations)){
    annotations.new <- annotations %>% dplyr::filter(.data$id %in%result_ind)
    cat("annotated", nrow(annotations.new), "entities with:",tag, " for ", column)
    annotations.new
  }else{
    warning("could not check annotations were added, in a few mins try: bancr:::with_banc(get_cave_table_data(table, live = 2))")
    NULL
  }
}

#' Read BANC-FlyWireCodex annotation table
#'
#' Provides access to centralised cell type annotations from the BANC core team,
#' which are the official annotations available on FlyWireCodex. These standardised
#' annotations ensure consistency across the dataset and serve as the authoritative
#' cell type classifications for the BANC connectome.
#'
#' @param rootids Character vector specifying one or more BANC rootids. As a
#'   convenience this argument is passed to \code{\link{banc_ids}} allowing you
#'   to pass in data.frames, BANC URLs or simple ids.
#' @param live logical, get the most recent data or pull from the latest materialisation
#' @param ... method passed to \code{\link{banc_cave_query}}.
#'
#' @return A \code{data.frame} describing that should be similar to what you find for BANC
#' in FlyWireCodex.
#'
#' @details
#' This function accesses centralised cell type annotations curated by the BANC core
#' team, in contrast to \code{\link{banc_cell_info}} which contains non-centralised
#' annotations from the broader research community. The centralised annotations provide
#' standardised cell type classifications that are displayed on FlyWireCodex and serve
#' as the official reference for BANC cell types.
#'
#' @seealso \code{\link{banc_cave_tables}}, \code{\link{banc_cell_info}}
#'
#' @export
#' @examples
#' \dontrun{
#' banc.meta <- banc_codex_annotations()
#' }
banc_codex_annotations <- function (rootids = NULL, live = TRUE, ...){
  table_name <- "codex_annotations"

  if (!is.null(rootids)) {
    # If rootids are specified, query for those specific rootids
    rootids <- banc_ids(rootids)
    if (length(rootids) < 200) {
      codex_annotations <- banc_cave_query(table_name,
                                           live = live,
                                           filter_in_dict = list(pt_root_id = rootids),
                                           ...)
    } else {
      # For large numbers of rootids, get all data and filter
      codex_annotations_part_1 <- banc_cave_query(table_name, live = live,
                                                  limit = 500000, ...)
      codex_annotations_part_2 <- banc_cave_query(table_name, live = live,
                                                  offset = 500000, limit = 350000, ...)
      codex_annotations_part_3 <- banc_cave_query(table_name, live = live,
                                                  offset = 850000, ...)
      codex_annotations_part_1 <- codex_annotations_part_1 %>%
        dplyr::mutate(cell_type = as.character(.data$cell_type))
      codex_annotations_part_2 <- codex_annotations_part_2 %>%
        dplyr::mutate(cell_type = as.character(.data$cell_type))
      codex_annotations_part_3 <- codex_annotations_part_3 %>%
        dplyr::mutate(cell_type = as.character(.data$cell_type))
      codex_annotations <- dplyr::bind_rows(codex_annotations_part_1,
                                            codex_annotations_part_2,
                                            codex_annotations_part_3) %>%
        dplyr::filter(.data$pt_root_id %in% rootids)
    }
  } else {
    # Get all data if no rootids specified
    codex_annotations_part_1 <- banc_cave_query(table_name, live = live,
                                                limit = 500000, ...)
    codex_annotations_part_2 <- banc_cave_query(table_name, live = live,
                                                offset = 500000, limit = 350000, ...)
    codex_annotations_part_3 <- banc_cave_query(table_name, live = live,
                                                offset = 850000, ...)
    codex_annotations_part_1 <- codex_annotations_part_1 %>%
      dplyr::mutate(cell_type = as.character(.data$cell_type))
    codex_annotations_part_2 <- codex_annotations_part_2 %>%
      dplyr::mutate(cell_type = as.character(.data$cell_type))
    codex_annotations_part_3 <- codex_annotations_part_3 %>%
      dplyr::mutate(cell_type = as.character(.data$cell_type))
    codex_annotations <- dplyr::bind_rows(codex_annotations_part_1,
                                          codex_annotations_part_2,
                                          codex_annotations_part_3)
  }

  codex_annotations <- codex_annotations %>%
    dplyr::mutate(cell_type = as.character(.data$cell_type))

  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package 'tidyr' is required for this function. Please install it with: install.packages('tidyr')")
  }

  if (live == 2) {
    codex_annotations_flat_table <- codex_annotations %>%
      dplyr::group_by(.data$target_id, .data$classification_system) %>%
      dplyr::summarise(cell_type_combined = paste(unique(.data$cell_type),
                                                  collapse = ", "), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = .data$classification_system,
                         values_from = .data$cell_type_combined, values_fill = NA_character_)
  }
  else {
    codex_annotations_flat_table <- codex_annotations %>%
      dplyr::group_by(.data$target_id, .data$classification_system,
                      .data$pt_supervoxel_id, .data$pt_root_id, .data$pt_position) %>%
      dplyr::summarise(cell_type_combined = paste(unique(.data$cell_type),
                                                  collapse = ", "), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = .data$classification_system,
                         values_from = .data$cell_type_combined, values_fill = NA_character_)
  }

  return(codex_annotations_flat_table)
}

#' @rdname banc_cave_tables
#' @export
banc_version <- function() {
  bcc=banc_cave_client()
  ver=bcc$materialize$version
  ver
}
