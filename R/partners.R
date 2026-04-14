#' Summarise the connectivity of BANC neurons
#'
#' Returns synaptically connected partners for specified neurons. Understanding
#' synaptic partnerships is crucial for analyzing neural circuits in the Brain And
#' Nerve Cord (BANC) connectome, revealing how distributed control architecture
#' coordinates behaviour across brain and ventral nerve cord regions.
#'
#' @details note that the rootids you pass in must be up to date. See example.
#'
#' @param rootids Character vector specifying one or more BANC rootids. As a
#'   convenience this argument is passed to \code{\link{banc_ids}} allowing you
#'   to pass in data.frames, BANC URLs or simple ids.
#' @param partners Character vector, either "outputs" or "inputs" to specify the direction of synaptic connections to retrieve.
#' @param threshold Integer threshold for minimum number of synapses (default 0).
#' @param remove_autapses Logical, whether to remove self-connections (default TRUE).
#' @param cleft.threshold Numeric threshold for cleft filtering (default 0).
#' @param datastack_name An optional CAVE \code{datastack_name}. If unset a
#'   sensible default is chosen.
#' @param synapse_table Character, the name of the synapse CAVE table you wish to use. Defaults to the latest.
#' @param ... Additional arguments passed to \code{\link[fafbseg]{flywire_partner_summary}}
#'
#' @return a data.frame
#' @seealso \code{\link{flywire_partner_summary}}, \code{\link{banc_latestid}}
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic connectivity analysis
#' sample_id=banc_latestid("720575941478275714")
#' head(banc_partner_summary(sample_id))
#' head(banc_partner_summary(sample_id, partners='inputs'))
#'
#' # Research application: Analyze descending neuron control circuits
#' library(dplyr)
#'
#' # Get DNa02 descending neurons that control walking behavior
#' dna02_annotations <- banc_codex_annotations() %>%
#'   filter(cell_type == "DNa02")
#' dna02_id <- dna02_annotations$pt_root_id[1]
#'
#' # Find their downstream targets in the VNC
#' dna02_outputs <- banc_partner_summary(dna02_id, partners='outputs') %>%
#'   slice_max(weight, n = 10)
#'
#' # Visualize the circuit in neuroglancer
#' banc_partner_summary(sample_id, partners='inputs') %>%
#'   slice_max(weight, n = 20) %>%
#'   banc_scene(open=TRUE)
#' }
#' @rdname banc_partner_summary
banc_partner_summary <- function(rootids,
                                 partners = c("outputs", "inputs"),
                                 synapse_table = c("synapses_v2","synapses_v3","synapses_v1"),
                                 threshold = 0,
                                 remove_autapses = TRUE,
                                 cleft.threshold = 0,
                                 datastack_name=NULL,
                                 ...) {
  if(is.null(datastack_name))
    datastack_name = banc_datastack_name()
  synapse_table=match.arg(synapse_table)
  with_banc(
    fafbseg::flywire_partner_summary(
      rootids,
      threshold = threshold,
      partners=partners,
      datastack_name = datastack_name,
      remove_autapses = remove_autapses,
      cleft.threshold = cleft.threshold,
      synapse_table=synapse_table,
      method = "cave",
      ...
    )
  )
}

### TODO ####
banc_datastack_name <- function() {
  if (requireNamespace("memoise", quietly = TRUE)) {
    return(memoise::memoise(banc_datastack_name_impl)())
  } else {
    return(banc_datastack_name_impl())
  }
}

banc_datastack_name_impl <- function() {
  banc_name <- tryCatch({
    cac=fafbseg::flywire_cave_client(NULL)
    datastacks=cac$info$get_datastacks()
    seldatastack=grep("brain_and_nerve_cord*", datastacks, value = T)
    if(length(seldatastack)==0)
      stop("Could not identify a banc production datastack amongst: ",
           paste(datastacks, collapse=','),
           "
Have you been granted access to banc production?")
    if(length(seldatastack)>1)
      warning("Multiple banc datastacks available; ",
              paste(seldatastack, collapse = ","),"
",
              "choosing: ", seldatastack[1])
    seldatastack[1]
  }, error = function(e){
    warning("using default setting")
    "brain_and_nerve_cord"
  })
  banc_name
}

#' @description \code{banc_partners} returns details of each unitary synaptic
#' connection (including its xyz location).
#'
#' @param details Logical. If \code{TRUE} and \code{synapse_table="synapses_v3"},
#'   additionally fetch \code{mean_score} and \code{median_score} from the
#'   reference tables \code{synapses_v3_mean_score} and
#'   \code{synapses_v3_median_score} and merge them onto the returned
#'   data.frame by synapse id. Default \code{FALSE}. Note: the median-score
#'   join is substantially slower than the mean-score join (on the order of
#'   10x), so only request \code{details=TRUE} when you need the per-synapse
#'   scores. Silently ignored for v2/v1.
#'
#' @examples
#' \dontrun{
#' # plot input and output synapses of a neuron
#' nclear3d()
#' fpi=banc_partners(banc_latestid("720575941478275714"), partners='in')
#' points3d(banc_raw2nm(fpi$post_pt_position), col='cyan')
#' fpo=banc_partners(banc_latestid("720575941478275714"), partners='out')
#' points3d(banc_raw2nm(fpo$pre_pt_position), col='red')
#'
#' # Compare results between the v2 (default) and v3 synapse tables
#' id <- banc_latestid("720575941478275714")
#' fpi_v2 <- banc_partners(id, partners='input', synapse_table="synapses_v2")
#' fpi_v3 <- banc_partners(id, partners='input', synapse_table="synapses_v3")
#' nrow(fpi_v2); nrow(fpi_v3)
#' # partner overlap
#' length(intersect(fpi_v2$pre_pt_root_id, fpi_v3$pre_pt_root_id))
#'
#' # Pull v3 synapses with mean and median scores attached (slower)
#' fpi_v3d <- banc_partners(id, partners='input', synapse_table="synapses_v3",
#'                          details=TRUE)
#' head(fpi_v3d[, c("id", "mean_score", "median_score")])
#' }
#' @export
#' @rdname banc_partner_summary
banc_partners <- function(rootids,
                          partners=c("input", "output"),
                          synapse_table = c("synapses_v2","synapses_v3","synapses_v1"),
                          details = FALSE,
                          ...) {
  partners=match.arg(partners)
  synapse_table=match.arg(synapse_table)
  rootids=banc_ids(rootids)
  fcc=banc_cave_client()
  pyids=fafbseg:::rids2pyint(rootids)
  res=if(partners=='input') {
    reticulate::py_call(fcc$materialize$synapse_query, post_ids=pyids, synapse_table=synapse_table, ...)
  } else {
    reticulate::py_call(fcc$materialize$synapse_query, pre_ids=pyids, synapse_table=synapse_table, ...)
  }
  # Flatten array-valued position columns to strings before R conversion
  # to avoid REAL()/character type mismatch in Arrow/reticulate
  tryCatch({
    cols <- res$columns$tolist()
    pos_cols <- grep("_position$", cols, value = TRUE)
    for (col in pos_cols) {
      res[[col]] <- res[[col]]$astype("str")
    }
  }, error = function(e) NULL)
  df <- fafbseg:::pandas2df(res)

  if (isTRUE(details) && identical(synapse_table, "synapses_v3") && nrow(df) > 0) {
    df <- .banc_attach_v3_scores(df, fcc)
  }
  df
}


# Fetch per-synapse mean and median scores from the v3 reference tables and
# merge onto `df` (which must have an `id` column from synapses_v3). CAVE
# rejects joining both reference tables in one query
# ("KeyError: 'synapses_v3_median_score'"), so we issue one reference-table
# query per score and join on target_id == synapse id. The median-score
# query is substantially slower than the mean-score query.
.banc_attach_v3_scores <- function(df, fcc) {
  syn_ids <- unique(df$id)
  syn_ids_py <- reticulate::r_to_py(as.list(as.character(syn_ids)))
  for (tbl in c("synapses_v3_mean_score", "synapses_v3_median_score")) {
    colnm <- sub("^synapses_v3_", "", tbl)
    out <- tryCatch({
      ref <- reticulate::py_call(
        fcc$materialize$query_table,
        table = tbl,
        filter_in_dict = reticulate::dict(target_id = syn_ids_py))
      ref_df <- fafbseg:::pandas2df(ref)
      if (!nrow(ref_df) || !all(c("target_id", "value") %in% colnames(ref_df)))
        return(NULL)
      data.frame(id = ref_df$target_id,
                 score = ref_df$value,
                 stringsAsFactors = FALSE)
    }, error = function(e) {
      warning(sprintf("Failed to fetch %s: %s", tbl, conditionMessage(e)))
      NULL
    })
    if (is.null(out)) next
    names(out)[2] <- colnm
    df <- merge(df, out, by = "id", all.x = TRUE, sort = FALSE)
  }
  df
}









