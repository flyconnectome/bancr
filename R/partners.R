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
                                 synapse_table = NULL,
                                 threshold = 0,
                                 remove_autapses = TRUE,
                                 cleft.threshold = 0,
                                 datastack_name=NULL,
                                 ...) {
  if(is.null(datastack_name))
    datastack_name = banc_datastack_name()
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
#'
#' @examples
#' \dontrun{
#' # plot input and output synapses of a neuron
#' nclear3d()
#' fpi=banc_partners(banc_latestid("720575941478275714"), partners='in')
#' points3d(banc_raw2nm(fpi$post_pt_position), col='cyan')
#' fpo=banc_partners(banc_latestid("720575941478275714"), partners='out')
#' points3d(banc_raw2nm(fpo$pre_pt_position), col='red')
#' }
#' @export
#' @rdname banc_partner_summary
banc_partners <- function(rootids,
                          partners=c("input", "output"),
                          synapse_table = NULL,
                          ...) {
  partners=match.arg(partners)
  rootids=banc_ids(rootids)
  fcc=banc_cave_client()
  pyids=fafbseg:::rids2pyint(rootids)
  res=if(partners=='input') {
    reticulate::py_call(fcc$materialize$synapse_query, post_ids=pyids, synapse_table=synapse_table, ...)
  } else {
    reticulate::py_call(fcc$materialize$synapse_query, pre_ids=pyids, synapse_table=synapse_table, ...)
  }
  fafbseg:::pandas2df(res)
}









