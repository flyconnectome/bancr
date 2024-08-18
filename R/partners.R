#' Summarise the connectivity of banc neurons
#'
#' @details note that the rootids you pass in must be up to date. See example.
#'
#' @param rootids Character vector specifying one or more BANC rootids. As a
#'   convenience this argument is passed to \code{\link{banc_ids}} allowing you
#'   to pass in data.frames, BANC URLs or simple ids.
#' @param datastack_name An optional CAVE \code{datastack_name}. If unset a
#'   sensible default is chosen.
#' @inheritParams fafbseg::flywire_partner_summary
#'
#' @return a data.frame
#' @seealso \code{\link{flywire_partner_summary}}, \code{\link{banc_latestid}}
#' @export
#'
#' @examples
#' \dontrun{
#' # NB id must be up to date
#' sample_id=banc_latestid("720575941478275714")
#' head(banc_partner_summary(sample_id))
#' head(banc_partner_summary(sample_id, partners='inputs'))
#' # get the latest id for an outdate
#' banc_partner_summary(banc_latestid("720575941478275714"))
#'
#' ## open banc/flywire scene containing top partners
#' library(dplyr)
#' banc_partner_summary(banc_latestid("720575941478275714"), partners='inputs') %>%
#'   slice_max(weight, n = 20) %>%
#'   banc_scene(open=TRUE)
#' }
#' @rdname banc_partner_summary
banc_partner_summary <- function(rootids,
                                 partners = c("outputs", "inputs"),
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
      method = "cave",
      datastack_name = datastack_name,
      remove_autapses = remove_autapses,
      cleft.threshold = cleft.threshold,
      ...
    )
  )
}

### TODO ####
banc_datastack_name <- memoise::memoise(function() {
  banc_name <- tryCatch({
    cac=fafbseg::flywire_cave_client(NULL)
    datastacks=cac$info$get_datastacks()
    seldatastack=grep("brain_and_nerve_cord*", datastacks, value = T)
    if(length(seldatastack)==0)
      stop("Could not identify a banc production datastack amongst: ",
           paste(datastacks, collapse=','),
           "\nHave you been granted access to banc production?")
    if(length(seldatastack)>1)
      warning("Multiple banc datastacks available; ",
              paste(seldatastack, collapse = ","),"\n",
              "choosing: ", seldatastack[1])
    seldatastack[1]
  }, error = function(e){
    warning("using default setting")
    "brain_and_nerve_cord"
  })
  banc_name
})

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
banc_partners <- function(rootids, partners=c("input", "output"), ...) {
  partners=match.arg(partners)
  rootids=banc_ids(rootids)
  fcc=banc_cave_client()
  pyids=fafbseg:::rids2pyint(rootids)
  res=if(partners=='input') {
    reticulate::py_call(fcc$materialize$synapse_query, post_ids=pyids, ...)
  } else {
    reticulate::py_call(fcc$materialize$synapse_query, pre_ids=pyids, ...)
  }
  fafbseg:::pandas2df(res)
}









