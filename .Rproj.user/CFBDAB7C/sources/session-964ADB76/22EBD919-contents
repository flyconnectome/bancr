#' Fetch change log information for one or more neurons
#'
#' @details As of August 2021 this is a simple wrapper of
#'   \code{fabseg::\link{flywire_change_log}}. For now the old (and less
#'   convenient format) available from the zetta API can be obtained with the
#'   private \code{bancr:::banc_change_log_zetta} function.
#'
#' @param x One or more banc ids in any format understandable by
#'   \code{\link[fafbseg]{ngl_segments}}
#' @return a \code{data.frame} See \code{fabseg::\link{flywire_change_log}} for
#'   details
#' @export
#' @importFrom pbapply pbsapply
#' @inheritParams fafbseg::flywire_change_log
#' @examples
#' \donttest{
#' banc_change_log("648518346473954669")
#' }
banc_change_log <- function(x, tz="UTC", filtered=TRUE, OmitFailures=TRUE, ...) {
  x=banc_ids(x)
  with_banc(flywire_change_log(
    x = x,
    tz = tz,
    filtered = filtered,
    OmitFailures = OmitFailures,
    ...
  ))
}

banc_change_log_zetta <- function(x, ...) {
  baseurl=banc_api_url(endpoint = "root/")
  x=fafbseg::ngl_segments(x, as_character = T)
  if(length(x)>1) {
    res=pbapply::pbsapply(x, banc_change_log_zetta, ..., simplify = FALSE)
    user_info=dplyr::bind_rows(sapply(res, "[[", "user_info", simplify = F),
                               .id = 'root_id')
    res2=sapply(res, function(x) x[setdiff(names(x), 'user_info')], simplify = F)
    res2df=do.call(rbind,res2)
    res2df=cbind(data.frame(root_id=rownames(res2df)), res2df)
    attr(user_info, 'summary')=res2df[1:2]
    attr(user_info, 'operations')=res2df[[3]]
    attr(user_info, 'past_ids')=res2df[[4]]
    return(user_info=user_info)
  }
  url=paste0(baseurl, x, "/change_log")
  res=banc_fetch(url, ...)
  res$user_info
  userdf=data.frame(user=names(res$user_info), stringsAsFactors = FALSE)
  userdf$n_splits=sapply(res$user_info, function(x) {splits=x$n_splits; ifelse(is.null(splits), 0L, as.integer(splits))})
  userdf$n_mergers=sapply(res$user_info, function(x) {mergers=x$n_mergers; ifelse(is.null(mergers), 0L, as.integer(mergers))})
  res$operations_ids=scan(text = chartr("[]","  ", res$operations_ids),
                          what=integer(), quiet = T)
  res$user_info=userdf
  res
}
