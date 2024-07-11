#' Query banc tables in the CAVE annotation system
#'
#' @param ... Additional arguments passed to \code{\link{flywire_cave_query}}
#'   \code{\link[fafbseg]{flywire_cave_query}}
#' @inheritParams banc_partner_summary
#' @inheritParams fafbseg::flywire_cave_query
#'
#' @return A data.frame
#'
#' @family banc-cave
#' @export
#' @seealso \code{\link[fafbseg]{flywire_cave_query}}
#' @examples
#' \donttest{
#' library(dplyr)
#' cell_info=banc_cave_query('cell_info')
#' cell_info %>%
#'   filter(tag2=='anterior-posterior projection pattern') %>%
#'   count(tag)
#' }
banc_cave_query <- function(table, datastack_name = NULL, live=TRUE, ...) {
  if(is.null(datastack_name)) datastack_name=banc_datastack_name()
  fafbseg::flywire_cave_query(table = table, datastack_name = datastack_name, live=live, ...)
}

#' Low level access to banc's CAVE annotation infrastructure
#'
#' @return A reticulate R object wrapping the python CAVEclient.
#' @export
#'
#' @examples
#' \dontrun{
#' fcc=banc_cave_client()
#' tables=fcc$annotation$get_tables()
#' fcc$materialize$get_table_metadata(tables[1])
#' }
banc_cave_client <- function() {
  with_banc(flywire_cave_client())
}

