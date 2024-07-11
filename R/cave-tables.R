#' Read BANC CAVE-tables, good sources of metadata
#'
#' @param rootids #' @param rootids Character vector specifying one or more BANC rootids. As a
#'   convenience this argument is passed to \code{\link{banc_ids}} allowing you
#'   to pass in data.frames, BANC URLs or simple ids.
#' @param nucleus_ids Character vector specifying one or more BANC nucleus ids.
#' @param rawcoords Logical, whether or not yto convert from raw coordinates into nanometers. Default is `FALSE`.
#' @param select A regex term for the name of the table you want
#' @param datastack_name  Defaults to "brain_and_nerve_cord". See https://global.daf-apis.com/info/ for other options.
#' @param table Possible alternative tables for the sort of data frame the function returns. One must be chosen.
#' @param ... Additional arguments passed to
#'   \code{fafbseg::\link{flywire_cave_query}}
#'
#' @return A \code{data.frame} describing a CAVE-table related to the BANC project.
#' In the case of \code{banc_cave_tables}, a vector is returned containing the names of
#' all queryable cave tables.
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
  fac <- flywire_cave_client(datastack_name = datastack_name)
  dsinfo <- fac$info$get_datastack_info()
  if (!is.null(dsinfo$soma_table))
    return(dsinfo$soma_table)
  tt <- fac$annotation$get_tables()
  if(!is.null(select)){
    chosen_tables <- grep(select, tt)
    if (length(chosen_tables) == 0)
      stop(sprintf("I cannot find a '%s' table for datastack: ", select),
           datastack_name, "\nPlease ask for help on #annotation_infrastructure https://flywire-forum.slack.com/archives/C01M4LP2Y2D")
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

#' @rdname banc_cave_tables
#' @export
banc_nuclei <- function (rootids = NULL,
                         nucleus_ids = NULL,
                         table = c("somas_v1a","somas_v1b"),
                         rawcoords = FALSE,
                         ...) {
  table <- match.arg(table)
  if (!is.null(rootids) & !is.null(nucleus_ids))
    stop("You must supply only one of rootids or nucleus_ids!")
  res <- if (is.null(rootids) && is.null(nucleus_ids))
    banc_cave_query(table = table, ...)
  else if (!is.null(rootids)) {
    rootids <- banc_ids(rootids)
    nuclei <- if (length(rootids) < 200)
      banc_cave_query(table =  table, filter_in_dict = list(pt_root_id=rootids),
                      ...)
    else
      banc_cave_query(table =  table, live = F, ...)
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
                                filter_in_dict = list(id=nucleus_ids), ...)
    nuclei %>%
      dplyr::right_join(data.frame(id = as.integer64(nucleus_ids)), by = "id") %>%
      dplyr::select(colnames(nuclei))
  }
  res
  #  apply coordinate transform
  # res <- standard_nuclei(res)
  if (isTRUE(rawcoords))
    res
  else {
    res %>% dplyr::mutate(dplyr::across(dplyr::ends_with("position"), function(x)
      nat::xyzmatrix2str(banc_raw2nm(x))))
  }
}

#' @rdname banc_cave_tables
#' @export
banc_cell_info <- function(rootids = NULL, rawcoords = FALSE){
  table <- "cell_info"
  res <- get_cave_table_data(table)
  if (isTRUE(rawcoords))
    res
  else {
    res %>% dplyr::mutate(dplyr::across(dplyr::ends_with("position"), function(x)
      nat::xyzmatrix2str(banc_raw2nm(x))))
  }
}

#' @rdname banc_cave_tables
#' @export
banc_cell_ids <- function(rootids = NULL){
  get_cave_table_data('cell_ids', rootids)
}

#' @rdname banc_cave_tables
#' @export
banc_neck_connective_neurons <- function(rootids = NULL,
                                         table = c("neck_connective_y92500", "neck_connective_y121000")){
  table <- match.arg(table)
  get_cave_table_data(table, rootids)
}

#' @rdname banc_cave_tables
#' @export
banc_peripheral_nerves <- function(rootids = NULL){
  get_cave_table_data("peripheral_nerves", rootids)
}

#' @rdname banc_cave_tables
#' @export
banc_backbone_proofread <- function(rootids = NULL){
  get_cave_table_data("backbone_proofread", rootids)
}

# hidden
get_cave_table_data <- function(table, rootids = NULL, ...){
  if (!is.null(rootids)) {
    rootids <- flywire_ids(rootids)
    df <- if (length(rootids) < 200) {
      rid <- paste(rootids, collapse = ",")
      ridq <- reticulate::py_eval(sprintf("{\"pt_root_id\": [%s]}",
                                         rid), convert = F)
      fafbseg::flywire_cave_query(table =  table,
                         filter_in_dict = ridq, ...)
    } else {
      fafbseg::flywire_cave_query(table =  table,
                         live = F, ...)
    }
  }else{
    df <- fafbseg::flywire_cave_query(table =  table , ...)
  }
  df
}
