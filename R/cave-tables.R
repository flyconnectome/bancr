#' Read BANC CAVE-tables, good sources of metadata
#'
#' @param rootids #' @param rootids Character vector specifying one or more BANC rootids. As a
#'   convenience this argument is passed to \code{\link{banc_ids}} allowing you
#'   to pass in data.frames, BANC URLs or simple ids.
#' @param nucleus_ids Character vector specifying one or more BANC nucleus ids.
#' @param rawcoords Logical, whether or not to convert from raw coordinates into nanometers. Default is `FALSE`.
#' @param select A regex term for the name of the table you want
#' @param datastack_name  Defaults to "brain_and_nerve_cord". See https://global.daf-apis.com/info/ for other options.
#' @param table Possible alternative tables for the sort of data frame the function returns. One must be chosen.
#' @param ... Additional arguments passed to
#'   \code{fafbseg::\link{flywire_cave_query}}
#'
#' @return A \code{data.frame} describing a CAVE-table related to the BANC project.
#' In the case of \code{banc_cave_tables}, a vector is returned containing the names of
#' all query-able cave tables.
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
  if (isFALSE(rawcoords)) {
    # res <- res %>%
    #   dplyr::mutate(dplyr::across(dplyr::ends_with("position"), function(x)
    #   nat::xyzmatrix2str(banc_raw2nm(x))))
    res$nucleus_position_nm <- apply(banc_raw2nm(res$nucleus_position),1,paste_coords)
    res$nucleus_position_nm  <- gsub("\\(|\\)","",res$nucleus_position_nm)
  }
  res
}

#' @rdname banc_cave_tables
#' @export
#' @importFrom dplyr mutate ends_with across
#' @importFrom nat xyzmatrix2str
banc_cell_info <- function(rootids = NULL, rawcoords = FALSE){
  table <- "cell_info"
  res <- with_banc(get_cave_table_data(table))
  if (isTRUE(rawcoords))
    res
  else {
    res %>% mutate(across(ends_with("position"),
                          function(x) xyzmatrix2str(banc_raw2nm(x))))
  }
}

#' @rdname banc_cave_tables
#' @export
banc_cell_ids <- function(rootids = NULL){
  with_banc(get_cave_table_data('cell_ids', rootids))
}

#' @rdname banc_cave_tables
#' @export
banc_neck_connective_neurons <- function(rootids = NULL,
                                         table = c("neck_connective_y92500", "neck_connective_y121000")){
  table <- match.arg(table)
  with_banc(get_cave_table_data(table, rootids))
}

#' @rdname banc_cave_tables
#' @export
banc_peripheral_nerves <- function(rootids = NULL){
  with_banc(get_cave_table_data("peripheral_nerves", rootids))
}

#' @rdname banc_cave_tables
#' @export
banc_backbone_proofread <- function(rootids = NULL){
  with_banc(get_cave_table_data("backbone_proofread", rootids))
}

# hidden
get_cave_table_data <- function(table, rootids = NULL, ...){
  if (!is.null(rootids)) {
    rootids <- flywire_ids(rootids)
    df <- if (length(rootids) < 200) {
      fafbseg::flywire_cave_query(table =  table,
                         filter_in_dict = list(pt_root_id=rootids), ...)
    } else {
      fafbseg::flywire_cave_query(table =  table, live = TRUE, ...)
    }
  } else {
    df <- fafbseg::flywire_cave_query(table =  table , ...)
  }
  df
}

# hidden
banc_cave_cell_types <- function(user_id = NULL){
  banc.cell.info <- banc_cell_info(rawcoords = TRUE)
  banc.cell.info$pt_position <- sapply(banc.cell.info$pt_position, paste, collapse=", ")
  banc.cell.info.mod <- banc.cell.info %>%
    dplyr::filter(valid == 't') %>%
    dplyr::rowwise() %>%
    dplyr::mutate(pt_position = paste0(pt_position,collapse=",")) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(pt_root_id) %>%
    dplyr::arrange(pt_position, tag2, tag) %>%
    dplyr::mutate(side =  dplyr::case_when(
      grepl("^soma side",tag2) ~ gsub("soma on |soma on ","",tag),
      TRUE ~ NA
    )) %>%
    dplyr::mutate(cell_type = dplyr::case_when(
      grepl("neuron identity", tag2) ~ tag,
      grepl("^DN|^AMMC|^PDN|^LH|^il|^T1|^T5|^T4|^TY4|^IN|^il|^HS|^Mi|^PS|^CB|^FB|^AL|
            ^FET|^bCS|SEZ-NSC-Hugin|^MDN|^OA|^PS|^ovi|giant fiber|^m-NSC|^l-NSC-ITP
            |^OA|^LH|^CSD|^BDN|^AN|^AL|^AV|^AN|^MN|^SA|^Mi|^LH|^L1|^BDN|^LAL",tag) ~ tag,
      TRUE ~ NA
    )) %>%
    dplyr::mutate(user_id = dplyr::case_when(
      !is.na(cell_type) ~ user_id,
      TRUE ~ NA
    )) %>%
    dplyr::mutate(cell_type = gsub("\\\n.*|\\*.*","",cell_type)) %>%
    dplyr::mutate(cell_class = dplyr::case_when(
      grepl("ascending|descending|descending|ascending", tag) ~ tag,
      grepl("sensory neuron|motor neuron|^trachea|^glia|^endocrine", tag) ~ tag,
      grepl("sensory neuron|motor neuron|^trachea|^glia|^endocrine", tag) ~ tag2,
      grepl("motor neuron", tag) ~ "motor",
      grepl("endocrine", tag) ~ "endocrine",
      grepl("central neuron", tag2) ~ tag,
      grepl("^innervates|^intersegmental", tag) ~ tag,
      TRUE ~ NA
    )) %>%
    dplyr::mutate(super_class = dplyr::case_when(
      grepl("ascend", cell_class) ~ "ascending",
      grepl("descend", cell_class) ~ "descending",
      grepl("sensory neuron", cell_class) ~ "sensory",
      grepl("motor neuron", cell_class) ~ "motor",
      grepl("endocrine", cell_class) ~ "endocrine",
      grepl("efferent", cell_class) ~ "efferent",
      grepl("optic", cell_class) ~ "optic",
      grepl("optic", tag) ~ "optic",
      grepl("optic", tag2) ~ "optic",
      grepl("central", cell_class) ~ "intrinsic",
      grepl("glia", cell_class) ~ "glia",
      TRUE ~ NA
    )) %>%
    dplyr::mutate(notes = paste(unique(na.omit(sort(tag))), collapse = ", "),
                  cell_class = paste(unique(na.omit(sort(cell_class))), collapse = ", "),
                  super_class = paste(unique(na.omit(sort(super_class))), collapse = ", "),
                  cell_type = paste(unique(na.omit(sort(cell_type))), collapse = ", "),
                  side = paste(unique(na.omit(sort(side))), collapse = ", "),
                  user_id = paste(unique(na.omit(sort(user_id))), collapse = ", ")) %>%
    dplyr::ungroup() %>%
    dplyr::rename(cell_id = id, root_id = pt_root_id, supervoxel_id = pt_supervoxel_id, position = pt_position) %>%
    dplyr::distinct(root_id, supervoxel_id, side, super_class, cell_class, cell_type, .keep_all = TRUE) %>%
    dplyr::select(cell_id, root_id, supervoxel_id, position, side, super_class, cell_class, cell_type, user_id,notes) %>%
    dplyr::left_join(banc_users %>% dplyr::distinct(pi_lab,cave_id) %>% dplyr::mutate(cave_id=as.character(cave_id)),
                     by=c("user_id"="cave_id")) %>%
    dplyr::rename(cell_type_source = pi_lab)
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






