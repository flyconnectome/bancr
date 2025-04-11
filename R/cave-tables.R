#' Read BANC CAVE-tables, good sources of metadata
#'
#' @param rootids #' @param rootids Character vector specifying one or more BANC rootids. As a
#'   convenience this argument is passed to \code{\link{banc_ids}} allowing you
#'   to pass in data.frames, BANC URLs or simple ids.
#' @param nucleus_ids Character vector specifying one or more BANC nucleus ids.
#' @param rawcoords Logical, whether or not to convert from raw coordinates into nanometers. Default is `FALSE`.
#' @param select A regex term for the name of the table you want
#' @param datastack_name  Defaults to "brain_and_nerve_cord". See https://global.daf-apis.com/info/ for other options.
#' @param table Character, possible alternative tables for the sort of data frame the function returns. One must be chosen.
#' @param edgelist_view Character, name of prepared CAVE view that computes the proofread-neuron edgelist.
#' @param ... Additional arguments passed to
#'   \code{fafbseg::\link{flywire_cave_query}} or \code{bancr:::get_cave_table_data}.
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
banc_edgelist <- function(edgelist_view = c("synapses_250226_backbone_proofread_and_peripheral_nerves_counts",
                                            "synapses_250226_backbone_proofread_counts",
                                            "synapses_v1_backbone_proofread_counts"),
                          ...){
  edgelist_view <- match.arg(edgelist_view)
  el <- with_banc(cave_view_query(edgelist_view, fetch_all_rows= TRUE, ...))
  el <- el %>%
    dplyr::arrange(dplyr::desc(n))
  el
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
           datastack_name, "\nPlease ask for help on #annotation_infrastructure https://flywire-forum.slack.com/archives/C01M4LP2Y2D")
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
banc_cell_info <- function(rootids = NULL, rawcoords = FALSE, ...){
  table <- "cell_info"
  res <- with_banc(get_cave_table_data(table, ...))
  if (isTRUE(rawcoords))
    res
  else {
    res %>% mutate(across(ends_with("position"),
                          function(x) xyzmatrix2str(banc_raw2nm(x))))
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
banc_cave_cell_types <- function(cave_id = NULL, invert = FALSE, ...){
  banc.cell.info <- banc_cell_info(rawcoords = TRUE, ...)
  if(!is.null(cave_id)){
    if(invert){
      banc.cell.info <- banc.cell.info %>%
        dplyr::filter(!user_id %in% cave_id)
    }else{
      banc.cell.info <- banc.cell.info %>%
        dplyr::filter(user_id %in% cave_id)
    }
  }
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
      !grepl(",",tag) ~ tag,
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


# hidden
cave_view_query <- function(table,
                            datastack_name = getOption("fafbseg.cave.datastack_name","flywire_fafb_production"),
                            version = NULL,
                            timestamp = NULL,
                            live = is.null(version),
                            timetravel = FALSE,
                            filter_in_dict = NULL, filter_out_dict = NULL, filter_regex_dict = NULL,
                            filter_equal_dict=NULL, filter_greater_dict=NULL, filter_less_dict=NULL,
                            filter_greater_equal_dict=NULL, filter_less_equal_dict=NULL, filter_spatial_dict=NULL,
                            select_columns = NULL, offset = 0L, limit = NULL, fetch_all_rows = FALSE,
                            ...) {
  if (isTRUE(live) && !is.null(version))
    warning("live=TRUE so ignoring materialization version")
  if (is.null(live) && !is.null(timestamp))
    live = TRUE
  if (isFALSE(live) && is.null(version)) {
    warning("Defaulting to latest materialisation version since live=FALSE\n",
            "Specify `version='latest' instead to avoid this warning")
    version = flywire_version("latest", datastack_name = datastack_name)
  }
  if (!is.null(timestamp) && !is.null(version))
    stop("You can only supply one of timestamp and materialization version")
  check_package_available("arrow")
  fac = flywire_cave_client(datastack_name = datastack_name)
  offset = checkmate::asInt(offset, lower = 0L)
  if (!is.null(limit))
    limit = checkmate::asInt(limit, lower = 0L)
  is_view = table %in% fafbseg:::cave_views(fac)
  version = fafbseg:::flywire_version(version, datastack_name = datastack_name)
  if (!is.null(version)) {
    available = version %in% fafbseg:::flywire_version("available", datastack_name = datastack_name)
    if (!available) {
      if (is_view)
        stop("Sorry! Views only work with unexpired materialisation versions.\n",
             "See https://flywire-forum.slack.com/archives/C01M4LP2Y2D/p1697956174773839 for info.")
      timestamp = fafbseg:::flywire_timestamp(version, datastack_name = datastack_name,
                                    convert = F)
      message("Materialisation version no longer available. Falling back to (slower) timestamp!")
      if (isFALSE(live))
        live = TRUE
      version = NULL
    }
  }
  now = fafbseg:::flywire_timestamp(timestamp = "now", convert = FALSE)
  if(timetravel) {
    timestamp2 = fafbseg:::flywire_timestamp(version, timestamp = timestamp,
                                   datastack_name = datastack_name)
    timestamp = now
    version = NULL
    live = 2L
  }
  filter_in_dict = fafbseg:::cavedict_rtopy(filter_in_dict, wrap_table = if (isTRUE(live == 2))
    table
    else NULL)
  filter_out_dict = fafbseg:::cavedict_rtopy(filter_out_dict)
  filter_equal_dict = fafbseg:::cavedict_rtopy(filter_equal_dict)
  filter_greater_dict = fafbseg:::cavedict_rtopy(filter_greater_dict)
  filter_less_dict = fafbseg:::cavedict_rtopy(filter_less_dict)
  filter_greater_equal_dict = fafbseg:::cavedict_rtopy(filter_greater_equal_dict)
  filter_less_equal_dict = fafbseg:::cavedict_rtopy(filter_less_equal_dict)
  filter_spatial_dict = fafbseg:::cavedict_rtopy(filter_spatial_dict)
  if (!is.null(filter_regex_dict)) {
    was_char = is.character(filter_regex_dict)
    if (was_char) {
      filter_regex_dict = as.list(filter_regex_dict)
      if (isTRUE(live == 2)) {
        warning("When live==2 / timetravel=T filter_regex_dict should be a list of form: ",
                "`list(<table_name>=c(<colname>='<regex>'))`",
                "\n", "I'm going to try and format your input correctly.")
        filter_regex_dict = list(filter_regex_dict)
        names(filter_regex_dict) = table
      }
    }
  }
  if (!is.null(select_columns)) {
    if (isTRUE(live == 2) && is.character(select_columns)) {
      warning("When live==2 / timetravel=T select_columns should be a list of form: ",
              "`list(<table_name>=c('col1', 'col2'))`", "\n",
              "I'm going to try and format your input correctly.")
      select_columns = list(select_columns)
      names(select_columns) = table
    }
  }
  annotdfs = list()
  while (offset >= 0) {
    pymsg <- reticulate::py_capture_output({
      annotdf <- if (is_view) {
        if (!is.null(timestamp))
          stop("Sorry! You cannot specify a timestamp when querying a view.\n",
               "You can specify older timepoints by using unexpired materialisation versions.\n",
               "See https://flywire-forum.slack.com/archives/C01M4LP2Y2D/p1697956174773839 for info.")
        reticulate::py_call(fac$materialize$query_view,
                            view_name = table, materialization_version = version,
                            filter_in_dict = filter_in_dict, filter_out_dict = filter_out_dict,
                            filter_equal_dict=filter_equal_dict,
                            #filter_greater_dict=filter_equal_dict,
                            #filter_less_dict=filter_less_dict,filter_greater_equal_dict=filter_greater_equal_dict,
                            #filter_less_equal_dict=filter_less_equal_dict, filter_spatial_dict=filter_spatial_dict,
                            filter_regex_dict = filter_regex_dict, select_columns = select_columns,
                            offset = offset, limit = limit, ...)
      }else{
        stop("Sorry, no valid view specified")
      }
      annotdf <- fafbseg:::pandas2df(annotdf)
    })
    annotdfs[[length(annotdfs) + 1]] = annotdf
    limited_query = isTRUE(grepl("Limited query to", pymsg))
    if (limited_query && is.null(limit) && !fetch_all_rows)
      warning(paste(pymsg, "\nUse fetch_all_rows=T or set an explicit limit to avoid warning!"))
    else if (!limited_query && nzchar(pymsg)) {
      warning(pymsg)
    }
    offset <- if (fetch_all_rows && limited_query)
      offset + nrow(annotdf)
    else -1L
  }
  res <- if (length(annotdfs) == 1)
    annotdfs[[1]]
  else dplyr::bind_rows(annotdfs)
  if (timetravel) {
    if (!all(c("pt_supervoxel_id", "pt_root_id") %in% colnames(res)))
      stop("Sorry I do not know how to time travel dataframes without `pt_supervoxel_id`, `pt_root_id` columns!",
           if (is.null(select_columns))
             ""
           else "\nPlease review your value of `select_columns`!")
    res$pt_root_id = flywire_updateids(res$pt_root_id, svids = res$pt_supervoxel_id,
                                       timestamp = timestamp2, cache = T, Verbose = F)
  }
  res
}


# hidden
banc_service_account <- function(datastack_name=banc_datastack_name()){
  if(is.null(datastack_name)){
    datastack_name=banc_datastack_name()
  }
  cavec = fafbseg:::check_cave()
  client = try(cavec$CAVEclient(datastack_name=datastack_name, auth_token_key='banc_service_account'))
  if (inherits(client, "try-error")) {
    stop("There seems to be a problem connecting to datastack as banc_service_account: ",
         datastack_name)
  }
  client
}

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
  if(!is.null(nrow(positions))){
    positions <- unlist(c(positions))
  }else if(nrow(positions)==1){
    positions <- unlist(c(positions))
  }

  # convert
  if(units=="nm"){
    positions = banc_nm2raw(positions)
  }
  positions
}

# hidden
#' @examples
#' # Add an annotation to a point in raw voxel space
#' banc_annotate_backbone_proofread(c(117105, 240526, 5122), user_id = 355, units = "raw")
#'
#' # Add an annotation to a point in nm
#' banc_annotate_backbone_proofread(c(468420, 962104 ,230490), user_id = 355, units = "nm")
#'
#' # deannotate a point, only from points added with given user_id. Use user_id = NULL to remove from full pool
#' banc_deannotate_backbone_proofread(c(468420, 962104, 230490), user_id = 355, units = "nm")

banc_annotate_backbone_proofread <- function(positions,
                                             user_id,
                                             units = c("raw","nm"),
                                             proofread = TRUE,
                                             datastack_name = NULL){

  # Validate positions
  positions <- banc_validate_positions(positions=positions, units=units)

  # get table
  cavec = fafbseg:::check_cave()
  np = reticulate::import("numpy")
  pd = reticulate::import("pandas")
  client = banc_service_account(datastack_name)

  # Current state
  annotations <- banc_backbone_proofread(live=2) %>%
    dplyr::filter(proofread==eval(proofread))
  if(!nrow(annotations)){
    stop("no annotations collected")
  }
  curr.positions <- do.call(rbind,annotations$pt_position)
  curr.positions <- as.data.frame(curr.positions)
  colnames(curr.positions) <- c("X","Y","Z")
  curr.positions$id <- annotations$id

  # Stage
  stage <- client$annotation$stage_annotations('backbone_proofread')
  if(is.data.frame(positions)){

    # Validate root IDs
    positions.orig <- positions
    positions <- dplyr::anti_join(positions, as.data.frame(curr.positions), by = c("X","Y","Z"))
    cat("given positions already in back_bone proofread:",nrow(positions.orig)-nrow(positions),"\n")
    if(!nrow(positions)){
      stop("all positions already marked:", nrow(positions.orig))
    }
    valid_ids = banc_xyz2id(positions, rawcoords = TRUE)
    valid_ids_not_0 = valid_ids[valid_ids!="0"]
    positions = positions[valid_ids!="0",]
    if(sum(valid_ids=="0")){
      warning("given positions with invalid root_id: ", sum(valid_ids=="0"))
    }
    if(!nrow(positions)){
      stop("no valid positions given")
    }

    # Create a pandas dataframe with all required columns
    pos_list <- lapply(1:nrow(positions), function(i) {
      np$array(positions[i,])
    })
    n_points <- nrow(positions)
    annotation_df <- pd$DataFrame(data = list(
      pt_position = reticulate::r_to_py(pos_list),
      valid = rep(TRUE, n_points),
      user_id = rep(user_id, n_points),
      valid_id = as.numeric(valid_ids_not_0),
      proofread = rep(proofread, n_points)
    ))
    stage$add_dataframe(annotation_df)
  }else{
    matching_rows <- which(apply(curr.positions[,1:3], 1, function(row) all(row == positions)))
    point_exists <- length(matching_rows) > 0
    if(point_exists){
      stop("given position already marked")
    }
    valid_id = banc_xyz2id(positions, rawcoords = TRUE)
    if(valid_id=="0"){
      stop("given position does not return a valid root_id")
    }
    stage$add(valid = TRUE,
              pt_position = np$array(positions),
              user_id = as.integer(user_id),
              valid_id = as.numeric(valid_id),
              proofread = proofread)
  }

  # Upload the staged annotations
  result <- client$annotation$upload_staged_annotations(stage)

  # Read table and check annotations are added
  annotations <- banc_backbone_proofread(live=2)
  annotations.new <- annotations %>%
    dplyr::filter(id %in% result)
  cat("annotated", nrow(annotations.new), "entities with backbone proofread:", proofread,"\n")
  return(annotations.new)
}

# hidden
banc_deannotate_backbone_proofread <- function(positions,
                                               user_id = NULL,
                                               units = c("raw","nm")){

  # Validate positions
  positions <- banc_validate_positions(positions=positions, units=units)

  # Read table and check annotations are added
  annotations <- banc_backbone_proofread(live=2)
  if(!is.null(user_id)){
    annotations <- annotations %>%
      dplyr::filter(user_id %in% !!user_id)
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
  cat("pt_positions in backbone_proofread match to", length(annotation_ids), "given points \n")
  if(length(annotation_ids)){
    # get table
    client <- banc_service_account(datastack_name=datastack_name)

    # Delete specified IDs
    result <- client$annotation$delete_annotation("backbone_proofread", annotation_ids)

    # Read table and check annotations are added
    annotations <- banc_backbone_proofread(live=2)
    annotations.new <- annotations %>%
      dplyr::filter(id %in% annotation_ids)
    if(nrow(annotations.new)){
      warning('not all given positions removed from : missing annotation_ids')
    }
    cat("deannotated", length(result), "entities, valid set to FALSE \n")
    return(result)
  }else{
    invisible()
  }
}









