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
#' \dontrun{
#' library(dplyr)
#' cell_info=banc_cave_query('cell_info')
#' cell_info %>%
#'   filter(tag2=='anterior-posterior projection pattern') %>%
#'   count(tag)
#' }
banc_cave_query <- function(table, live=TRUE, ...) {
  with_banc(fafbseg::flywire_cave_query(table = table, live=live, ...))
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
  with_banc(fafbseg::flywire_cave_client())
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
    warning("Defaulting to latest materialisation version since live=FALSE
",
            "Specify `version='latest' instead to avoid this warning")
    version = flywire_version("latest", datastack_name = datastack_name)
  }
  if (!is.null(timestamp) && !is.null(version))
    stop("You can only supply one of timestamp and materialization version")
  check_package_available("arrow")
  fac = flywire_cave_client(datastack_name = datastack_name)
  if (!requireNamespace("checkmate", quietly = TRUE)) {
    stop("Package 'checkmate' is required for this function. Please install it with: install.packages('checkmate')")
  }
  offset = checkmate::asInt(offset, lower = 0L)
  if (!is.null(limit))
    limit = checkmate::asInt(limit, lower = 0L)
  is_view = table %in% fafbseg:::cave_views(fac)
  version = fafbseg:::flywire_version(version, datastack_name = datastack_name)
  if (!is.null(version)) {
    available = version %in% fafbseg:::flywire_version("available", datastack_name = datastack_name)
    if (!available) {
      if (is_view)
        stop("Sorry! Views only work with unexpired materialisation versions.
",
             "See https://flywire-forum.slack.com/archives/C01M4LP2Y2D/p1697956174773839 for info.")
      timestamp = fafbseg::flywire_timestamp(version, datastack_name = datastack_name,
                                              convert = F)
      message("Materialisation version no longer available. Falling back to (slower) timestamp!")
      if (isFALSE(live))
        live = TRUE
      version = NULL
    }
  }
  now = fafbseg::flywire_timestamp(timestamp = "now", convert = FALSE)
  if(timetravel) {
    timestamp2 = fafbseg::flywire_timestamp(version,
                                             timestamp = timestamp,
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
                "
", "I'm going to try and format your input correctly.")
        filter_regex_dict = list(filter_regex_dict)
        names(filter_regex_dict) = table
      }
    }
  }
  if (!is.null(select_columns)) {
    if (isTRUE(live == 2) && is.character(select_columns)) {
      warning("When live==2 / timetravel=T select_columns should be a list of form: ",
              "`list(<table_name>=c('col1', 'col2'))`", "
",
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
          stop("Sorry! You cannot specify a timestamp when querying a view.
",
               "You can specify older timepoints by using unexpired materialisation versions.
",
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
      warning(paste(pymsg, "
Use fetch_all_rows=T or set an explicit limit to avoid warning!"))
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
           else "
Please review your value of `select_columns`!")
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

# function banc_cave_schema
# function to get schema of specific table
banc_cave_schema <- function(table_name = NULL, datastack_name = NULL)
{
  if (is.null(datastack_name))
    datastack_name = banc_datastack_name()
  fac <- fafbseg::flywire_cave_client(datastack_name = datastack_name)
  dsinfo <- fac$info$get_datastack_info()
  tt <- fac$annotation$get_tables()
  if (!(table_name %in% tt))
  {
    stop(sprintf("I cannot find a '%s' table for datastack: ",table_name), datastack_name)
  }
  else {
    table_meta <- fac$annotation$get_table_metadata(table_name)
    this_schema <- table_meta$schema_type
    return(this_schema)
  }
}

# function banc_cave_valid_schemas
# function that returns all valid schemas
banc_cave_valid_schemas <- function(datastack_name = NULL)
{
  if (is.null(datastack_name))
    datastack_name = banc_datastack_name()
  fac <- fafbseg::flywire_cave_client(datastack_name = datastack_name)

  schema_types <- fac$schema$get_schemas()

  return(schema_types)
}


# function banc_annotate_bound_double_tag_user
# function to upload new annotations to CAVE tables of schema bound_double_tag_user

# hidden
banc_annotate_bound_double_tag_user <- function (data, units = c("raw", "nm"),
                                                 table_name = NULL, datastack_name = NULL, use_admin_creds = FALSE)
{
  if (is.null(datastack_name))
    datastack_name = banc_datastack_name()

  # check that table to write to exists
  all_tables <- banc_cave_tables(datastack_name = datastack_name)
  if (!(table_name %in% all_tables)) {
    stop(sprintf("I cannot find a '%s' table for datastack: ", table_name),
         datastack_name)
  }

  # check that table is of bound_double_tag_user schema
  this_schema <- banc_cave_schema(table_name = table_name,
                                  datastack_name = datastack_name)
  if (this_schema != "bound_double_tag_user") {
    stop(sprintf("'%s' is not of schema bound_double_tag_user", table_name))
  }

  # check that data is a dataframe that contains the columns pt_position,
  #  tag, tag2, and user_id
  # Check if data is a dataframe
  if (!is.data.frame(data)) {
    stop("data must be a dataframe")
  }

  # Check for required columns
  required_columns <- c("pt_position", "tag", "tag2", "user_id")
  missing_columns <- required_columns[!required_columns %in% colnames(data)]

  if (length(missing_columns) > 0) {
    stop(paste("data is missing required columns:",
               paste(missing_columns, collapse = ", ")))
  }

  # validate positions and convert to dataframe of X,Y,Z
  positions <- banc_validate_positions(positions = data$pt_position,
                                               units = units)

  # check for valid IDs; remove rows where pt_position corresponds to invalid ID
  valid_ids = banc_xyz2id(positions, rawcoords = TRUE)
  valid_ids_not_0 = valid_ids[valid_ids != "0"]
  positions = positions[valid_ids != "0", ]
  data = data[valid_ids != "0", ]
  if (sum(valid_ids == "0")) {
    warning("given positions with invalid root_id: ",
            sum(valid_ids == "0"))
  }
  if (!nrow(positions)) {
    stop("no valid positions given")
  }

  # check CAVE
  cavec = fafbseg:::check_cave()
  # initialize python packages
  np = reticulate::import("numpy")
  pd = reticulate::import("pandas")

  # instantiate CAVE client
  if (use_admin_creds) {
    client = banc_service_account(datastack_name)
  }
  else {
    client = fafbseg::flywire_cave_client(datastack_name = datastack_name)
  }

  stage <- client$annotation$stage_annotations(table_name)

  result_ind <- integer(0)
  if (!requireNamespace("progress", quietly = TRUE)) {
    stop("Package 'progress' is required for this function. Please install it with: install.packages('progress')")
  }
  pb <- progress::progress_bar$new(format = "[:bar] :percent | ETA: :eta | :current/:total rows | Elapsed: :elapsedfull",
                                   total = nrow(data), clear = FALSE, width = 80)
  for (i in 1:nrow(data)) {
    pb$tick()
    this_pos <- unlist(positions[i, ])
    this_tag <- data$tag[i]
    this_tag2 <- data$tag2[i]
    this_user_id <- data$user_id[i]
    stage$add(valid = TRUE, pt_position = np$array(this_pos),
              tag = this_tag, tag2 = this_tag2,
              user_id = as.integer(this_user_id))
    this_result <- client$annotation$upload_staged_annotations(stage)
    result_ind <- c(result_ind, this_result)
    stage$clear_annotations()
  }
  pause_seconds <- nrow(data) * 0.2
  Sys.sleep(pause_seconds)

  annotations <- banc_cave_query(table_name, live = 2)
  annotations.new <- annotations %>%
    dplyr::filter(id %in% result_ind)
  cat("added ", nrow(annotations.new), "new annotations to ", table_name, "
")
  return(annotations.new)
}

# function banc_cave_new_table
# function to create new CAVE table, non reference types only
banc_cave_new_table <- function(table_name, schema_name, description = NULL,
                                voxel_resolution = banc_voxdims(),
                                write_permission = c("PRIVATE", "GROUP", "PUBLIC"),
                                read_permission = c("PRIVATE", "GROUP", "PUBLIC"),
                                user_id = NULL, datastack_name = NULL)
{
  if (is.null(datastack_name))
    datastack_name = banc_datastack_name()

  # check that table_name does not already exist
  all_tables <- banc_cave_tables(datastack_name = datastack_name)
  if ((table_name %in% all_tables)) {
    stop(sprintf("'%s' table already exists for datastack: ", table_name),
         datastack_name)
  }

  # check that schema_name is valid and not of reference kind
  all_schemas <- banc_cave_valid_schemas()
  all_nonRef_schemas <- all_schemas[!grepl("reference", all_schemas)]

  if (!(schema_name %in% all_nonRef_schemas)) {
    stop(sprintf("'%s' is not a valid schema type", schema_name))
  }

  # instantiate CAVE client
  client = fafbseg::flywire_cave_client(datastack_name = datastack_name)

  # get read and write permissions
  write_permission <- match.arg(write_permission)
  read_permission <- match.arg(read_permission)

  # create table
  if (is.null(user_id)) {
    client$annotation$create_table(table_name = table_name, schema_name = schema_name,
                                   description = description, voxel_resolution = voxel_resolution,
                                   write_permission = write_permission,
                                   read_permission = read_permission)
  }
  else {
    client$annotation$create_table(table_name = table_name, schema_name = schema_name,
                                   description = description, voxel_resolution = voxel_resolution,
                                   user_id = user_id, write_permission = write_permission,
                                   read_permission = read_permission)
  }

  # check that table was created
  all_tables <- banc_cave_tables(datastack_name = datastack_name)
  if ((table_name %in% all_tables)) {
    cat(table_name, "created 
")
  }
}


# function banc_delete_cave_table
# deletes named CAVE table. Only works on CAVE tables you have created
banc_delete_cave_table <- function(table_name = NULL, datastack_name = NULL)
{
  if (is.null(datastack_name))
    datastack_name = banc_datastack_name()

  # check that table to delete to exists
  all_tables <- banc_cave_tables(datastack_name = datastack_name)
  if (!(table_name %in% all_tables)) {
    stop(sprintf("I cannot find a '%s' table for datastack: ", table_name),
         datastack_name)
  }

  # have user confirm deletion
  cat(sprintf("Are you sure you want to delete '%s' (y/n)? ", table_name))
  first_confirm <- tolower(readLines(n = 1))

  if (first_confirm == "y") {
    cat(sprintf("To confirm deletion, please type the exact table name '%s': ", table_name))
    second_confirm <- readLines(n = 1)

    if (second_confirm == table_name) {
      # Execute the deletion code here
      cat(sprintf("Deleting table '%s'...
", table_name))

      client <- fafbseg::flywire_cave_client(datastack_name = datastack_name)

      client$annotation$delete_table(table_name = table_name)

      all_tables <- banc_cave_tables(datastack_name = datastack_name)
      if (!(table_name %in% all_tables)) {
        cat("Table deleted successfully.
")
      }
      else {
        cat("Table marked for deletion.
")
      }
    } else {
      cat("Table name confirmation failed. Deletion cancelled.
")
    }
  } else {
    cat("Deletion cancelled.
")
  }
}


