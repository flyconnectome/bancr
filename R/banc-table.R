#' @title Read and write to the seatable for draft BANC annotations
#'
#' @description These functions use the logic and wrap some code
#' from the `flytable_.*` functions in the `fafbseg` R package.
#' \code{banctable_set_token} will obtain and store a permanent
#'   seatable user-level API token.
#'   \code{banctable_query} performs a SQL query against a banctable
#'   database. You can omit the \code{base} argument unless you have tables of
#'   the same name in different bases.
#'   \code{banctable_base} returns a \code{base} object (equivalent to
#'   a mysql database) which allows you to access one or more tables, logging in
#'   to the service if necessary. The returned base object give you full access
#'   to the Python
#'   \href{https://seatable.github.io/seatable-scripts/python/base/}{\code{Base}}
#'    API allowing a range of row/column manipulations.
#'    \code{banctable_update_rows} updates existing rows in a table, returning TRUE on success.
#'
#' @param sql A SQL query string. See examples and
#'   \href{https://seatable.github.io/seatable-scripts/python/query/}{seatable
#'   docs}.
#' @param limit An optional limit, which only applies if you do not specify a
#'   limit directly in the \code{sql} query. By default seatable limits SQL
#'   queries to 100 rows. We increase the limit to 100000 rows by default.
#' @param convert Expert use only: Whether or not to allow the Python seatable
#'   module to process raw output from the database. This is is principally for
#'   debugging purposes. NB this imposes a requirement of seatable_api >=2.4.0.
#' @param python Logical. Whether to return a Python pandas DataFrame. The default of FALSE returns an R data.frame
#' @param base Character vector specifying the \code{base}
#' @param table Character vector specifying a table foe which you want a
#'   \code{base} object.
#' @param workspace_id A numeric id specifying the workspace. Advanced use only
#   since we can normally figure this out from \code{base_name}.
# @param cached Whether to use a cached base object
# @param token normally retrieved from \code{BANCTABLE_TOKEN} environment variable.
#' @param user,pwd banctable user and password used by \code{banctable_set_token}
#'   to obtain a token
#' @param url Optional URL to the server
#' @param ac A seatable connection object as returned by \code{banctable_login}.
#' @param df A data.frame containing the data to upload including an `_id`
#' column that can identify each row in the remote table.
#' @param append_allowed Logical. Whether rows without row identifiers can be appended.
#' @param chunksize To split large requests into smaller ones with max this many rows.
#' @param token_name The name of the token in your .Renviron file, should be \code{BANCTABLE_TOKEN}.
#' @param where Optional SQL-like where clause to filter rows (default: NULL moves all rows)
#' @param bigdata logical, if `TRUE` new rows are added to the bigdata archive rather than the 'normal' seatable.
#' @param invert whether to send the specified rows (`where`) to big data storage (`FALSE`) or from storage to the 'normal' table (`FALSE`.)
#' @param table.max the maximum number of rows to read from the seatable at one time, which is capped at 10000L by seatable.
#' @param row_ids Character, seatable row IDs
#' @param ... Additional arguments passed to the underlying parallel processing functions which might include cl=2 to specify a number of parallel jobs to run.
#' @param retries if a request to the seatable API fails, the number of times to re-try with a 0.1 second pause.
#' @return a \code{data.frame} of results. There should be 0 rows if no rows
#'   matched query.
#'
#' @seealso \code{fafbseg::\link{flytable_query}}
#' @examples
#' \dontrun{
#' # Do this once
#' banctable_set_token(user="MY_EMAIL_FOR_SEATABLE.com",
#'                     pwd="MY_SEATABLE_PASSWORD",
#'                     url="https://cloud.seatable.io/")
#'
#' # Thereafter:
#' banc.meta <- banctable_query()
#' }
#' @export
#' @rdname banctable_query
banctable_query <- function (sql = "SELECT * FROM banc_meta",
                             limit = 200000L,
                             base = NULL,
                             python = FALSE,
                             convert = TRUE,
                             ac = NULL,
                             token_name = "BANCTABLE_TOKEN",
                             workspace_id = "57832",
                             retries = 10,
                             table.max = 10000L){
  if(is.null(ac)) ac <- banctable_login(token_name=token_name)
  table.max <- 10000L
  if(limit>table.max){
    offset <- 0
    df <- data.frame()
    while(offset<limit){
      cat("reading from row: ", offset,"
")
      sql.new <- sprintf("%s LIMIT %d OFFSET %d", sql, table.max, offset)
      tries <- retries
      bc <- data.frame()
      while(tries>0&&!nrow(bc)){
        bc <- banctable_query(sql=sql.new,
                              limit=FALSE,
                              base=base,
                              python=python,
                              convert=convert,
                              ac=ac,
                              token_name=token_name,
                              workspace_id=workspace_id)
        tries <- tries -1
        Sys.sleep(0.1)
      }
      df <- rbind(df,bc)
      offset <- offset+nrow(bc)
      if(!length(bc)|nrow(bc)<table.max){
        cat("read rows: ",nrow(df), " read columns:", ncol(df),"
")
        return(df)
      }
    }
    cat("read rows: ",nrow(df), " read columns:", ncol(df),"
")
    return(df)
  }
  if (!requireNamespace("checkmate", quietly = TRUE)) {
    stop("Package 'checkmate' is required for this function. Please install it with: install.packages('checkmate')")
  }
  if (!requireNamespace("stringr", quietly = TRUE)) {
    stop("Package 'stringr' is required for this function. Please install it with: install.packages('stringr')")
  }
  checkmate::assert_character(sql, len = 1, pattern = "select",
                              ignore.case = T)
  res = stringr::str_match(sql, stringr::regex("\\s+FROM\\s+[']{0,1}([^, ']+).*",
                                               ignore_case = T))
  if (any(is.na(res)[, 2]))
    stop("Cannot identify a table name in your sql statement!
")
  table = res[, 2]
  if (is.null(base)) {
    base = try(banctable_base(table = table, workspace_id = workspace_id, token_name = token_name))
    if (inherits(base, "try-error"))
      stop("I inferred table_name: ", table, " from your SQL query but couldn't connect to a base with this table!")
  }
  else if (is.character(base))
    base = banctable_base(base_name = base, workspace_id = workspace_id, token_name = token_name)
  if (!isTRUE(grepl("\\s+limit\\s+\\d+", sql)) && !isFALSE(limit)) {
    if (!is.finite(limit))
      limit = .Machine$integer.max
    sql = paste(sql, "LIMIT", limit)
  }
  pyout <- reticulate::py_capture_output(ll <- try(reticulate::py_call(base$query,
                                                                       sql, convert = convert), silent = T))
  if (inherits(ll, "try-error")) {
    warning(paste("No rows returned by banctable", pyout,
                  collapse = "
"))
    return(NULL)
  }
  pd = reticulate::import("pandas")
  reticulate::py_capture_output(pdd <- reticulate::py_call(pd$DataFrame,
                                                           ll))
  if (python)
    pdd
  else {
    colinfo = fafbseg::flytable_columns(table, base)
    df = banctable2df(fafbseg:::pandas2df(pdd, use_arrow = F), tidf = colinfo)
    fields = fafbseg:::sql2fields(sql)
    if (length(fields) == 1 && fields == "*") {
      toorder = intersect(colinfo$name, colnames(df))
    }
    else {
      toorder = intersect(fafbseg:::sql2fields(sql), colnames(df))
    }
    rest = setdiff(colnames(df), toorder)
    df[c(toorder, rest)]
  }
}

#' @export
#' @rdname banctable_query
banctable_set_token <- function(user,
                                pwd,
                                url = "https://cloud.seatable.io/",
                                token_name = "BANCTABLE_TOKEN"){
  st <- fafbseg:::check_seatable()
  ac <- reticulate::py_call(st$Account, login_name = user,
                            password = pwd, server_url = url)
  ac$auth()
  Sys.setenv(banctable_TOKEN = ac$token)
  cat(token_name,"='", ac$token, "'
", sep = "", append = TRUE,
      file = path.expand("~/.Renviron"))
  return(invisible(NULL))
}

#' @export
#' @rdname banctable_query
banctable_login <- function(url = "https://cloud.seatable.io/",
                            token_name = "BANCTABLE_TOKEN"){
  token = Sys.getenv(token_name, unset = NA_character_)
  fafbseg::flytable_login(url=url, token=token)
}


#' @export
#' @rdname banctable_query
banctable_update_rows <- function (df,
                                   table,
                                   base = NULL,
                                   append_allowed = FALSE,
                                   chunksize = 1000L,
                                   workspace_id = "57832",
                                   token_name = "BANCTABLE_TOKEN",
                                   ...) {
  df <- as.data.frame(df)
  if (is.character(base) || is.null(base))
    base = banctable_base(base_name = base, table = table, workspace_id = workspace_id, token_name = token_name)
  nx <- nrow(df)
  if (!isTRUE(nx > 0)) {
    warning("No rows to update in `df`!")
    return(TRUE)
  }
  tablecols = fafbseg::flytable_columns(table,base)
  df = fafbseg:::df2flytable(df, append = ifelse(append_allowed, NA,FALSE))
  newrows = is.na(df[["row_id"]])
  if (any(newrows)) {
    stop("Adding new rows not yet implemented")
    banctable_append_rows(df[newrows, , drop = FALSE], table = table,
                         base = base, chunksize = chunksize, ...)
    df = df[!newrows, , drop = FALSE]
    nx = nrow(df)
  }
  if (!isTRUE(nx > 0))
    return(TRUE)
  if (nx > chunksize) {
    nchunks = ceiling(nx/chunksize)
    chunkids = rep(seq_len(nchunks), rep(chunksize, nchunks))[seq_len(nx)]
    chunks = split(df, chunkids)
    if (!requireNamespace("pbapply", quietly = TRUE)) {
      stop("Package pbapply is required for this function. Please install it with: install.packages('pbapply')")
    }
    oks = pbapply::pbsapply(chunks, banctable_update_rows,
                            table = table, base = base, chunksize = Inf, append_allowed = FALSE,
                            ...)
    return(all(oks))
  }
  multi = tablecols$name[tablecols$type=="multiple-select"]
  if(length(multi)){
    i = intersect(colnames(df),multi)
    if(length(i)){
      for(j in i){
        df[[j]][is.na(df[[j]])] = ''
        l = sapply(df[[j]], strsplit, split = ",|, ")
        l = unname(l)
        df[[j]] = l
      }
    }
  }
  pyl = banc_df2updatepayload(df, via_json = TRUE)
  res = base$batch_update_rows(table_name = table, rows_data = pyl)
  ok = isTRUE(all.equal(res, list(success = TRUE)))
  return(ok)
}

# hidden
banctable_base <- function(base_name = "banc_meta",
                            table = NULL,
                            url = "https://cloud.seatable.io/",
                            token_name = "BANCTABLE_TOKEN",
                            workspace_id = "57832",
                            cached = TRUE,
                            ac = NULL) {
  if(is.null(ac)) ac <- banctable_login(token_name=token_name)
  if (!cached) {
    if (requireNamespace("memoise", quietly = TRUE)) {
      memoise::forget(banctable_base_impl)
    }
  }
  base = try({
    banctable_base_impl(table = table, base_name = base_name,
                       url = url, workspace_id = workspace_id)
  }, silent = TRUE)
  stale_token <- isTRUE(try(difftime(base$jwt_exp, Sys.time(),
                                     units = "hours") < 1, silent = T))
  retry = (cached && inherits(base, "try-error")) || stale_token
  if (!retry)
    return(base)
  if (requireNamespace("memoise", quietly = TRUE)) {
    memoise::forget(banctable_base_impl)
  }
  banctable_base_impl(table = table,
                      base_name = base_name,
                      url = url,
                      workspace_id = workspace_id,
                      token_name = token_name)
}

# hidden
banctable_base_impl <- function (base_name = "banc_meta",
                                 table = NULL,
                                 url = "https://cloud.seatable.io/",
                                 workspace_id = "57832",
                                 token_name = "BANCTABLE_TOKEN",
                                 ac = NULL){
    if(is.null(ac)) ac <- banctable_login(token_name=token_name)
    if (is.null(base_name) && is.null(table))
      stop("you must supply one of base or table name!")
    if (is.null(base_name)) {
      base = fafbseg:::flytable_base4table(table, ac = ac, cached = F)
      return(invisible(base))
    }
    if (is.null(workspace_id)) {
      wsdf = fafbseg:::flytable_workspaces(ac = ac)
      wsdf.sel = subset(wsdf, wsdf$name == base_name)
      if (nrow(wsdf.sel) == 0)
        stop("Unable to find a workspace containing basename:",
             base_name, "
Check basename and/or access permissions.")
      if (nrow(wsdf.sel) > 1)
        stop("Multiple workspaces containing basename:",
             base_name, "
You must use banctable_base() specifying a workspace_id to resolve this ambiguity.")
      workspace_id = wsdf.sel[["workspace_id"]]
    }
    base = reticulate::py_call(ac$get_base, workspace_id = workspace_id,
                               base_name = base_name)
    base
}

#' @export
#' @rdname banctable_query
banctable_move_to_bigdata <- function(table = "banc_meta",
                                      base = "banc_meta",
                                      url = "https://cloud.seatable.io/",
                                      workspace_id = "57832",
                                      token_name = "BANCTABLE_TOKEN",
                                      where = "`region` = 'optic'",
                                      invert = FALSE,
                                      row_ids = NULL){

  # get base
  ac <- banctable_login(token_name=token_name)
  base <- banctable_base_impl(table = table,
                              base_name = base,
                              url = url,
                              workspace_id = workspace_id)
  base_uuid <- base$dtable_uuid
  token <- base$jwt_token

  # Remove any protocol prefix if present
  server <- gsub("^https?://", "", base$server_url)
  server <- gsub("/$", "", server)

  # Construct the URL
  if(invert){
    movement <- "unarchive"
  }else{
    movement <- "archive-view"
  }
  endpoint <- sprintf("https://%s/api-gateway/api/v2/dtables/%s/%s/", server, base_uuid, movement)

  # Prepare the request body
  body <- list(table_name = table)
  if(invert){
    if (!is.null(row_ids)) {
      body$row_ids <- as.list(row_ids)
    }
  }else{
    # Add where clause if provided
    if (!is.null(where)) {
      body$where <- where
    }
  }

  # Make the request
  response <- httr2::request(endpoint) %>%
    httr2::req_headers(
      "Authorization" = sprintf("Bearer %s", token),
      "Accept" = "application/json",
      "Content-Type" = "application/json"
    ) %>%
    httr2::req_body_json(body) %>%
    httr2::req_error(is_error = function(resp) FALSE) %>%  # This allows us to handle errors manually
    httr2::req_perform()

  # Check for successful response
  if (httr2::resp_status(response) != 200) {
      # Try to get error message from response body
      error_msg <- tryCatch({
        if (httr2::resp_content_type(response) == "application/json") {
          error_content <- httr2::resp_body_json(response)
        } else {
          # If not JSON, get the raw text
          httr2::resp_body_string(response)
        }
      }, error = function(e) {
        "Could not parse error message"
    })
   stop(error_msg)
  }

  # Return the response
  invisible()
}

# import requests
#
# url = "https://cloud.seatable.io/api-gateway/api/v2/dtables/cc271335-227d-4bc7-94bb-1b5ae8816bd0/unarchive/"
#
# payload = {
#   "row_ids": ["FoDxhChYQSycLm88JZ11RA"],
#   "table_id": "franken_meta"
# }
# headers = {
#   "content-type": "application/json",
#   "authorization": "Bearer XX"
#
# response = requests.post(url, json=payload, headers=headers)
#
# print(response.text)

# ## in python:
# url = "https://cloud.seatable.io/api-gateway/api/v2/dtables/397da290-5aec-44dc-8a05-e2f58254d84a/archive-view/"
# headers = {
#   "accept": "application/json",
#   "content-type": "application/json",
#   "authorization": "Bearer MY_TOKEN"
# }
# body = {
#   "table_name": "banc_meta",
#   "where": "`cell_class` = 'glia'"
# }
# response = requests.post(url, headers=headers, json=body)
# print(response.text)

#' @export
#' @rdname banctable_query
franken_meta <- function(sql = "SELECT * FROM franken_meta",
                         base = "cns_meta", ...){
  # df <- banctable_query(sql=sql, base=base, ...) %>%
  #   dplyr::select(-dplyr::starts_with(c("FAFB_", "MANC_")))
  df <- banctable_query(sql=sql, base=base, ...)
  df
}

#' @export
#' @rdname banctable_query
banctable_append_rows <- function (df,
                                   table,
                                   bigdata = FALSE,
                                   base = NULL,
                                   chunksize = 1000L,
                                   workspace_id = "57832",
                                   token_name = "BANCTABLE_TOKEN",
                                   ...) {
  if (is.character(base) || is.null(base)){
    base <- banctable_base(base_name = base, table = table, workspace_id = workspace_id, token_name = token_name)
  }
  nx = nrow(df)
  if (!isTRUE(nx > 0)) {
    warning("No rows to append in `df`!")
    return(TRUE)
  }
  df = fafbseg:::df2flytable(df, append = TRUE)
  if (nx > chunksize) {
    nchunks = ceiling(nx/chunksize)
    chunkids = rep(seq_len(nchunks), rep(chunksize, nchunks))[seq_len(nx)]
    chunks = split(df, chunkids)
    if (!requireNamespace("pbapply", quietly = TRUE)) {
      stop("Package pbapply is required for this function. Please install it with: install.packages('pbapply')")
    }
    oks = pbapply::pbsapply(chunks, banctable_append_rows,
                            table = table, base = base, chunksize = Inf, bigdata = bigdata,
                            ...)
    return(all(oks))
  }
  pyl = fafbseg:::df2appendpayload(df)
  if(!bigdata){
    res = base$batch_append_rows(table_name = table, rows_data = pyl)
  }else{
    res = base$big_data_insert_rows(table_name = table, rows_data = pyl)
  }
  ok = isTRUE(all.equal(res[["inserted_row_count"]], nx))
  return(ok)
}

# modified to enable list uploads to multi-select columns
banc_df2updatepayload <- function(x, via_json = TRUE){
  if (via_json) {
    othercols <- setdiff(colnames(x), "row_id")
    listcols <- names(x)[sapply(x, is.list)]
    listcols <- intersect(othercols, listcols)
    updates <- list()
    for(i in 1:nrow(x)){
      updates[[i]] <- list(row_id = x[i, "row_id"], row = as.list(x[i,othercols]))
      for(col in listcols){
        if(length((x[i,][[col]][[1]]))==1){
          updates[[i]]$row[[col]] <- x[i,][[col]]
        }else{
          updates[[i]]$row[[col]] <- x[i,][[col]][[1]]
        }
      }
    }
    js <- jsonlite::toJSON(updates, auto_unbox = TRUE)
    pyjson <- reticulate::import("json")
    pyl <- reticulate::py_call(pyjson$loads, js)
    return(pyl)
  }
  pdf = reticulate::r_to_py(x)
  pyfun = fafbseg:::df2updatepayload_py()
  reticulate::py_call(pyfun$pdf2list, pdf)
}

# hidden, modified to enable working with list columns
banctable2df <- function (df, tidf = NULL) {
  if (!isTRUE(ncol(df) > 0))
    return(df)
  nr = nrow(df)
  listcols = sapply(df, is.list)
  for (i in which(listcols)) {
    li = lengths(df[[i]])
    if (isTRUE(all(li == 1))) {
      ul = unlist(df[[i]])
      if (!isTRUE(length(ul) == nr))
        ul = sapply(ul,paste,collapse=",")
      else df[[i]] = ul
    }
    else if (isTRUE(all(li %in% 0:1))) {
      df[[i]][!nzchar(df[[i]])] = NA
      df[[i]] = fafbseg:::null2na(df[[i]])
    }
    else df[[i]] = sapply(df[[i]],paste,collapse=",")
  }
  if (is.null(tidf))
    df
  else {
    if (is.character(tidf))
      tidf = fafbseg::flytable_columns(tidf)
    fafbseg:::flytable_fix_coltypes(df, tidf = tidf)
  }
}

# hidden, helper function to update status column
banc_update_status <- function(df,
                               update,
                               col = "status",
                               wipe = FALSE){
  if(wipe){
    df$status <- ""
  }else{
    df$status[is.na(df$status)] <- ""
    df$status[df$status%in%c("NA","NaN")] <- ""
  }
  update.col <- sapply(df$status, function(x){
    x=paste(c(x,update),collapse=",")
    paste(sort(unique(unlist(strsplit(x,split=",|, ")))),collapse=",")
  }
  )
  update.col <- gsub("^,","",update.col)
  df[[col]] <- update.col
  df
}

# # # Example of adding a labels to the status column
# bc <- banctable_query()
# sizes <- as.numeric(bc$l2_cable_length_um)
# sizes[is.na(sizes)] <- 0
# tadpoles <- bc[sizes>1&sizes<10,]
# tadpoles <- banc_update_status(tadpoles,update="TOO_SMALL")
# banctable_update_rows(base = 'banc_meta',
#                       table = "banc_meta",
#                       df = tadpoles[,c("_id","super_class","status")],
#                       append_allowed = FALSE,
#                       chunksize = 100)

# Update the BANC IDs
banctable_updateids <- function(){

  # Get cell info table
  cat('reading cell info cave table...
')
  info <- banc_cell_info(rawcoords = TRUE)  %>%
    dplyr::mutate(pt_position = xyzmatrix2str(.data$pt_position)) %>%
    dplyr::select(.data$pt_root_id, .data$pt_supervoxel_id, .data$pt_position) %>%
    rbind(banc_backbone_proofread() %>%
            dplyr::select(.data$pt_root_id, .data$pt_supervoxel_id, .data$pt_position) %>%
            dplyr::mutate(pt_position = xyzmatrix2str(.data$pt_position))) %>%
    rbind(banc_neck_connective_neurons() %>%
            dplyr::select(.data$pt_root_id, .data$pt_supervoxel_id, .data$pt_position) %>%
            dplyr::mutate(pt_position = xyzmatrix2str(.data$pt_position))) %>%
    dplyr::mutate(pt_root_id=as.character(.data$pt_root_id),
                  pt_supervoxel_id=as.character(.data$pt_supervoxel_id)) %>%
    dplyr::distinct(.data$pt_supervoxel_id, .keep_all = TRUE) %>%
    dplyr::rowwise()

  # Get current table
  cat('reading banc meta seatable...
')
  bc <- banctable_query(sql = 'select _id, root_id, supervoxel_id, position, banc_match, banc_match_supervoxel_id, banc_png_match, banc_png_match_supervoxel_id, banc_nblast_match, banc_nblast_match_supervoxel_id from banc_meta') %>%
    dplyr::select(.data$root_id, .data$supervoxel_id, .data$position,
                  .data$banc_match, .data$banc_match_supervoxel_id, .data$banc_png_match, .data$banc_png_match_supervoxel_id, .data$banc_nblast_match, .data$banc_nblast_match_supervoxel_id,
                  .data$`_id`)
  bc[bc=="0"] <- NA
  bc[bc==""] <- NA

  # Update
  cat('updating column: root_id ...
')
  bc.new <- bc %>%
    dplyr::left_join(info,
                     by = c("supervoxel_id"="pt_supervoxel_id")) %>%
    dplyr::mutate(root_id = ifelse(is.na(.data$pt_root_id), .data$root_id, .data$pt_root_id)) %>%
    dplyr::mutate(position = ifelse(is.na(.data$position), .data$pt_position, .data$position)) %>%
    dplyr::select(-.data$pt_root_id, -.data$pt_position)

  # Update root IDs directly where needed
  bc.new <- banc_updateids(bc.new,
                           root.column = "root_id",
                           supervoxel.column = "supervoxel_id",
                           position.column = "position")

  # Make sure supervoxel and root position information that is missing, is filled in
  bc.new <- bc.new %>%
    dplyr::left_join(info %>% dplyr::distinct(.data$pt_root_id, .keep_all = TRUE),
                     by = c("root_id"="pt_root_id")) %>%
    dplyr::mutate(supervoxel_id = ifelse(is.na(.data$supervoxel_id), .data$pt_supervoxel_id, .data$supervoxel_id)) %>%
    dplyr::mutate(position = ifelse(is.na(.data$position), .data$pt_position, .data$position)) %>%
    dplyr::select(-.data$pt_supervoxel_id, -.data$pt_position)

  # Update match columns
  lookup <- bc.new %>%
    dplyr::select(.data$root_id, .data$supervoxel_id) %>%
    dplyr::rename(lookup_root_id=.data$root_id,
                  lookup_supervoxel_id=.data$supervoxel_id) %>%
    dplyr::filter(!is.na(.data$lookup_supervoxel_id), .data$lookup_supervoxel_id!="0",
                  !is.na(.data$lookup_root_id), .data$lookup_root_id!="0") %>%
    dplyr::distinct(.data$lookup_root_id, .data$lookup_supervoxel_id)
  bc.new <- bc.new %>%
    dplyr::left_join(lookup, by = c("banc_match_supervoxel_id"="lookup_supervoxel_id")) %>%
    dplyr::mutate(banc_match = dplyr::case_when(
      !is.na(.data$lookup_root_id) ~ .data$lookup_root_id,
      TRUE ~ .data$banc_match
    )) %>%
    dplyr::select(-.data$lookup_root_id) %>%
    dplyr::left_join(lookup, by = c("banc_png_match_supervoxel_id"="lookup_supervoxel_id")) %>%
    dplyr::mutate(banc_png_match = dplyr::case_when(
      !is.na(.data$lookup_root_id) ~ .data$lookup_root_id,
      TRUE ~ .data$banc_png_match
    )) %>%
    dplyr::select(-.data$lookup_root_id) %>%
    dplyr::left_join(lookup, by = c("banc_nblast_match_supervoxel_id"="lookup_supervoxel_id")) %>%
    dplyr::mutate(banc_nblast_match = dplyr::case_when(
      !is.na(.data$lookup_root_id) ~ .data$lookup_root_id,
      TRUE ~ .data$banc_nblast_match
    )) %>%
    dplyr::select(-.data$lookup_root_id)

  # Update directly
  cat('updating column: banc_match ...
')
  bc.new[!is.na(bc.new$banc_match),] <- banc_updateids(bc.new[!is.na(bc.new$banc_match),],
                                                       root.column = "banc_match",
                                                       supervoxel.column = "banc_match_supervoxel_id",
                                                       position.column = "banc_match_position")
  cat('updating column: banc_png_match ...
')
  bc.new[!is.na(bc.new$banc_png_match),] <- banc_updateids(bc.new[!is.na(bc.new$banc_png_match),],
                                                       root.column = "banc_png_match",
                                                       supervoxel.column = "banc_png_match_supervoxel_id",
                                                       position.column = "banc_png_match_position")
  cat('updating column: banc_nblast_match ...
')
  bc.new[!is.na(bc.new$banc_nblast_match),] <- banc_updateids(bc.new[!is.na(bc.new$banc_nblast_match),],
                                                       root.column = "banc_nblast_match",
                                                       supervoxel.column = "banc_nblast_match_supervoxel_id",
                                                       position.column = "banc_nblast_match_position")
  bc.new <- bc.new %>%
    dplyr::left_join(lookup %>%dplyr::distinct(.data$lookup_root_id, .keep_all=TRUE),
                     by = c("banc_match"="lookup_root_id")) %>%
    dplyr::mutate(banc_match_supervoxel_id = dplyr::case_when(
      is.na(.data$banc_match_supervoxel_id)&!is.na(.data$lookup_supervoxel_id) ~ .data$lookup_supervoxel_id,
      TRUE ~ .data$banc_match_supervoxel_id
    )) %>%
    dplyr::select(-.data$lookup_supervoxel_id) %>%
    dplyr::left_join(lookup %>%dplyr::distinct(.data$lookup_root_id, .keep_all=TRUE),
                     by = c("banc_png_match"="lookup_root_id")) %>%
    dplyr::mutate(banc_png_match_supervoxel_id = dplyr::case_when(
      is.na(.data$banc_png_match_supervoxel_id)&!is.na(.data$lookup_supervoxel_id) ~ .data$lookup_supervoxel_id,
      TRUE ~ .data$banc_png_match_supervoxel_id
    )) %>%
    dplyr::select(-.data$lookup_supervoxel_id) %>%
    dplyr::left_join(lookup %>%dplyr::distinct(.data$lookup_root_id, .keep_all=TRUE),
                     by = c("banc_nblast_match"="lookup_root_id")) %>%
    dplyr::mutate(banc_nblast_match_supervoxel_id = dplyr::case_when(
      is.na(.data$banc_nblast_match_supervoxel_id)&!is.na(.data$lookup_supervoxel_id) ~ .data$lookup_supervoxel_id,
      TRUE ~ .data$banc_nblast_match_supervoxel_id
    )) %>%
    dplyr::select(-.data$lookup_supervoxel_id)

  # Update
  cat('updating banc meta seatable...
')
  bc.new <- bc.new %>%
    dplyr::filter(!is.na(.data$`_id`)) %>%
    dplyr::distinct(.data$`_id`, .keep_all = TRUE)
  bc.new[is.na(bc.new)] <- ''
  bc.new[bc.new=="0"] <- ''
  banctable_update_rows(df = bc.new,
                        base = "banc_meta",
                        table = "banc_meta",
                        append_allowed = FALSE,
                        chunksize = 1000)
  cat('done.')

  # Return
  invisible()
}

banctable_annotate <- function(root_ids,
                               update,
                               overwrite = FALSE,
                               append = FALSE,
                               column="notes"){


  # Get current table
  cat('reading banc meta seatable...
')
  bc <- banctable_query(sql = sprintf('select _id, root_id, supervoxel_id, %s from banc_meta',column)) %>%
    dplyr::filter(.data$root_id %in% root_ids)
  if(!nrow(bc)){
    message("root_ids not in BANC meta")
    return(invisible())
  }
  bc[bc=="0"] <- NA
  bc[bc==""] <- NA

  # Update
  cat('updating column: root_id ...
')
  bc.new <- bc
  if(overwrite){
    bc.new[[column]] <- NA
  }
  if(append){
    bc.new <- bc.new %>%
      dplyr::rowwise() %>%
      dplyr::mutate(update = dplyr::case_when(
        is.na(.data[[column]]) ~ update,
        TRUE ~ paste(.data[[column]], update, sep = ", ", collapse = ", "),
      )
      )
  }else{
    bc.new <- bc.new %>%
      dplyr::rowwise() %>%
      dplyr::mutate(update = dplyr::case_when(
        is.na(.data[[column]]) ~ update,
        TRUE ~ .data[[column]],
      )
    )
  }
  changed <- sum(bc.new[[column]] != bc.new$update, na.rm = TRUE)
  bc.new[[column]] <- bc.new$update
  bc.new$update <- NULL

  # Summarise update
  message("Changed ", changed, " rows")
  cat(sprintf("%s before update:
",column))
  if(nrow(bc.new)==1){
    cat(bc[[column]])
  }else{
    cat(sort(table(bc[[column]])))
  }
  cat(sprintf("
 %s after update:
",column))
  if(nrow(bc.new)==1){
    cat(bc.new[[column]])
  }else{
    cat(sort(table(bc.new[[column]])))
  }

  # Update
  cat('updating banc meta seatable...
')
  bc.new <- as.data.frame(bc.new)
  bc.new[is.na(bc.new)] <- ''
  bc.new[bc.new=="0"] <- ''
  banctable_update_rows(df = bc.new,
                        base = "banc_meta",
                        table = "banc_meta",
                        append_allowed = FALSE,
                        chunksize = 1000)
  cat('done.')

  # Return
  invisible()

}


# hidden
banctable_snapshots <- function(token_name = "BANCTABLE_TOKEN",
                                base_name = "banc_meta",
                                workspace_id = "57832"){
  # Build the URL
  url <- sprintf(
    "https://cloud.seatable.io/api/v2.1/workspace/%s/dtable/%s/snapshots/",
    workspace_id, base_name
  )

  # Set headers
  token <- Sys.getenv(token_name, unset = NA_character_)
  headers <- add_headers(
    accept = "application/json",
    authorization = paste("Bearer", token)
  )

  # Perform the GET request
  response <- httr::GET(url, headers)

  # Given content_text:
  content_text <- content(response, as = "text", encoding = "UTF-8")
  json_result <- jsonlite::fromJSON(content_text)

  # To get just the snapshot list as a data.frame:
  snapshot_list <- json_result$snapshot_list

  # See the first few rows:
  snapshot_list
}




