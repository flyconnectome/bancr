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
# @param workspace_id A numeric id specifying the workspace. Advanced use only
#   since we can normally figure this out from \code{base_name}.
# @param cached Whether to use a cached base object
#' @param token normally retrieved from \code{BANCTABLE_TOKEN} environment
#'   variable.
#' @param user,pwd banctable user and password used by \code{banctable_set_token}
#'   to obtain a token
#' @param url Optional URL to the server
#' @param ac A seatable connection object as returned by \code{banctable_login}.
#' @param df A data.frame containing the data to upload including an `_id`
#' column that can identify each row in the remote table.
#' @param append_allowed Logical. Whether rows without row identifiers can be appended.
#' @param chunksize To split large requests into smaller ones with max this many rows.
#' @param ... Additional arguments passed to pbsapply which might include cl=2 to specify a number of parallel jobs to run.
#'
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
                             ac = NULL){
  if(is.null(ac)) ac <- banctable_login()
  seatable.max <- 10000L
  if(limit>seatable.max){
    offset <- 0
    df <- data.frame()
    while(offset<limit){
      cat("reading from row: ", offset,"\n")
      sql.new <- sprintf("%s LIMIT %d OFFSET %d", sql, seatable.max, offset)
      bc <- banctable_query(sql=sql.new,
                            limit=FALSE,
                            base=base,
                            python=python,
                            convert=convert,
                            ac=ac)
      df <- rbind(df,bc)
      offset <- offset+nrow(bc)
      if(!length(bc)|nrow(bc)<seatable.max){
        cat("read rows: ",nrow(df), " read columns:", ncol(df),"\n")
        return(df)
      }
    }
    cat("read rows: ",nrow(df), " read columns:", ncol(df),"\n")
    return(df)
  }
  checkmate::assert_character(sql, len = 1, pattern = "select",
                              ignore.case = T)
  res = stringr::str_match(sql, stringr::regex("\\s+FROM\\s+[']{0,1}([^, ']+).*",
                                               ignore_case = T))
  if (any(is.na(res)[, 2]))
    stop("Cannot identify a table name in your sql statement!\n")
  table = res[, 2]
  if (is.null(base)) {
    base = try(banctable_base(table = table))
    if (inherits(base, "try-error"))
      stop("I inferred table_name: ", table, " from your SQL query but couldn't connect to a base with this table!")
  }
  else if (is.character(base))
    base = banctable_base(base_name = base)
  if (!isTRUE(grepl("\\s+limit\\s+\\d+", sql)) && !isFALSE(limit)) {
    if (!is.finite(limit))
      limit = .Machine$integer.max
    sql = paste(sql, "LIMIT", limit)
  }
  pyout <- reticulate::py_capture_output(ll <- try(reticulate::py_call(base$query,
                                                                       sql, convert = convert), silent = T))
  if (inherits(ll, "try-error")) {
    warning(paste("No rows returned by banctable", pyout,
                  collapse = "\n"))
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
banctable_set_token <- function(user, pwd, url = "https://cloud.seatable.io/"){
  st <- fafbseg:::check_seatable()
  ac <- reticulate::py_call(st$Account, login_name = user,
                            password = pwd, server_url = url)
  ac$auth()
  Sys.setenv(banctable_TOKEN = ac$token)
  cat("BANCTABLE_TOKEN='", ac$token, "'\n", sep = "", append = TRUE,
      file = path.expand("~/.Renviron"))
  return(invisible(NULL))
}

#' @export
#' @rdname banctable_query
banctable_login <- function(url = "https://cloud.seatable.io/",
                            token = Sys.getenv("BANCTABLE_TOKEN", unset = NA_character_)){
  fafbseg::flytable_login(url=url, token=token)
}


#' @export
#' @rdname banctable_query
banctable_update_rows <- function (df, table, base = NULL, append_allowed = FALSE, chunksize = 1000L,  ...) {
  df <- as.data.frame(df)
  if (is.character(base) || is.null(base))
    base = banctable_base(base_name = base, table = table)
  nx = nrow(df)
  if (!isTRUE(nx > 0)) {
    warning("No rows to update in `df`!")
    return(TRUE)
  }
  tablecols = fafbseg:::flytable_columns(table,base)
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
banctable_base <- function (base_name = "banc_meta",
                            table = NULL,
                            url = "https://cloud.seatable.io/",
                            workspace_id = "57832",
                            cached = TRUE,
                            ac = NULL) {
  if(is.null(ac)) ac <- banctable_login()
  if (!cached)
    memoise::forget(banctable_base_impl)
  base = try({
    banctable_base_impl(table = table, base_name = base_name,
                       url = url, workspace_id = workspace_id)
  }, silent = TRUE)
  stale_token <- isTRUE(try(difftime(base$jwt_exp, Sys.time(),
                                     units = "hours") < 1, silent = T))
  retry = (cached && inherits(base, "try-error")) || stale_token
  if (!retry)
    return(base)
  memoise::forget(banctable_base_impl)
  banctable_base_impl(table = table, base_name = base_name,
                     url = url, workspace_id = workspace_id)
}

# hidden
banctable_base_impl <- function (base_name = "banc_meta",
                                 table = NULL,
                                 url = "https://cloud.seatable.io/",
                                 workspace_id = "57832",
                                 ac = NULL){
    if(is.null(ac)) ac <- banctable_login()
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
             base_name, "\nCheck basename and/or access permissions.")
      if (nrow(wsdf.sel) > 1)
        stop("Multiple workspaces containing basename:",
             base_name, "\nYou must use banctable_base() specifying a workspace_id to resolve this ambiguity.")
      workspace_id = wsdf.sel[["workspace_id"]]
    }
    base = reticulate::py_call(ac$get_base, workspace_id = workspace_id,
                               base_name = base_name)
    base
}

#' @export
#' @rdname banctable_query
banctable_append_rows <- function (df, table, base = NULL, chunksize = 1000L, ...) {
  if (is.character(base) || is.null(base)){
    base = banctable_base(base_name = base, table = table)
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
    oks = pbapply::pbsapply(chunks, banctable_append_rows,
                            table = table, base = base, chunksize = Inf, ...)
    return(all(oks))
  }
  pyl = fafbseg:::df2appendpayload(df)
  res = base$batch_append_rows(table_name = table, rows_data = pyl)
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
      tidf = fafbseg:::flytable_columns(tidf)
    fafbseg:::flytable_fix_coltypes(df, tidf = tidf)
  }
}

# hidden, helper function to update status column
banc_update_status <- function(df, update, col = "status", wipe = FALSE){
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
  cat('reading cell info cave table...\n')
  info <- banc_cell_info(rawcoords = TRUE)  %>%
    dplyr::mutate(pt_position = xyzmatrix2str(pt_position)) %>%
    dplyr::select(pt_root_id, pt_supervoxel_id,pt_position) %>%
    rbind(banc_backbone_proofread() %>%
            dplyr::select(pt_root_id, pt_supervoxel_id,pt_position) %>%
            dplyr::mutate(pt_position = xyzmatrix2str(pt_position))) %>%
    rbind(banc_neck_connective_neurons() %>%
            dplyr::select(pt_root_id, pt_supervoxel_id,pt_position) %>%
            dplyr::mutate(pt_position = xyzmatrix2str(pt_position))) %>%
    dplyr::mutate(pt_root_id=as.character(pt_root_id),
                  pt_supervoxel_id=as.character(pt_supervoxel_id)) %>%
    dplyr::distinct(pt_supervoxel_id, .keep_all = TRUE) %>%
    dplyr::rowwise()

  # Get current table
  cat('reading banc meta seatable...\n')
  bc <- banctable_query(sql = 'select _id, root_id, supervoxel_id, position, banc_match, banc_match_supervoxel_id, banc_png_match, banc_png_match_supervoxel_id, banc_nblast_match, banc_nblast_match_supervoxel_id from banc_meta') %>%
    dplyr::select(root_id, supervoxel_id, position,
                  banc_match, banc_match_supervoxel_id, banc_png_match, banc_png_match_supervoxel_id, banc_nblast_match, banc_nblast_match_supervoxel_id,
                  `_id`)
  bc[bc=="0"] <- NA
  bc[bc==""] <- NA

  # Update
  cat('updating column: root_id ...\n')
  bc.new <- bc %>%
    dplyr::left_join(info,
                     by = c("supervoxel_id"="pt_supervoxel_id")) %>%
    dplyr::mutate(root_id = ifelse(is.na(pt_root_id),root_id,pt_root_id)) %>%
    dplyr::mutate(position = ifelse(is.na(position),pt_position,position)) %>%
    dplyr::select(-pt_root_id,-pt_position)

  # Update root IDs directly where needed
  bc.new <- banc_updateids(bc.new, root.column = "root_id", supervoxel.column = "supervoxel_id")

  # Make sure supervoxel and root position information that is missing, is filled in
  bc.new <- bc.new %>%
    dplyr::left_join(info %>% dplyr::distinct(pt_root_id, .keep_all = TRUE),
                     by = c("root_id"="pt_root_id")) %>%
    dplyr::mutate(supervoxel_id = ifelse(is.na(supervoxel_id),pt_supervoxel_id,supervoxel_id)) %>%
    dplyr::mutate(position = ifelse(is.na(position),pt_position,position)) %>%
    dplyr::select(-pt_supervoxel_id,-pt_position)

  # Update match columns
  lookup <- bc.new %>%
    dplyr::select(root_id, supervoxel_id) %>%
    dplyr::rename(lookup_root_id=root_id,
                  lookup_supervoxel_id=supervoxel_id) %>%
    dplyr::filter(!is.na(lookup_supervoxel_id), lookup_supervoxel_id!="0",
                  !is.na(lookup_root_id), lookup_root_id!="0") %>%
    dplyr::distinct(lookup_root_id,lookup_supervoxel_id)
  bc.new <- bc.new %>%
    dplyr::left_join(lookup, by = c("banc_match_supervoxel_id"="lookup_supervoxel_id")) %>%
    dplyr::mutate(banc_match = dplyr::case_when(
      !is.na(lookup_root_id) ~ lookup_root_id,
      TRUE ~ banc_match
    )) %>%
    dplyr::select(-lookup_root_id) %>%
    dplyr::left_join(lookup, by = c("banc_png_match_supervoxel_id"="lookup_supervoxel_id")) %>%
    dplyr::mutate(banc_png_match = dplyr::case_when(
      !is.na(lookup_root_id) ~ lookup_root_id,
      TRUE ~ banc_png_match
    )) %>%
    dplyr::select(-lookup_root_id) %>%
    dplyr::left_join(lookup, by = c("banc_nblast_match_supervoxel_id"="lookup_supervoxel_id")) %>%
    dplyr::mutate(banc_nblast_match = dplyr::case_when(
      !is.na(lookup_root_id) ~ lookup_root_id,
      TRUE ~ banc_nblast_match
    )) %>%
    dplyr::select(-lookup_root_id)

  # Update directly
  cat('updating column: banc_match ...\n')
  bc.new <- banc_updateids(bc.new, root.column = "banc_match", supervoxel.column = "banc_match_supervoxel_id")
  cat('updating column: banc_png_match ...\n')
  bc.new <- banc_updateids(bc.new, root.column = "banc_png_match", supervoxel.column = "banc_png_match_supervoxel_id")
  cat('updating column: banc_nblast_match ...\n')
  bc.new <- banc_updateids(bc.new, root.column = "banc_nblast_match", supervoxel.column = "banc_nblast_match_supervoxel_id")
  bc.new <- bc.new %>%
    dplyr::left_join(lookup %>%dplyr::distinct(lookup_root_id, .keep_all=TRUE),
                     by = c("banc_match"="lookup_root_id")) %>%
    dplyr::mutate(banc_match_supervoxel_id = dplyr::case_when(
      is.na(banc_match_supervoxel_id)&!is.na(lookup_supervoxel_id) ~ lookup_supervoxel_id,
      TRUE ~ banc_match_supervoxel_id
    )) %>%
    dplyr::select(-lookup_supervoxel_id) %>%
    dplyr::left_join(lookup %>%dplyr::distinct(lookup_root_id, .keep_all=TRUE),
                     by = c("banc_png_match"="lookup_root_id")) %>%
    dplyr::mutate(banc_png_match_supervoxel_id = dplyr::case_when(
      is.na(banc_png_match_supervoxel_id)&!is.na(lookup_supervoxel_id) ~ lookup_supervoxel_id,
      TRUE ~ banc_png_match_supervoxel_id
    )) %>%
    dplyr::select(-lookup_supervoxel_id) %>%
    dplyr::left_join(lookup %>%dplyr::distinct(lookup_root_id, .keep_all=TRUE),
                     by = c("banc_nblast_match"="lookup_root_id")) %>%
    dplyr::mutate(banc_nblast_match_supervoxel_id = dplyr::case_when(
      is.na(banc_nblast_match_supervoxel_id)&!is.na(lookup_supervoxel_id) ~ lookup_supervoxel_id,
      TRUE ~ banc_nblast_match_supervoxel_id
    )) %>%
    dplyr::select(-lookup_supervoxel_id)

  # Update
  cat('updating banc meta seatable...\n')
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




