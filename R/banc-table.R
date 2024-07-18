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
#' @param sql A SQL query string. See examples and
#'   \href{https://seatable.github.io/seatable-scripts/python/query/}{seatable
#'   docs}.
#' @param limit An optional limit, which only applies if you do not specify a
#'   limit directly in the \code{sql} query. By default seatable limits SQL
#'   queries to 100 rows. We increase the limit to 100000 rows by default.
#' @param convert Expert use only: Whether or not to allow the Python seatable
#'   module to process raw output from the database. This is is principally for
#'   debugging purposes. NB this imposes a requirement of seatable_api >=2.4.0.
#' @param base_name Character vector specifying the \code{base}
#' @param table Character vector specifying a table foe which you want a
#'   \code{base} object.
#' @param workspace_id A numeric id specifying the workspace. Advanced use only
#'   since we can normally figure this out from \code{base_name}.
#' @param cached Whether to use a cached base object
#' @param token normally retrieved from \code{BANCTABLE_TOKEN} environment
#'   variable.
#' @param user,pwd banctable user and password used by \code{banctable_set_token}
#'   to obtain a token
#' @param url Optional URL to the server
#' @param ac A seatable connection object as returned by \code{banc_table_login}.
#'
#' @return a \code{data.frame} of results. There should be 0 rows if no rows
#'   matched query.
#'
#' @seealso \code{\link{fafbseg::flytable_query}}
#' @examples
#' \donttest{
#' banc_set_token(user="alexander.shakeel.bates@gmail.com",
#               pwd="MY_PASSWORD",
#               url="https://cloud.seatable.io/")
#' banc.meta <- banctable_query()
#' }
#' @export
#' @rdname banctable_query
banctable_query <- function (sql = "SELECT * FROM banc_meta",
                             limit = 100000L,
                             base = NULL,
                             python = FALSE,
                             convert = TRUE,
                             ac = banc_table_login()){
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
    colinfo = fafbseg:::flytable_columns(table, base)
    df = fafbseg:::flytable2df(fafbseg:::pandas2df(pdd, use_arrow = F), tidf = colinfo)
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
banc_set_token <- function (user, pwd, url = "https://cloud.seatable.io/"){
  st <- check_seatable()
  ac <- reticulate::py_call(st$Account, login_name = user,
                            password = pwd, server_url = url)
  ac$auth()
  Sys.setenv(banctable_TOKEN = ac$token)
  cat("banctable_TOKEN='", ac$token, "'\n", sep = "", append = TRUE,
      file = path.expand("~/.Renviron"))
  return(invisible(NULL))
}

#' @export
#' @rdname banctable_query
banc_table_login <- function(url = "https://cloud.seatable.io/",
                             token = Sys.getenv("BANCTABLE_TOKEN", unset = NA_character_)){
  fafbseg::banctable_login(url=url, token=token)
}

# hidden
banctable_base <- function (base_name = "banc_meta",
                            table = NULL,
                            url = "https://cloud.seatable.io/",
                            workspace_id = "57832",
                            cached = TRUE,
                            ac = banc_table_login()) {
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
                                 ac = banc_table_login()){
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
