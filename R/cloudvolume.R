# hidden
banc_cloudvolume <- function(...) {
  cv<-fafbseg::flywire_cloudvolume(cloudvolume.url = banc_cloudvolume_url(), ...)
  cv
}

# hidden
banc_cloudvolume_url <- function() {
  rr=with_banc(getOption("fafbseg.cloudvolume.url"))
  sub("graphene://middleauth+", "graphene://", rr, fixed = TRUE)
}

# hidden
banc_api_url <- function(endpoint="") {
  fafbseg:::flywire_api_url(endpoint=endpoint,
                            cloudvolume.url = banc_cloudvolume_url())
}

#' Set the token to be used to authenticate to banc autosegmentation resources
#'
#' @param token An optional token string. When missing you are prompted to
#'   generate a new token via your browser.
#'
#' @return The path to the token file (invisibly)
#' @export
banc_set_token <- function(token=NULL) {
  # path=fafbseg::flywire_set_token(token=token, domain='cave.banc-fly.com')
  path=fafbseg::flywire_set_token(token=token, domain = NULL)
  # clear the token cache so the new one is immediately available
  banc_token(cached=FALSE)
  invisible(path)
}

# hidden
banc_token <- function(cached=TRUE) {
  # fafbseg::chunkedgraph_token(url='cave.banc-fly.com', cached = cached)
  fafbseg::chunkedgraph_token(cached = cached)
}

# hidden
banc_token_available <- function() {
  !inherits(try(banc_token(), silent = TRUE), 'try-error')
}

#' Print information about your BANC setup including tokens and python modules
#'
#' @seealso \code{\link{dr_fafbseg}}
#' @examples
#' \dontrun{
#' dr_banc()
#' }
#' @export
dr_banc <- function() {
  banc_api_report()
  cat("

")
  res = fafbseg:::py_report()
  cat("
")
  try(fafbseg:::check_cloudvolume_reticulate(min_version = "8.32.1"))
  invisible(res)
}

# hidden
banc_api_report <- function(with_dataset = bancr::with_banc) {
  message("Neuroglancer / CAVE API access
----")
  token=try(banc_token(cached = F), silent = FALSE)
  if(inherits(token, "try-error")) {
    FUN=if(requireNamespace('usethis', quietly = T)) usethis::ui_todo else message
    FUN(paste('No valid banc API token found. Set your token by doing:
',
                  "{ui_code('*_set_token()')}"))
  } else{
    cat("Valid banc API ChunkedGraph token is set!
")
  }
  ff=dir(fafbseg:::cv_secretdir(), pattern = '-secret\\.json$')
  if(length(ff)){
    cat(length(ff), "CloudVolume credential files available at
",
        fafbseg:::cv_secretdir(),"
")
    print(ff)
  }
  u=with_dataset(fafbseg:::check_cloudvolume_url(set = F))
  cat("
Zetta cloudvolume URL:", u)
}

