#' Return a sample Neuroglancer scene URL for banc dataset
#'
#' @details See
#'   \href{https://banc-reconstruction.slack.com/archives/C01RZP5JH9C/p1616522511001900}{banc
#'    slack} for details.
#'
#' @param ids A set of root ids to include in the scene. Can also be a
#'   data.frame.
#' @param open Whether to open the URL in your browser (see
#'   \code{\link{browseURL}})
#' @return A character vector containing a single Neuroglancer URL (invisibly
#'   when \code{open=TRUE}).
#' @export
#' @importFrom utils browseURL
#' @examples
#' \dontrun{
#' browseURL(banc_scene())
#' banc_scene(open=T)
#' banc_scene("648518346498254576", open=T)
#' }
banc_scene <- function(ids=NULL, open=FALSE) {
  url="https://spelunker.cave-explorer.org/#!middleauth+https://global.daf-apis.com/nglstate/api/v1/6283844278812672"
  url=sub("#!middleauth+", "?", url, fixed = T)
  parts=unlist(strsplit(url, "?", fixed = T))
  json=try(fafbseg::flywire_fetch(parts[2], token=banc_token(), return = 'text', cache = TRUE))
  if(inherits(json, 'try-error')) {
    badtoken=paste0("You have a token but it doesn't seem to be authorised for banc.\n",
                    "Have you definitely used `banc_set_token()` to make a token for the banc dataset?")
    if(grepl(500, json))
      stop("There seems to be a (temporary?) problem with the zetta server!")
    else if(grepl(401, json))
      stop(badtoken)

    token=try(banc_token(), silent = T)
    if(inherits(token, 'try-error'))
      stop("It looks like you do not have a stored token. Please use `banc_set_token()` to make one.")
    else
      stop(badtoken)
  }
  u=ngl_encode_url(json, baseurl = parts[1])
  if(!is.null(ids)){
    fafbseg::ngl_segments(u) <- banc_ids(ids)
  }
  if(open) {
    browseURL(u)
    invisible(u)
  } else (u)
}

#' Choose or (temporarily) use the banc autosegmentation
#'
#' @details \code{bancr} inherits a significant amount of infrastructure from
#'   the \code{\link{fafbseg}} package. This has the concept of the
#'   \emph{active} autosegmentation, which in turn defines one or more R options
#'   containing URIs pointing to voxel-wise segmentation, mesh etc data. These
#'   are normally contained within a single neuroglancer URL which points to
#'   multiple data layers. For banc this is the neuroglancer scene returned by
#'   \code{\link{banc_scene}}.
#' @param set Whether or not to permanently set the banc autosegmentation as the
#'   default for \code{\link{fafbseg}} functions.

#'
#' @return If \code{set=TRUE} a list containing the previous values of the
#'   relevant global options (in the style of \code{\link{options}}. If
#'   \code{set=FALSE} a named list containing the option values.
#' @export
#'
#' @examples
#' \dontrun{
#' choose_banc()
#' options()[grep("^fafbseg.*url", names(options()))]
#' }
choose_banc <- function(set=TRUE) {
  fafbseg::choose_segmentation(
    banc_scene(),
    set=set,
    moreoptions=list(
      fafbseg.cave.datastack_name=banc_datastack_name()
      ))
}

#' @param expr An expression to evaluate while banc is the default
#'   autosegmentation
#' @rdname choose_banc
#' @export
#' @examples
#' \donttest{
#' with_banc(fafbseg::flywire_islatest('648518346498254576'))
#' }
#' \dontrun{
#' with_banc(fafbseg::flywire_latestid('648518346498254576'))
#' with_banc(fafbseg::flywire_latestid('648518346494405175'))
#' }
with_banc <- function(expr) {
  op <- choose_banc(set = TRUE)
  on.exit(options(op))
  force(expr)
}

banc_fetch <- function(url, token=banc_token(), ...) {
  flywire_fetch(url, token=token, ...)
}


