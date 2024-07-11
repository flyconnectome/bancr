#' Find the root identifier of a banc neuron
#'
#' @inheritParams fafbseg::flywire_rootid
#' @param ... Additional arguments passed to \code{pbapply::pbsapply} and
#'   eventually to Python \code{cv$CloudVolume} object.
#' @return A vector of root ids (by default character)
#' @export
#' @family banc-ids
#' @seealso \code{\link{flywire_rootid}}
#' @examples
#' \donttest{
#' banc_rootid("73186243730767724")
#' }
banc_rootid <- function(x, integer64 = FALSE, ...) {
  rid = flywire_rootid(
    x = x,
    integer64 = integer64,
    method = "cloudvolume",
    # agglomerate = T,
    cloudvolume.url = banc_cloudvolume_url(),
    ...
  )
  rid
}

#' Find the supervoxel identifiers of a banc neuron
#'
#' @param ... additional arguments passed to \code{\link{flywire_leaves}}
#' @inheritParams fafbseg::flywire_leaves
#'
#' @return A vector of supervoxel ids
#' @family banc-ids
#' @seealso \code{\link{flywire_leaves}}
#' @export
#'
#' @examples
#' \dontrun{
#' svids=banc_leaves("720575941650432785")
#' head(svids)
#' }
banc_leaves <- function(x, integer64=TRUE, ...) {
  svids=flywire_leaves(x=x, integer64 = integer64, cloudvolume.url = banc_cloudvolume_url(), ...)
  svids
}

#' Convert xyz locations to root or supervoxel ids
#'
#' @details This used to be very slow because we do not have a supervoxel
#'   field on spine.
#'
#'   I am somewhat puzzled by the voxel dimensions for banc. Neuroglancer
#'   clearly shows voxel coordinates of 4.3x4.3x45. But in this function, the
#'   voxel coordinates must be set to 4.25 in x-y to give the correct answers.
#'
#' @param voxdims The voxel dimensions (in nm). See details.
#' @inheritParams fafbseg::flywire_xyz2id
#'
#' @return A character vector of segment ids, NA when lookup fails.
#' @family banc-ids
#' @seealso \code{\link{flywire_xyz2id}}
#' @export
#' @importFrom nat xyzmatrix
#' @examples
#' \dontrun{
#' # a point from neuroglancer, should map to 648518346498932033
#' banc_xyz2id(cbind(34495, 82783, 1954), rawcoords=TRUE)
#' }
banc_xyz2id <- function(xyz,
                        rawcoords=FALSE,
                        voxdims=c(4, 4, 45),
                        root=TRUE,
                        ...){
  if(is.numeric(xyz) && !is.matrix(xyz) && length(xyz)==3)
    xyz=matrix(xyz, ncol=3)
  if(rawcoords)
    xyz=scale(xyzmatrix(xyz), center = F, scale = 1/voxdims)
  svids=banc_supervoxels(xyz, voxdims=voxdims)
  if(isTRUE(root)) banc_rootid(svids) else svids
}

# rawxyz=cbind(34496, 82782, 1954)
# nmxyz=cbind(34496, 82782, 1954)*c(4.3,4.3,45)
banc_supervoxels <- function(x, voxdims=c(4,4,45)) {
  pts=scale(xyzmatrix(x), center = F, scale = voxdims)
  nas=rowSums(is.na(pts))>0
  if(any(nas)) {
    svids=rep("0", nrow(pts))
    svids[!nas]=banc_supervoxels(pts[!nas,,drop=F], voxdims = c(1,1,1))
    return(svids)
  }
  # URL is wrong
  u="https://services.itanna.io/app/transform-service/query/dataset/banc_v1/s/2/values_array_string_response"
  body=jsonlite::toJSON(list(x=pts[,1], y=pts[,2], z=pts[,3]))
  res=httr::POST(u, body = body)
  httr::stop_for_status(res)
  j=httr::content(res, as='text', encoding = 'UTF-8')
  svids=unlist(jsonlite::fromJSON(j, simplifyVector = T), use.names = F)
  svids
}

#' Check if a banc root id is up to date
#'
#' @inheritParams fafbseg::flywire_islatest
#' @param ... Additional arguments passed to \code{\link{flywire_islatest}}
#'
#' @export
#' @family banc-ids
#' @examples
#' banc_islatest("648518346473954669")
banc_islatest <- function(x, timestamp=NULL, ...) {
  with_banc(flywire_islatest(x=x, timestamp = timestamp, ...))
}


#' Find the latest id for a banc root id
#'
#' @inheritParams fafbseg::flywire_latestid
#' @param ... Additional arguments passed to \code{\link{flywire_latestid}}
#'
#' @export
#' @seealso \code{\link{banc_islatest}}
#' @family banc-ids
#' @examples
#' \dontrun{
#' banc_latestid("648518346473954669")
#' }
banc_latestid <- function(rootid, sample=1000L, cloudvolume.url=NULL, Verbose=FALSE, ...) {
  with_banc(fafbseg::flywire_latestid(rootid=rootid, sample = sample, Verbose=Verbose, ...))
}


#' Return a vector of banc root ids from diverse inputs
#'
#' @param x A data.frame, URL or vector of ids
#' @param integer64 Whether to return ids as \code{bit64::integer64} or
#'   character vectors. Default value of NA leaves the ids unmodified.
#'
#' @return A vector of ids
#' @export
#' @family banc-ids
#' @examples
#' banc_ids(data.frame(rootid="648518346474360770"))
banc_ids <- function(x, integer64=NA) {
  if(is.data.frame(x)) {
    colstocheck=c("rootid", "id", "pre_id", "post_id")
    for(col in colstocheck) {
      if(col %in% colnames(x))
        return(banc_ids(x[[col]], integer64 = integer64))
    }
    i64=sapply(x, bit64::is.integer64)
    if(sum(i64)==1)
      return(banc_ids(x[[which(i64)]], integer64 = integer64))
    stop("Unable to find a column containing ids!")
  }
  if(!all(fafbseg:::valid_id(x)))
    stop("Some ids are invalid!")
  if(isTRUE(integer64)) bit64::as.integer64(x)
  else if(isFALSE(integer64)) as.character(x)
  else x
}

#' Convert between banc cell ids and root ids
#'
#' @description Converts between BANC cell ids (should survive most edits) and
#'   root ids (guaranteed to match just one edit state). See details.
#'
#' @details CAVE/PyChunkedGraph assigns a 64 bit integer root id to all bodies
#'   in the segmentation. These root ids are persistent in a computer science
#'   sense, which is often the exact opposite of what neuroscientists might
#'   imagine. Specifically, a given root id is matched to a single edit state of
#'   a neuron. If the neuron is edited, then root id changes. In contrast, cell
#'   ids do not change even in the face of edits. However, it is important to
#'   understand that they correspond to a specific point on a neuron, commonly
#'   the nucleus. If the nucleus is edited away from a the rest of a neuron to
#'   which it previously belonged, then the cell id and any associated edits
#'   will effectively with move it.
#'
#'   For further details see
#'   \href{https://banc-reconstruction.slack.com/archives/CLDH21J4U/p1690755500802509}{banc
#'   slack} and
#'   \href{https://github.com/htem/banc_auto_recon/wiki/Neuron-annotations#neuron_information}{banc
#'   wiki}.
#'
#' @param rootids banc root ids in any form understood by
#'   \code{\link{banc_ids}}. The default value of NULL will return all cell ids.
#' @param cellids Integer cell ids between between 1 and around 20000 that
#'   \emph{should} uniquely identify each cell in the dataset.
#' @param timestamp An optional time stamp. You should give only one of
#'   \code{version} or \code{timestamp}. When both are missing, ids should match
#'   the live materialisation version including up to the second edits.
#' @param version An optional integer CAVE materialisation version. You should
#'   give only one of \code{version} or \code{timestamp}. When both are missing,
#'   ids should match the live materialisation version including up to the
#'   second edits.
#' @param rval Whether to return the cell ids or the whole of the CAVE table
#'   with additional columns.
#' @param cellid_table Optional name of cell id table (the default value of
#'   \code{NULL} should find the correct table).
#' @return Either a vector of ids or a data.frame depending on \code{rval}. For
#'   cell ids the vector will be an integer for root ids (segment ids), a
#'   character vector or an \code{bit64::integer64} vector depending on the
#'   \code{integer64} argument.
#' @inheritParams banc_ids
#' @family banc-ids
#' @export
#'
#' @examples
#' \donttest{
#' banc_cellid_from_segid(banc_latestid("720575941480769421"))
#' }
banc_cellid_from_segid <- function(rootids=NULL, timestamp=NULL, version=NULL, cellid_table = NULL, rval=c("ids", 'data.frame')) {
  rval=match.arg(rval)
  if(is.null(cellid_table))
    cellid_table=banc_cellid_table()
  if(!is.null(rootids)) {
    rootids=banc_ids(rootids, integer64=F)
  idlist=list(pt_root_id=rootids)
  } else idlist=NULL
  live=is.null(timestamp) && is.null(version)
  res=banc_cave_query(table = cellid_table, timestamp=timestamp,
                      version=version, filter_in_dict=idlist, live=live)
  if(is.null(rootids)) {
    if(rval=='ids') {
      banc_ids(res[['id']], integer64 = F)
    } else res
  }
  ids64=banc_ids(rootids, integer64=T)
  if(!all(found <- ids64 %in% res$pt_root_id)) {
    warning(sum(!found), "/", length(rootids), " could not be found!")
  }
  if(rval=='ids') {
    res[['id']][match(rootids, res[['pt_root_id']])]
  } else res
}

#' @rdname banc_cellid_from_segid
#' @export
#'
#' @examples
#' \donttest{
#' banc_cellid_from_segid(banc_latestid("720575941480769421"))
#' }
banc_segid_from_cellid <- function(cellids=NULL, timestamp=NULL, version=NULL, rval=c("ids", 'data.frame'), integer64=FALSE, cellid_table = NULL) {
  rval=match.arg(rval)
  if(is.null(cellid_table))
    cellid_table=banc_cellid_table()
  if(!is.null(cellids)) {
    cellids <- checkmate::assert_integerish(cellids, coerce = T)
    idlist=list(id=cellids)
  } else idlist=NULL
  live=is.null(timestamp) && is.null(version)
  res=banc_cave_query(table = cellid_table, timestamp=timestamp,
                      version=version, filter_in_dict=idlist, live=live)
  if(is.null(cellids)) {
    if(rval=='ids') {
      banc_ids(res[['pt_root_id']], integer64 = F)
    } else res
  }
  if(!all(found <- cellids %in% res[['id']])) {
    warning(sum(!found), "/", length(cellids), " could not be found!")
  }
  if(rval=='ids') {
    banc_ids(res[['pt_root_id']][match(cellids, res[['id']])], integer64 = integer64)
  } else res
}

# private function to return the latest cellids table
# this is a configurable option in the python package
banc_cellid_table <- memoise::memoise(function(fac=banc_cave_client()) {
  tables=fac$materialize$tables
  tablenames=names(tables)
  seltable=rev(sort(grep("cell_ids", tablenames, value = T)))[1]
  return(seltable)
})
