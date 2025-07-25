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
banc_rootid <- function(x,
                        integer64 = FALSE,
                        version = NULL,
                        timestamp = NULL,
                        ...) {
  if(!is.null(version)){
    timestamp <- with_banc(fafbseg:::flywire_timestamp(version = version,
                                                      convert = FALSE))
  }
  rid = flywire_rootid(
    x = x,
    integer64 = integer64,
    method = "cloudvolume",
    # agglomerate = T,
    timestamp = timestamp,
    cloudvolume.url = banc_cloudvolume_url(),
    ...
  )
  rid
}
# timestamp = with_banc(fafbseg:::flywire_timestamp(version = "612",
#                                                   timestamp = NULL,
#                                                   convert = FALSE))
# with_banc(fafbseg:::flywire_roots_helper("76562155939424493", method = "cloudvolume",
#                                          cloudvolume.url =  banc_cloudvolume_url(),
#                                          integer64 = FALSE,
#                                          timestamp = timestamp,
#                                          stop_layer = NULL))

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
#' svids=banc_leaves("720575941478275714")
#' head(svids)
#' }
banc_leaves <- function(x, integer64=TRUE, ...) {
  svids=with_banc(flywire_leaves(x=x, integer64 = integer64, ...))
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
#' @inheritParams fafbseg::flywire_xyz2id
#'
#' @return A character vector of segment ids, NA when lookup fails.
#' @family banc-ids
#' @seealso \code{\link{flywire_xyz2id}}
#' @export
#' @importFrom nat xyzmatrix
#' @examples
#' \dontrun{
#' # a point from neuroglancer, should map to 720575941623125868
#' banc_xyz2id(cbind(438976,985856,215955),  version="282", rawcoords=FALSE)
#'
#' # Get root ID for an older materialization
#' banc_xyz2id(cbind(462572, 370000, 134955), version="612", root=TRUE)
#'
#' # Get the most recent root ID
#' banc_xyz2id(cbind(462572, 370000, 134955),root=TRUE)
#' }
banc_xyz2id <- function(xyz,
                        rawcoords = FALSE,
                        root = TRUE,
                        integer64 = FALSE,
                        fast_root = TRUE,
                        method = c("cloudvolume", "spine"),
                        version = NULL,
                        ...) {
  fafbseg:::check_cloudvolume_reticulate()
  voxdims <- banc_voxdims()
  method <- match.arg(method)
  if (isTRUE(is.numeric(xyz) && is.vector(xyz) && length(xyz) ==
             3)) {
    xyz <- matrix(xyz, ncol = 3)
  }
  else {
    xyz <- xyzmatrix(xyz)
  }
  if (isTRUE(rawcoords)) {
    xyz <- scale(xyz, scale = 1/voxdims, center = FALSE)
  }
  checkmate::assertNumeric(xyz)
  if (method %in% c("spine")) {
    res <- banc_supervoxels(xyz, voxdims=voxdims)
  }
  else {
    cv <- banc_cloudvolume()
    pycode <- sprintf("\nfrom cloudvolume import Vec\n\ndef py_flywire_xyz2id(cv, xyz, agglomerate):\n  pt = Vec(*xyz) // cv.meta.resolution(0)\n  img = cv.download_point(pt, mip=0, size=1, agglomerate=agglomerate)\n  return str(img[0,0,0,0])\n")
    pydict <- reticulate::py_run_string(pycode)
    safexyz2id <- function(pt) {
      tryCatch(pydict$py_flywire_xyz2id(cv, pt, agglomerate = root &&
                                          !fast_root), error = function(e) {
                                            warning(e)
                                            NA_character_
                                          })
    }
    res = pbapply::pbapply(xyz, 1, safexyz2id, ...)
  }
  if (fast_root && root) {
    res = banc_rootid(res,
                     integer64 = integer64,
                     version = version,
                     ...)
  }
  if (isFALSE(rawcoords) && sum(res == 0) > 0.25 * length(res)) {
    if (all(banc_is_rawcoord(xyz))) {
      warning("It looks like you may be passing in raw coordinates. If so, use rawcoords=TRUE")
    }
  }
  if (integer64)
    bit64::as.integer64(res)
  else as.character(res)
}

# hidden
banc_is_rawcoord <- function (xyz){
  rawbb = structure(c(19848.22, 241544.98, 8881.11, 282792.41, -1.4, 7010.37
  ), dim = 2:3, class = "boundingbox")
  pointsinside(xyz, rawbb)
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
#' banc_islatest("720575941520182775")
banc_islatest <- function(x, timestamp=NULL, ...) {
  with_banc(flywire_islatest(x=x, timestamp = timestamp, ...))
}

#' Find the latest id for a banc root id
#'
#' @inheritParams fafbseg::flywire_latestid
#' @param x a `data.frame` with at least one of: `root_id`, `pt_root_id`, `supervoxel_id` and/or `pt_supervoxel_id`.
#' Supervoxels will be preferentially used to update the `root_id` column.
#' Else a vector of `BANC` root IDs.
#' @param ... Additional arguments passed to \code{\link{flywire_latestid}}
#' @param root.column when `x` is a `data.frame`, the `root_id` column you wish to update
#' @param supervoxel.column when `x` is a `data.frame`, the `supervoxel_id` column you wish to use to update `root.column`
#' @param position.column when `x` is a `data.frame`, the `position` column with xyz values you wish to use to update `supervoxel.column`
#' @param use.cave read from the best established CAVE tables and join by `pt_supervoxel_id` to update `root_id`
#' @param serial if TRUE and x is a vector, calls `banc_updateids` on each ID in sequence to bufffer against connection failures. Slower.
#' @export
#' @seealso \code{\link{banc_islatest}}
#' @family banc-ids
#' @examples
#' \dontrun{
#' banc_latestid("720575941520182775")
#' }
banc_latestid <- function(rootid, sample=1000L, cloudvolume.url=NULL, Verbose=FALSE, ...) {
  with_banc(fafbseg::flywire_latestid(rootid=rootid, sample = sample, Verbose=Verbose, ...))
}

#' @export
#' @rdname banc_latestid
banc_updateids <- function(x,
                           root.column = "root_id",
                           supervoxel.column = "supervoxel_id",
                           position.column = "position",
                           use.cave = TRUE,
                           serial = FALSE,
                           ...){
  if(is.data.frame(x)){

    # Use CAVE tables to join by supervoxel_id
    if(use.cave&!is.null(supervoxel.column)){
      if(supervoxel.column%in%colnames(x)){
        cat('joining to CAVE tables ...\n')
        proofed <- banc_backbone_proofread() %>% dplyr::distinct(pt_root_id, pt_supervoxel_id)
        info <- banc_cell_info()  %>% dplyr::distinct(pt_root_id, pt_supervoxel_id)
        nuclei <- banc_nuclei()  %>% dplyr::distinct(pt_root_id = root_id, pt_supervoxel_id)
        nerves <- banc_peripheral_nerves() %>% dplyr::distinct(pt_root_id, pt_supervoxel_id)
        seeds <- banc_neck_connective_neurons() %>% dplyr::distinct(pt_root_id, pt_supervoxel_id)
        cave.tables <- proofed %>%
          rbind(info) %>%
          rbind(nuclei) %>%
          rbind(nerves) %>%
          rbind(seeds) %>%
          dplyr::mutate(cave_pt_root_id=as.character(pt_root_id),
                        cave_pt_supervoxel_id=as.character(pt_supervoxel_id)) %>%
          dplyr::distinct(cave_pt_root_id, cave_pt_supervoxel_id) %>%
          as.data.frame()
        x$cave_pt_supervoxel_id = x[[supervoxel.column]]
        x <- x %>%
          dplyr::left_join(cave.tables, by = "cave_pt_supervoxel_id")
        x[[root.column]] <- ifelse(is.na(x$cave_pt_root_id)|x$cave_pt_root_id=="0",x[[root.column]],x$cave_pt_root_id)
        x <- x %>%
          dplyr::select(-cave_pt_root_id, -cave_pt_supervoxel_id)
      }
    }

    # Update supervoxel IDs
    if(!is.null(position.column)){
      if(all(c(position.column,supervoxel.column)%in%colnames(x))){
        no.sp <- is.na(x[[supervoxel.column]])|x[[supervoxel.column]]=="0"
        if(sum(no.sp)){
          cat('determining missing supervoxel_ids ...\n')
          x[no.sp,][[supervoxel.column]] <- unname(pbapply::pbsapply(x[no.sp,][[position.column]], function(row){
            tryCatch(quiet_function(banc_xyz2id(row,rawcoords = TRUE, root = FALSE, ...)),
                     error = function(e) NA)
          }))
        }
      }
    }

    # what needs updating?
    if(!length(root.column)){
      root.column <- "root_id"
    }
    if(root.column%in%colnames(x)){
      cat('determining old root_ids...\n')
      old <- !banc_islatest(x[[root.column]], ...)
    }else{
      old <- rep(TRUE,nrow(x))
    }
    old[is.na(old)] <- TRUE
    message("old root_ids: ",sum(old),"\n")
    if(!sum(old)){
      return(x)
    }

    # update based on supervoxels
    if(!is.null(supervoxel.column)){
      if(supervoxel.column%in%colnames(x)){
        cat('updating root_ids with a supervoxel_id...\n')
        update <- unname(pbapply::pbsapply(x[old,][[supervoxel.column]], banc_rootid, ...))
        bad <- is.na(update)|update=="0"
        update <- update[!bad]
        if(length(update)) x[old,][[root.column]][!bad] <- update
        old[!bad] <- FALSE
      }
      old[is.na(old)] <- TRUE
    }

    # update based on position
    if(!is.null(position.column)){
      if(any(position.column%in%colnames(x)) && sum(old)){
        cat('updating root_ids with a position ...\n')
        update <- unname(pbapply::pbsapply(x[old,][[position.column]], function(row){
          tryCatch(quiet_function(banc_xyz2id(row,rawcoords = TRUE, root = TRUE, ...)),
                   error = function(e) NA)
        }))
        bad <- is.na(update)|update=="0"
        update <- update[!bad]
        if(length(update)) x[old,][[root.column]][!bad] <- update
        old[!bad] <- FALSE
      }
      old[is.na(old)] <- TRUE
    }

    # update based on root Ids
    if(root.column%in%colnames(x) && sum(old)){
      cat('updating root_ids without a supervoxel_id...\n')
      update <- banc_latestid(x[old,][[root.column]], ...)
      bad <- is.na(update)|update=="0"
      update <- update[!bad]
      if(length(update)) x[old,][[root.column]][!bad] <- update
      old[old][!bad] <- FALSE
    }
    old[is.na(old)] <- FALSE

  }else{
    if(serial){
      x <- pbapply::pbsapply(x, function(x) try(quiet_function(banc_updateids, serial = FALSE)))
    }else{
      cat('updating root_ids directly ...\n')
      old <- !banc_islatest(x, ...)
      old[is.na(old)] <- TRUE
      update <- banc_latestid(x[old], ...)
      bad <- is.na(update)|update=="0"
      update <- update[!bad]
      if(length(update)) x[old][!bad]  <- update
      old[!bad] <- FALSE
    }
  }

  # return
  if(sum(old)){
    warning("failed to update: ", sum(old),"\n")
  }
  x
}

# hidden
quiet_function <- function(...) {
  suppressMessages(
    suppressWarnings(
      capture.output(
        original_function(...),
        file = nullfile()
      )
    )
  )
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
    colstocheck=c("rootid", "id", "pre_id", "post_id", "root_id")
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

#' Convert between BANC cell ids and root ids
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
#' banc_cellid_from_segid(banc_latestid("720575941626035769"))
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
  seltable=rev(sort(grep("^cell_info", tablenames, value = T)))[1]
  return(seltable)
})

# hideen
banc_nucelus_id_to_rootid <- function(nucleus_ids){
  nuclei <- banc_nuclei()
  nuclei$root_id <- as.character(nuclei$root_id)
  nuclei$nucleus_id <- as.character(nuclei$nucleus_id)
  ids <- nuclei$root_id[match(nucleus_ids, nuclei$nucleus_id)]
  ids[is.na(ids)] <- "0"
  ids
}
