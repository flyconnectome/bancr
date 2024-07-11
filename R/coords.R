#' Handle raw and nm calibrated banc coordinates
#'
#' @description \code{banc_voxdims} returns the image voxel dimensions which
#'   are normally used to scale between \bold{raw} and \bold{nm} coordinates.
#'
#' @param url Optional neuroglancer URL containing voxel size. Defaults to
#'   \code{getOption("fafbseg.sampleurl")} as set by
#'   \code{\link{choose_banc}}.
#'
#' @return For \code{banc_voxdims} A 3-vector
#' @export
#'
#' @examples
#' banc_voxdims()
banc_voxdims <- memoise::memoise(function(url=choose_banc(set=FALSE)[['fafbseg.sampleurl']]) {
  fafbseg::flywire_voxdims(url)
})
# 4  4 40

#' @param x 3D coordinates in any form compatible with \code{\link{xyzmatrix}}
#'
#' @return for \code{banc_raw2nm} and \code{banc_nm2raw} an Nx3 matrix of
#'   coordinates
#' @param vd The voxel dimensions in nm. Expert use only. Normally found
#'   automatically.
#' @export
#' @rdname banc_voxdims
#' @details relies on nat >= 1.10.4
#' @examples
#' banc_raw2nm(c(159144, 22192, 3560))
#' banc_raw2nm('159144 22192 3560')
#' \dontrun{
#' banc_nm2raw(clipr::read_clip())
#' }
banc_nm2raw <- function(x, vd=banc_voxdims()) fancr::fanc_nm2raw(x, vd=vd)

#' @export
#' @rdname banc_voxdims
banc_raw2nm <- function(x, vd=banc_voxdims()) fancr::fanc_raw2nm(x, vd=vd)
