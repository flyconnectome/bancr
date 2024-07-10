#' Simplified tissue surface of BANC
#'
#' @name banc.surf
#' @docType data
#' @description This is unsymmetrical and not a normalised version of the mesh.
#' It is the outline of the dataset in native voxel space.
#'
#' @examples
#' \dontrun{
#' # Depends on nat
#' library(nat)
#' rgl::wire3d(banc.surf)
#' }
"banc.surf"

#' @rdname banc.surf
"banc_neuropil.surf"

#' @rdname banc.surf
"banc_brain_neuropil.surf"

#' @rdname banc.surf
"banc_vnc_neuropil.surf"

#' @rdname banc.surf
"banc_neck_connective.surf"

## How it was obtained:
# res <- httr::GET("https://www.googleapis.com/storage/v1/b/zetta_lee_fly_cns_001_kisuk/o/final%2Fv2%2Fvolume_meshes%2Fmeshes%2F1%3A0.drc?alt=media&neuroglancer=610000b05b6497edcf20b78f29516970")
# httr::stop_for_status(res)
# bytes <- httr::content(res, as = "raw")
# banc.mesh <- malevnc:::decode_neuroglancer_mesh(bytes)
# banc.surf <- as.hxsurf(banc.mesh)
# banc.surf$Vertices[,"Z"] <- banc.surf$Vertices[,"Z"]*0.9462 # scaling Jasper worked out
# usethis::use_data(banc.surf, overwrite = TRUE)

#' Thin-Plate Spline Registration from BANC to JRC2018F templatebrain
#'
#' @description
#' A thin-plate spline (TPS) registration object that transforms 3D points
#' from the BANC nanometer coordinate system
#' to the D. melanogaster template brain JRC2018F coordinate system.
#'
#' @format An object of class \code{tpsCoeff} created using \code{Morpho::computeTransform}.
#'   It contains the following components:
#'   \describe{
#'     \item{Lw}{Matrix of TPS coefficients}
#'     \item{refmat}{Reference matrix (source landmarks)}
#'     \item{tarmat}{Target matrix (target landmarks)}
#'     \item{lattice}{3D array representing the deformation grid}
#'     \item{lambda}{Smoothing parameter used in the TPS computation}
#'     \item{scale}{Logical indicating whether scaling was used}
#'     \item{reflection}{Logical indicating whether reflection was allowed}
#'   }
#'
#' @details
#' This TPS registration was computed based on landmark correspondences
#' derived from an Elastix registration between the BANC and JRC2018F spaces.
#' It provides a smooth, interpolated transformation for any point in the BANC space
#' to its corresponding location in the JRC2018F space.
#'
#' @usage data(banc_to_jrc2018f_tpsreg)
#'
#' @source
#' Derived from Elastix registration results using the \code{banc_to_JRC2018F}
#' function and landmark correspondences extracted from that registration.
#'
#' @examples
#' \dontrun{
#' data(banc_to_jrc2018f_tpsreg)
#'
#' # Transform BANC points to JRC2018F using the TPS registration
#' banc_points <- matrix(rnorm(300), ncol=3)  # Example BANC points
#' jrc2018f_points <- Morpho::tps3d(banc_points, banc_to_jrc2018f_tpsreg)
#' }
#'
#' @seealso
#' \code{\link{banc_to_JRC2018F}} for the function that performs transformations
#' between BANC and JRC2018F spaces.
#' \code{\link[Morpho]{computeTransform}} for details on creating tpsCoeff objects.
#'
#' @name banc_to_jrc2018f_tpsreg
#' @docType data
"banc_to_jrc2018f_tpsreg"

#' @rdname banc_to_jrc2018f_tpsreg
"jrc2018f_to_banc_tpsreg"

#' Thin-Plate Spline Registration for Mirroring in BANC Space
#'
#' @description
#' A thin-plate spline (TPS) registration object that mirrors 3D points
#' within the BANC (Buhmann et al. Adult Neural Connectome) coordinate system.
#'
#' @format An object of class \code{tpsCoeff} created using \code{Morpho::computeTransform}.
#'   It contains the following components:
#'   \describe{
#'     \item{Lw}{Matrix of TPS coefficients}
#'     \item{refmat}{Reference matrix (source landmarks)}
#'     \item{tarmat}{Target matrix (mirrored landmarks)}
#'     \item{lattice}{3D array representing the deformation grid}
#'     \item{lambda}{Smoothing parameter used in the TPS computation}
#'     \item{scale}{Logical indicating whether scaling was used}
#'     \item{reflection}{Logical indicating whether reflection was allowed}
#'   }
#'
#' @details
#' This TPS registration was computed to allow mirroring of points directly within
#' the BANC coordinate system. It provides a smooth, interpolated transformation
#' for any point in the BANC space to its mirrored counterpart, accounting for
#' any asymmetries in the BANC reference brain.
#'
#' @usage data(banc_mirror_tpsreg)
#'
#' @source
#' Derived from landmark correspondences between original and mirrored points
#' in the BANC space, possibly utilizing transformations to and from the
#' JRC2018F space for accurate mirroring.
#'
#' @examples
#' \dontrun{
#' data(banc_mirror_tpsreg)
#'
#' # Mirror BANC points using the TPS registration
#' banc_points <- matrix(rnorm(300), ncol=3)  # Example BANC points
#' mirrored_points <- Morpho::tps3d(banc_points, banc_mirror_tpsreg)
#' }
#'
#' @seealso
#' \code{\link{banc_mirror}} for the function that performs mirroring
#' of points in BANC space.
#' \code{\link{banc_to_jrc2018f_tpsreg}} for the TPS registration between
#' BANC and JRC2018F spaces.
#' \code{\link[Morpho]{computeTransform}} for details on creating tpsCoeff objects.
#'
#' @name banc_mirror_tpsreg
#' @docType data
"banc_mirror_tpsreg"




