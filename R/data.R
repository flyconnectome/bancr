#' Simplified tissue and neuropil surfaces for BANC
#'
#' @name banc.surf
#' @docType data
#' @description `banc.surf` is unsymmetrical and not a normalised version of the mesh.
#' It is the outline of the dataset in nanometers.`banc_neuropil.surf` represents the synaptic neuropil.
#' Built from the BANC synapse cloud, but not optimised to include 100% of bona fide presynapses.
#' `banc_neuropils.surf` contains the standard Ito et al., 2014 brain neuropil volumes transformed into BANC space.
#' `banc_al.surf` contains the standard Bates and Schlegel et al, 2021 right antennal lobe glomeruli brain neuropil volumes transformed into BANC space.
#' These neuropils may also be seen in neuroglancer,
#' \href{https://spelunker.cave-explorer.org/#!middleauth+https://global.daf-apis.com/nglstate/api/v1/6237404072509440}{here}.
#'
#' @examples
#' \dontrun{
#' # Depends on nat
#' library(nat)
#' rgl::wire3d(banc.surf)
#' }
"banc.surf"

#' @docType data
#' @rdname banc.surf
"banc_neuropil.surf"

#' @docType data
#' @rdname banc.surf
"banc_brain_neuropil.surf"

#' @docType data
#' @rdname banc.surf
"banc_vnc_neuropil.surf"

#' @docType data
#' @rdname banc.surf
"banc_neck_connective.surf"

#' @docType data
#' @rdname banc.surf
"banc_brain_neuropils.surf"

#' @docType data
#' @rdname banc.surf
"banc_al.surf"

#' @docType data
#' @rdname banc.surf
"banc_vnc_neuropils.surf"

#' @docType data
#' @rdname banc.surf
"banc_vnc_tracts.surf"

#' @docType data
#' @rdname banc.surf
"banc_vnc_nerves.surf"

#' @docType data
#' @rdname banc.surf
"banc_hemibrain.surf"

## How it was obtained:
# res <- httr::GET("https://www.googleapis.com/storage/v1/b/zetta_lee_fly_cns_001_kisuk/o/final%2Fv2%2Fvolume_meshes%2Fmeshes%2F1%3A0.drc?alt=media&neuroglancer=610000b05b6497edcf20b78f29516970")
# httr::stop_for_status(res)
# bytes <- httr::content(res, as = "raw")
# banc.mesh <- malevnc:::decode_neuroglancer_mesh(bytes)
# banc.surf <- as.hxsurf(banc.mesh)
# banc.surf$Vertices[,"Z"] <- banc.surf$Vertices[,"Z"]*0.9462 # scaling Jasper worked out
# usethis::use_data(banc.surf, overwrite = TRUE)

#' Thin-Plate Spline Registration from BANC to JRC2018F template brain
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

#' @docType data
#' @rdname banc_to_jrc2018f_tpsreg
"jrc2018f_to_banc_tpsreg"

#' @docType data
#' @rdname banc_to_jrc2018f_tpsreg
"jrcvnc2018f_to_banc_tpsreg"

#' @docType data
#' @rdname banc_to_jrc2018f_tpsreg
"banc_to_jrcvnc2018f_tpsreg"

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

#' BANC neuropil name to number correspondence for neuroglancer
#'
#' @name banc.surf
#' @docType data
#' @description A BANC neuroglaner scene can be directed to a google cloud storage
#' location, where BANC-transformed standard neuropil meshes reside.
#' The source is
#' `precomputed://gs://lee-lab_brain-and-nerve-cord-fly-connectome/volume_meshes`
#' They can be
#' plotted in neuroglancer by adding
#' \href{https://spelunker.cave-explorer.org/#!middleauth+https://global.daf-apis.com/nglstate/api/v1/6237404072509440}{this}
#' location, entering the `Seg.` pane
#' and entering the number or name that corresponds to the correct mesh.
#' This data frame gives the mesh name to number correspondences.
#'
#' @seealso
#' \code{\link{banc.surf}} for the available neuropil objects for BANC.
#' These are `hxsruf` objects, names for subregions can be found as so:
#' `banc_brain_neuropil.surf$RegionList`
"banc_volumes.df"

#' User information (name + CAVE ID) for active BANC users
#'
#' @name banc_users
#' @docType data
#' @description The purpose of this table is to map CAVE users IDs to names, in order to credit annotation work done in BANC CAVE.
#' This information is based on \href{https://docs.google.com/spreadsheets/d/1UFmeWr2uF9jTLVMw3bD6nM3ejM-b-HDZz6sQBPTEoZ8/edit?gid=1163959922#gid=1163959922}{google sheet}.
#'
#' @examples
#' \dontrun{
#' View(banc_users)
#' }
"banc_users"

# banc_users <- googlesheets4::read_sheet("1UFmeWr2uF9jTLVMw3bD6nM3ejM-b-HDZz6sQBPTEoZ8")
# colnames(banc_users) <- snakecase::to_snake_case(colnames(banc_users))
# banc_users <- banc_users %>%
#   dplyr::select(name, pi_lab, cave_id) %>%
#   dplyr::mutate(name = gsub("\\(.*","",name)) %>%
#   dplyr::mutate(pi_lab = gsub("\\(.*","",pi_lab)) %>%
#   dplyr::mutate(pi_lab = ifelse(is.na(pi_lab),name,pi_lab)) %>%
#   dplyr::mutate(pi_lab = ifelse(pi_lab=="PI",paste0(name," lab"),pi_lab)) %>%
#   dplyr::arrange(pi_lab, name) %>%
#   dplyr::distinct(cave_id, .keep_all = TRUE) %>%
#   dplyr::filter(!is.na(pi_lab))
# usethis::use_data(banc_users, overwrite = TRUE)

