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

#' @export
#' @rdname banc.surf
"banc_neuropil.surf"

#' @export
#' @rdname banc.surf
"banc_brain_neuropil.surf"

#' @export
#' @rdname banc.surf
"banc_vnc_neuropil.surf"

#' @export
#' @rdname banc.surf
"banc_neck_connective.surf"

## How it was obtained:
# res <- httr::GET("https://www.googleapis.com/storage/v1/b/zetta_lee_fly_cns_001_kisuk/o/final%2Fv2%2Fvolume_meshes%2Fmeshes%2F1%3A0.drc?alt=media&neuroglancer=610000b05b6497edcf20b78f29516970")
# httr::stop_for_status(res)
# bytes <- httr::content(res, as = "raw")
# banc.mesh <- malevnc:::decode_neuroglancer_mesh(bytes)
# banc.surf <- as.hxsurf(banc.mesh)
# banc.surf$Vertices[,"Z"] <- banc.surf$Vertices[,"Z"]*0.9462 # scaling Jasper worked out
# save(banc.surf, file="data/BANC.surf.rda")
