#' Read one or more BANC neuron and nuclei meshes
#'
#' @param ids One or more root ids
#' @param savedir An optional location to save downloaded meshes. This acts as a
#'   simple but effective cache since flywire neurons change id whenever they
#'   are edited.
#' @param lod The level of detail (highest resolution is 0, default of 2 gives a
#' good overall morphology while 3 is also useful and smaller still).
#' @param method How to treat the mesh object returned from neuroglancer, i.e. as
#' a \code{mesh3d} object or a \code{ply} mesh.
#' @param ... Additional arguments passed to
#'   \code{fafbseg::\link{read_cloudvolume_meshes}}
#' @inheritParams fafbseg::save_cloudvolume_meshes
#'
#' @return A \code{\link[nat]{neuronlist}} containing one or more \code{mesh3d}
#'   objects. See \code{nat::\link[nat]{read.neurons}} for details.
#' @export
#' @seealso \code{fafbseg::\link{read_cloudvolume_meshes}}
#' @examples
#' \donttest{
#' neuron.mesh <- banc_read_neuron_meshes("720575941478275714")
#' plot3d(neuron.mesh, alpha = 0.1)
#' nucleus.mesh <- banc_read_nuclei_mesh("72903876004544795")
#' plot3d(nucleus.mesh, col = "black")
#' }
banc_read_neuron_meshes <- function(ids, savedir=NULL, format=c("ply", "obj"), ...) {
  format=match.arg(format)
  with_banc(read_cloudvolume_meshes(ids, savedir = savedir, format=format, ...))
}

#' @export
#' @rdname banc_read_neuron_meshes
banc_read_nuclei_mesh <- function(ids, lod = 0L, savedir=NULL,  method=c('vf', 'ply'), ...) {
  cvu <- "precomputed://gs://lee-lab_brain-and-nerve-cord-fly-connectome/nuclei/seg_v1"
  cv <- fafbseg::flywire_cloudvolume(cloudvolume.url = cvu)
  li <- reticulate::py_eval(ids, convert = F)
  lod <- as.integer(lod)
  cm <- cv$mesh$get(li, lod = lod)
  if (!any(ids %in% names(cm))) {
    stop("Failed to read segid: ", ids)
  }
  cmesh <- reticulate::py_get_item(cm, li)
  m <- cvmesh2mesh(cmesh, method = method)
  m
}

#' Subset points to be in the brain or in the VNC
#'
#' @param x an object with 3d points to be subsetted, e.g. an xyz matrix, a \code{neuron}, \code{neuronlist} or a \code{mesh3d} object.
#' Points must be in native BANC space, i.e. plottable inside \code{banc.surf}.
#' @param y.cut Numeric, the Y-axis cut point, in nanometers, in BANC space,
#'  that separates the head from the neck and ventral nerve cord. For fitting to the
#'  MANC data set, a cut height of `y.cut=5e05` seems good.
#' @param invert if \code{FALSE} returns brain points, if \code{TRUE} returns VNC points.
#' @param ... Additional arguments passed to \code{\link{nlapply}} and then \code{\link{prune_vertices}}
#'
#' @return Remove points above or below the midsection of the neck connective of BANC.
#' @seealso \code{\link{banc.surf}}
#' @export
#' @examples
#' \donttest{
#' # DNa02
#' m = banc_read_neuron_meshes("720575941478275714")
#' m.brain = banc_decapitate(m)
#' m.vnc = banc_decapitate(m, invert = TRUE)
#' }
#' \dontrun{
#' plot3d(m.brain, col = "red")
#' plot3d(m.vnc, col = "cyan")
#' plot3d(banc.surf, col = "grey", alpha = 0.1)
#' }
banc_decapitate <- function(x, y.cut = 325000, invert = FALSE, ...) UseMethod('banc_decapitate')

#' @export
#' @rdname banc_decapitate
banc_decapitate.NULL <- function(x, y.cut = 325000, invert = FALSE, ...){
  #warning("banc_decapitate given a NULL object, returning NULL")
  NULL
}

#' @export
#' @rdname banc_decapitate
banc_decapitate.neuron <- function(x, y.cut = 325000, invert = FALSE, ...){
  z <- as.data.frame(nat::xyzmatrix(x))
  rownames(z) <- x$d$PointNo
  z <- subset(z,z$Y<y.cut)
  if(!nrow(z)){
   if(invert){
     warning("no points remaining after banc_decapitate.neuron")
     NULL
   }else{
     x
   }
  }else if ("splitneuron"%in%class(x)|"synapticneuron"%in%class(x)){
    prune_vertices.synapticneuron(x, verticestoprune = rownames(z), invert = invert, ...)
  }else{
    nat::prune_vertices(x, verticestoprune = rownames(z), invert = invert, ...)
  }
}

#' @export
#' @rdname banc_decapitate
banc_decapitate.neuronlist <- function(x, y.cut = 325000, invert = FALSE, ...){
  nat::nlapply(x, banc_decapitate, y.cut = y.cut, invert = invert, ...)
}

#' @export
#' @rdname banc_decapitate
banc_decapitate.matrix <- function(x, y.cut = 325000, invert = FALSE, ...){
  z <- nat::xyzmatrix(x)
  if(invert){
    z[z[,2]<y.cut,]
  }else{
    z[z[,2]>y.cut,]
  }
}

#' @export
#' @rdname banc_decapitate
banc_decapitate.data.frame <- function(x, y.cut = 325000, invert = FALSE, ...){
  z <- nat::xyzmatrix(x)
  if(invert){
    x[z[,2]<y.cut,]
  }else{
    x[z[,2]>y.cut,]
  }
}

#' @export
#' @rdname banc_decapitate
banc_decapitate.mesh3d <- function(x, y.cut = 325000, invert = FALSE, ...){
  v1 <- c(0,y.cut,0)
  v2 <- c(1e10,y.cut,0)
  v3 <- c(0,y.cut,1e10)
  Morpho::cutMeshPlane(x,
                      v1=v1, v2=v2, v3=v3,
                      normal = NULL, keep.upper = invert)
}

#' @export
#' @rdname banc_decapitate
banc_decapitate.hxsurf <- function(x, y.cut = 325000, invert = FALSE, ...){
  x <- rgl::as.mesh3d(x)
  m <- banc_decapitate.mesh3d(x=x, y.cut=y.cut, invert=invert, ...)
  nat::as.hxsurf(m)
}

# hidden, for now
banc_in_neuropil.mesh3d <- function(x,
                                    surf = bancr::banc_neuropil.surf,
                                    facenum = 2500,
                                    invert = FALSE){

  # Step 0, slightly scale the meshes
  m <- rgl::as.mesh3d(surf)
  # m1.scaled <- Morpho::scalemesh(as.mesh3d(bancr::banc_brain_neuropil.surf),
  #                               size = 1.2,
  #                               center = "mean")
  # m2.scaled <- Morpho::scalemesh(as.mesh3d(bancr::banc_vnc_neuropil.surf),
  #                               size = 1.2,
  #                               center = "mean")
  # m <-merge(m1.scaled,m2.scaled)

  # Step 1: Calculate distances
  distance_field <- Rvcg::vcgClostKD(x, m, sign = TRUE)
  if(invert){
    outside_vertices <- which(distance_field$quality > 0)
  }else{
    outside_vertices <- which(distance_field$quality <= 0)
  }
  C <- x

  # Step 2: Identify vertices to keep
  vertices_to_keep <- setdiff(1:ncol(x$vb), outside_vertices)

  # Step 3: Update the vertices
  C$vb <- x$vb[, vertices_to_keep, drop = FALSE]

  # Step 4: Create a mapping from old vertex indices to new ones
  vertex_map <- stats::setNames(seq_along(vertices_to_keep), vertices_to_keep)

  # Step 5: Update the faces (triangles)
  message("calculating faces ...")
  valid_faces <- pbapply::pbapply(C$it, 2, function(face) all(face %in% vertices_to_keep))
  C$it <- C$it[, valid_faces, drop = FALSE]
  message("calculating vertex map ...")
  C$it <- pbapply::pbapply(C$it, 2, function(face) vertex_map[as.character(face)])

  # Step 6: Remove any unused attributes
  C$normals <- NULL
  C$material <- NULL
  C$texcoords <- NULL

  # Step 7: Ensure the class is set correctly
  class(C) <- "mesh3d"
  D <- Rvcg::vcgClean(C, sel=c(0:8))
  D <- Rvcg::vcgIsolated(D, facenum = facenum)
  D
}

#' Read BANC euroglancer meshes, e.g., ROI meshes
#'
#' @param x the numeric identifier that specifies the mesh to read, defaults to \code{1} the BANC outline mesh.
#' @param url the URL that directs \code{bancr} to where BANC meshes are stored.
#' @param ... additional arguments to \code{\link{GET}}
#' @return a mesh3d object for the specified mesh.
#' @export
#' @seealso \code{\link{banc_read_neuron_meshes}}
#' @examples
#' \dontrun{
#' banc.mesh  <- banc_read_neuroglancer_mesh()
#' }
banc_read_neuroglancer_mesh <- function(x = 1,
                                        url = paste0(
                                          "https://www.googleapis.com/storage/v1/b/",
                                          "zetta_lee_fly_cns_001_kisuk/o/final%2Fv2%2F",
                                          "volume_meshes%2Fmeshes%2F{x}%3A0.drc?alt=media",
                                          "&neuroglancer=610000b05b6497edcf20b78f29516970"
                                        ),
                                        ...){
  completed_url <- glue::glue(url, x=x)
  res <- httr::GET(completed_url, ...)
  httr::stop_for_status(res)
  bytes <- httr::content(res, as = "raw")
  malevnc:::decode_neuroglancer_mesh(bytes)
}

# convert cloudvolume python mesh to an R mesh3d object
# method vf just uses the vertex and face arrays
# ply writes out to Stanford ply format and reads back in again
# from hemibrainr
cvmesh2mesh <- function(x, method=c('vf', 'ply'), ...) {
  method=match.arg(method)
  if(method=='vf') {
    verts=t(x$vertices)
    stopifnot(nrow(verts)==3)
    faces=t(x$faces+1)
    stopifnot(max(faces)<=ncol(verts))
    stopifnot(min(faces)>0)
    stopifnot(nrow(faces)==3)
    m=rgl::tmesh3d(vertices = verts, indices = faces, homogeneous = F)
  } else {
    bytes=x$to_ply()
    tf=tempfile(fileext = paste0('.', method))
    on.exit(unlink(tf))
    writeBin(bytes, con = tf)
    m=Rvcg::vcgPlyRead(tf, ...)
  }
  m
}

