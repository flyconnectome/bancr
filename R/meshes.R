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
#' neuron.mesh <- read_banc_meshes("720575941650432785")
#' plot3d(neuron.mesh, alpha = 0.1)
#' nucleus.mesh <- read_banc_meshes("72903876004544795")
#' plot3d(nucleus.mesh, col = "black")
#' }
read_banc_meshes <- function(ids, savedir=NULL, format=c("ply", "obj"), ...) {
  format=match.arg(format)
  read_cloudvolume_meshes(ids, savedir = savedir, cloudvolume.url = banc_cloudvolume_url(set=FALSE), format=format, ...)
}

#' @export
banc_read_nuclei_mesh <- function(ids, lod = 1L, savedir=NULL,  method=c('vf', 'ply'), ...) {
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
#' @param invert if \code{FALSE} returns brain points, if \code{TRUE} returns VNC points.
#'
#' @return Remove points above or below the midsection of the neck connective of BANC.
#' @seealso \code{\link{banc.surf}}
#' @export
#' @examples
#' \donttest{
#' m = read_banc_meshes("720575941650432785")
#' m.brain = banc_decapitate(m)
#' m.vnc = banc_decapitate(m, invert = TRUE)
#' plot3d(m.brain, col = "red")
#' plot3d(m.brain, col = "cyan")
#' plot3d(banc.surf, col = "grey", alpha = 0.1)
#' }
banc_decapitate <- function(x, invert = FALSE, reference = "BANC"){
  y.cut <- 5e+05
  v1 <- c(0,y.cut,0)
  v2 <- c(1e10,y.cut,0)
  v3 <- c(0,y.cut,1e10)
  if (reference!="BANC"){
    v1 <- nat.templatebrains::xform_brain(v1, sample = "BANC", reference = reference)
    v2 <- nat.templatebrains::xform_brain(v2, sample = "BANC", reference = reference)
    v3 <- nat.templatebrains::xform_brain(v3, sample = "BANC", reference = reference)
    y.cut <- v1[,2]
  }
  ismesh <- any(class(x[[1]]), class(x)) %in% "mesh3d"
  if(!ismesh & any(class(x)%in%c("neuron","neuronlist","mesh3d"))){
    if (invert){
      z <- subset(x,x$y<y.cut)
    }else{
      z <- subset(x,x$y>y.cut)
    }
  }else if (ismesh){
    if (any(class(x)=="mesh3d")){
      z <- Morpho::cutMeshPlane(x,
                                v1=v1, v2=v2, v3=v3,
                                normal = NULL, keep.upper = invert)
    }else if (any(class(x)=="neuronlist")){
      z <- nat::nlapply(x, Morpho::cutMeshPlane,
                        v1=v1, v2=v2, v3=v3,
                        normal = NULL, keep.upper = invert)
    }else{
      z <- lapply(x, Morpho::cutMeshPlane,
                  v1=v1, v2=v2, v3=v3,
                  normal = NULL, keep.upper = invert)
    }
  }else{
    ydim <- intersect(colnames(x),c("Y","y"))[[1]]
    if (invert){
      z <- x[x[,ydim]<y.cut,]
    }else{
      z <- x[x[,ydim]>y.cut,]
    }
  }
  z
}

#' Read BANC euroglancer meshes, e.g., ROI meshes
#'
#' @param x the numeric identifier that specifies the mesh to read, defaults to \code{1} the BANC outline mesh.
#' @param url the URL that directs \code{bancr} to where BANC meshes are stored.
#' @return a mesh3d object for the specified mesh.
#' @export
#' @seealso \code{\link{read_banc_meshes}}
#' @examples
#' \dontrun{
#' banc.mesh  <- banc_read_neuroglancer_mesh()
#' }
banc_read_neuroglancer_mesh <- function(x = 1,
                                        url="https://www.googleapis.com/storage/v1/b/zetta_lee_fly_cns_001_kisuk/o/final%2Fv2%2Fvolume_meshes%2Fmeshes%2F{x}%3A0.drc?alt=media&neuroglancer=610000b05b6497edcf20b78f29516970",
                                        ...){
  completed_url <- glue::glue(url, x=x)
  res <- httr::GET(completed_url, ...)
  httr::stop_for_status(res)
  bytes <- httr::content(res, as = "raw")
  decode_neuroglancer_mesh(bytes)
}

# convert cloudvolume python mesh to an R mesh3d object
# method vf just uses the vertex and face arrays
# ply writes out to Stanford ply format and reads back in again
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

# from package: malevnc
decode_neuroglancer_mesh <- function (bytes, format = c("mesh3d", "raw")) {
  format = match.arg(format)
  con = rawConnection(bytes)
  on.exit(close(con))
  nverts = readBin(con, what = "int", size = 4, n = 1)
  verts = readBin(con, what = "numeric", n = nverts * 3, size = 4)
  nidxs = length(bytes)/4 - 1L - length(verts)
  idx = readBin(con, what = "int", n = nidxs, size = 4)
  if (format == "raw") {
    structure(list(v = matrix(verts, ncol = 3, byrow = T),
                   i = matrix(idx, ncol = 3, byrow = T)), class = "ngmesh")
  }
  else {
    rgl::tmesh3d(matrix(verts, nrow = 3, byrow = F), matrix(idx +
                                                              1L, nrow = 3, byrow = F), homogeneous = F)
  }
}




