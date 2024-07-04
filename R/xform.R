#' Set Default View for BANC EM Dataset
#'
#' @description
#' This function sets a default view for visualizing the 'BANC' Electron Microscopy (EM) dataset
#' using the rgl package. It adjusts the viewpoint to a specific orientation and zoom level
#' that is optimal for viewing this particular dataset.
#'
#' @details
#' The function uses `rgl::rgl.viewpoint()` to set a predefined user matrix and zoom level.
#' This matrix defines the rotation and translation of the view, while the zoom parameter
#' adjusts the scale of the visualization.
#'
#' @return
#' This function is called for its side effect of changing the rgl viewpoint.
#' It does not return a value.
#'
#' @examples
#' \dontrun{
#' # Assuming you have already plotted your BANC EM data
#' banc_view()
#' }
#'
#' @note
#' This function assumes that an rgl device is already open and that the BANC EM dataset
#' has been plotted. It will not create a new plot or open a new rgl device.
#'
#' @seealso
#' \code{\link[rgl]{rgl.viewpoint}} for more details on setting viewpoints in rgl.
#'
#' @export
banc_view <- function(){
  rgl::rgl.viewpoint(userMatrix  = banc_rotation_matrices[["main"]], zoom = 0.82)
}

# for nm
#' @export
#' @rdname banc_view
banc_side_view <- function(){
  rgl::rgl.viewpoint(userMatrix = banc_rotation_matrices[["side"]], zoom = 0.29)
}

# for nm
#' @export
#' @rdname banc_view
banc_front_view <- function(){
  rgl::rgl.viewpoint(userMatrix = banc_rotation_matrices[["front"]], zoom = 0.62)
}

# for nm
#' @export
#' @rdname banc_view
banc_vnc_view <- function(){
  rgl::rgl.viewpoint(userMatrix = banc_rotation_matrices[["vnc"]], zoom = 0.51)
}

# hidden
banc_rotation_matrices <- list(
  main = structure(c(0.961547076702118, 0.037275392562151,
                                  0.27209860086441, 0, 0.0369537360966206, -0.999296963214874,
                                  0.00630810856819153, 0, 0.272142440080643, 0.00398948788642883,
                                  -0.962248742580414, 0, 0, 0, 0, 1), dim = c(4L, 4L)),
  side = structure(c(0.188666880130768, 0.137750864028931,
                     -0.972331881523132, 0, 0.130992725491524, -0.98479551076889,
                     -0.114099271595478, 0, -0.97326534986496, -0.105841755867004,
                     -0.203842639923096, 0, 0, 0, 0, 1), dim = c(4L, 4L)),
  front = structure(c(0.99931389093399, 0.0139970388263464,
                      -0.0342894680798054, 0, -0.0321401171386242, -0.132316529750824,
                      -0.990686297416687, 0, -0.0184037387371063, 0.991108655929565,
                      -0.131775915622711, 0, 0, 0, 0, 1), dim = c(4L, 4L)),
  vnc = structure(c(0.159858450293541, -0.951453745365143,
                                 0.263022243976593, 0, -0.95634800195694, -0.0832427442073822,
                                 0.280123054981232, 0, -0.244629606604576, -0.296320915222168,
                                 -0.923228204250336, 0, 169877.109375, 8134.845703125, -597.831604003906,
                                 1), dim = c(4L, 4L)))



