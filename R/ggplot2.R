#' Compare two neurons from the BANC connectome dataset
#'
#' @description
#' This function creates a visual comparison of two neurons from the BANC (Bilateral Antennal lobe Neuron Connectome)
#' dataset. It generates four different anatomical views to facilitate comparison.
#'
#' @param neuron1 A neuronlist object representing the first neuron to be compared.
#' @param neuron2 A neuronlist object representing the second neuron to be compared.
#' @param neuron3 A neuronlist object representing the third neuron to be compared.
#' @param neuron1.info Character, a vector to be printed on the outputted ggplot in reference to neuron1.
#' @param neuron2.info Character, a vector to be printed on the outputted ggplot in reference to neuron2.
#' @param neuron3.info Character, a vector to be printed on the outputted ggplot in reference to neuron3.
#' @param banc_brain_neuropil A mesh object representing the brain neuropil. Default is banc_brain_neuropil.surf.
#' @param banc_vnc_neuropil A mesh object representing the ventral nerve cord (VNC) neuropil. Default is banc_vnc_neuropil.surf.
#' @param banc_neuropil A mesh object representing the entire neuropil. Default is banc_neuropil.surf.
#' @param filename Optional. If provided, the plot will be saved to this file. If NULL (default), the plot will be displayed but not saved.
#' @param width Numeric. The width of the output plot in inches. Default is 16.
#' @param height Numeric. The height of the output plot in inches. Default is 16.
#'
#' @return If filename is NULL, the function returns a ggplot object. If filename is provided, the function saves the plot and returns invisibly.
#'
#' @details
#' The function generates four views of the neurons:
#' 1. Main view
#' 2. Side view
#' 3. Front view (brain only)
#' 4. VNC (Ventral Nerve Cord) view
#'
#' Each view applies appropriate rotations and uses different neuropil meshes as backgrounds.
#' The two neurons are plotted in different colors (navy/turquoise for neuron1, red/darkred for neuron2) for easy comparison.
#'
#' @importFrom ggplot2 ggplot geom_point scale_color_gradient coord_fixed theme_void theme element_blank margin
#' @importFrom ggnewscale new_scale_colour new_scale_fill
#' @importFrom ggpubr ggarrange
#' @importFrom nat xyzmatrix
#'
#' @examples
#' \dontrun{
#' # Assuming neuron1 and neuron2 are valid neuron objects
#' banc_neuron_comparison_plot(neuron1, neuron2, filename = "neuron_comparison.png")
#' }
#'
#' @export
banc_neuron_comparison_plot <- function(neuron1,
                                        neuron2 = NULL,
                                        neuron3 = NULL,
                                        neuron1.info = NULL,
                                        neuron2.info = NULL,
                                        neuron3.info = NULL,
                                        banc_brain_neuropil = NULL,
                                        banc_vnc_neuropil = NULL,
                                        banc_neuropil = NULL,
                                        filename = NULL,
                                        width = 16,
                                        height = 16) {

  # Get 3D spatial points
  glist <- list()
  title.col <- "black"
  if(is.null(banc_brain_neuropil)) banc_brain_neuropil <- utils::data("banc_brain_neuropil.surf", envir = environment())
  if(is.null(banc_vnc_neuropil)) banc_vnc_neuropil <- utils::data("banc_vnc_neuropil.surf", envir = environment())
  if(is.null(banc_neuropil)) banc_neuropil <- utils::data("banc_neuropil.surf", envir = environment())
  for(view in names(banc_rotation_matrices)){

    # Choose mesh
    rotation_matrix <- banc_rotation_matrices[[view]]
    if(view=="front"){
      mesh <- banc_brain_neuropil
      decaptitate <- "brain"
    }else if(view=="vnc"){
      mesh <- banc_vnc_neuropil
      decaptitate <- "vnc"
    }else if (view=="main"){
      mesh <- banc_neuropil
      decaptitate <- "none"
      title.col <- "darkturquoise"
    }else if (view=="side"){
      mesh <- banc_neuropil
      decaptitate <- "none"
      title.col <- "firebrick"
    }

    # Get 3D points
    vertices <- nat::xyzmatrix(mesh)
    if(decaptitate=="brain"){
      neuron_pruned1 <- banc_decapitate(neuron1, invert = TRUE)
      neuron_pruned2 <- banc_decapitate(neuron2, invert = TRUE)
      neuron_pruned3 <- banc_decapitate(neuron3, invert = TRUE)
    }else if(decaptitate=="vnc"){
      neuron_pruned1 <- banc_decapitate(neuron1, invert = FALSE)
      neuron_pruned2 <- banc_decapitate(neuron2, invert = FALSE)
      neuron_pruned3 <- banc_decapitate(neuron3, invert = FALSE)
    }else{
      neuron_pruned1 <- neuron1
      neuron_pruned2 <- neuron2
      neuron_pruned3 <- neuron3
    }

    # Create the plot
    p <- ggplot2::ggplot() +
      geom_neuron.mesh3d(x = mesh, rotation_matrix = rotation_matrix, alpha = 0.05, low = "grey95", high = "grey50") +
      geom_neuron(x=neuron_pruned1[[1]], rotation_matrix = rotation_matrix, low = "turquoise", high = "navy", alpha = 0.5, linewidth = 0.3) +
      geom_neuron(x=neuron_pruned2[[1]], rotation_matrix = rotation_matrix, low = "red", high = "darkred", alpha = 0.5, linewidth = 0.3) +
      geom_neuron(x=neuron_pruned3[[1]], rotation_matrix = rotation_matrix, low = "chartreuse", high = "darkgreen", alpha = 0.5, linewidth = 0.3) +
      ggplot2::coord_fixed() +
      ggplot2::theme_void() +
      ggplot2::guides(fill="none",color="none") +
      ggplot2::theme(legend.position = "none",
            plot.title = ggplot2::element_text(hjust = 0, size = 8, face = "bold", colour = title.col),
            axis.title.x=ggplot2::element_blank(),
            axis.text.x=ggplot2::element_blank(),
            axis.ticks.x=ggplot2::element_blank(),
            axis.title.y=ggplot2::element_blank(),
            axis.text.y=ggplot2::element_blank(),
            axis.ticks.y=ggplot2::element_blank(),
            axis.line = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            plot.margin = ggplot2::margin(0, 0, 0, 0),
            panel.spacing = ggplot2::unit(0, "cm"),
            panel.border = ggplot2::element_blank(),
            panel.background = ggplot2::element_blank(), #gplot2::element_rect(fill = "grey95", color = NA),
            plot.background = ggplot2::element_blank()) + #gplot2::element_rect(fill = "grey95", color = NA))
      if(view=="main"){
        ggplot2::labs(title = neuron1.info)
      }else if (view=="side"){
        ggplot2::labs(title = neuron2.info)
      }else{
        ggplot2::labs(title = "")
      }

    # Add to list
    glist[[view]] <- p
  }

  # Arrange
  ga <- ggpubr::ggarrange(glist[["main"]], glist[["side"]], glist[["front"]], glist[["vnc"]],
            heights = c(1, 1, 1, 1), widths = c(1, 1, 1, 1),
            ncol = 2, nrow = 2) +
    ggplot2::theme(plot.margin = margin(0,0,0,0, "cm"))

  # # Annotate with some more information
  # if(!is.null(neuron1.info)&!is.null(neuron2.info)){
  #   ggpubr::annotate_figure(ga,
  #                           top = ggpubr::text_grob(neuron1.info, color = "blue", face = "bold", size = 12),
  #                           bottom = ggpubr::text_grob(neuron2.info, color = "coral", size = 12)
  #   )
  # }

  # Save
  if(!is.null(filename)){
    ggplot2::ggsave(filename, plot = ga, width = width, height = height, bg = "#fcfcfa")
    message("saved as:  ", filename)
    invisible()
  }else{
    plot(ga)
  }
}

#' Convert Neuron Objects to ggplot2-Compatible Data
#'
#' @description
#' This function converts 'neuron', 'mesh3d', or 'neuronlist' objects,
#' which represent 3D points linked by lines in space, into data frames
#' that describe paths compatible with ggplot2's geom_path, or geom_polygon
#' for mesh3d objects.
#'
#' @param x A 'neuron', 'neuronlist', or 'mesh3d' object to be converted.
#' @param rotation_matrix An optional 4x4 rotation matrix to apply to the neuron coordinates.
#' @param ... Additional arguments passed to methods.
#'
#' @return A data frame with columns X, Y, Z, and group, where each group
#' represents a continuous path in the neuron or a polygon in the mesh.
#'
#' @examples
#' \dontrun{
#' # Convert a neuron object
#' neuron_data <- ggplot2_neuron_path(my_neuron)
#'
#' # Convert a neuronlist object
#' neuronlist_data <- ggplot2_neuron_path(my_neuronlist)
#'
#' # Convert a mesh3d object
#' mesh_data <- ggplot2_neuron_path(my_mesh)
#' }
#'
#' @export
ggplot2_neuron_path <- function(x, rotation_matrix = NULL, ...) UseMethod('ggplot2_neuron_path')

#' @rdname ggplot2_neuron_path
#' @method ggplot2_neuron_path neuron
#' @export
ggplot2_neuron_path.neuron <- function(x, rotation_matrix = NULL, ...){
  x$d <- x$d[order(x$d$Parent),]
  x$d <- x$d[order(x$d$PointNo),]
  npoints <- as.data.frame(nat::xyzmatrix(x))
  if(!is.null(rotation_matrix)){
    npoints <- as.data.frame(t(rotation_matrix[,1:3] %*% t(nat::xyzmatrix(npoints))))
    npoints <- npoints[,-4]
    colnames(npoints) <- c("X","Y","Z")
  }
  ss <-nat::seglist(x)
  seglist <- ss[[1]]$SegList
  edges_df <- data.frame()
  for(s in 1:length(seglist)){
    g <- npoints[seglist[[s]],]
    g$group <- s
    edges_df <- rbind(edges_df,g)
  }
  edges_df
}

#' @rdname ggplot2_neuron_path
#' @method ggplot2_neuron_path neuronlist
#' @export
ggplot2_neuron_path.neuronlist <- function(x, rotation_matrix = NULL, ...){
  ll <- lapply(x, ggplot2_neuron_path, rotation_matrix = rotation_matrix, ...)
  max.group <- 0
  for(i in 1:length(ll)){
    ll[[i]]$group <- ll[[i]]$group + max.group
    ll[[i]]$id <- i
    max.group <- max(ll[[i]]$group, na.rm = TRUE)
  }
  do.call(rbind, ll)
}

#' @rdname ggplot2_neuron_path
#' @method ggplot2_neuron_path mesh3d
#' @export
ggplot2_neuron_path.mesh3d <- function(x, rotation_matrix = NULL, ...) {

  # Extract vertices
  vertices <- as.data.frame(t(x$vb[-4,]))
  if(!nrow(vertices)){
    warning("ggplot2_neuron_path.mesh3d given an invalid mesh3d object")
    return(data.frame(X=NA,Y=NA,Z=0,group=NA))
  }else if(nrow(vertices)<3){
    warning("ggplot2_neuron_path.mesh3d given an invalid mesh3d object")
    return(data.frame(X=NA,Y=NA,Z=0,group=NA))
  }
  colnames(vertices) <- c("X","Y","Z")

  # Apply rotation if specified
  if(!is.null(rotation_matrix)){
    vertices <- as.data.frame(t(rotation_matrix[,1:3] %*% t(nat::xyzmatrix(vertices))))
    vertices <- vertices[,-4]
    colnames(vertices) <- c("X","Y","Z")
  }

  # Extract faces
  faces.matrix <- t(x$it)
  faces <- as.vector(t(faces.matrix))

  # Create a data frame of triangles
  triangles <- data.frame(
    X = vertices$X[faces],
    Y = vertices$Y[faces],
    Z = vertices$Z[faces],
    group = rep(1:nrow(faces.matrix), each = 3)
  )
  triangles <- triangles %>%
    dplyr::arrange(dplyr::desc(.data$Z))

  # return: triangles
  triangles
}

#' Create ggplot2 Geom Layer for Neuron Visualization
#'
#' @description
#' This function creates a ggplot2 geom layer for visualizing neuron objects.
#' It supports 'neuron', 'neuronlist', and 'mesh3d' objects.
#'
#' @param x A 'neuron', 'neuronlist', or 'mesh3d' object to be visualized.
#' @param rotation_matrix An optional 4x4 rotation matrix to apply to the neuron coordinates.
#' @param low Color for the lowest Z values. The 'Z' axis is the axis perpendicular to the viewing plane.
#' @param high Color for the highest Z values.
#' @param stat The statistical transformation to use on the data for this layer.
#' @param position Position adjustment, either as a string, or the result of a call to a position adjustment function.
#' @param na.rm If FALSE, the default, missing values are removed with a warning. If TRUE, missing values are silently removed.
#' @param show.legend logical. Should this layer be included in the legends? NA, the default, includes if any aesthetics are mapped.
#' @param inherit.aes If FALSE, overrides the default aesthetics, rather than combining with them.
#' @param ... Other arguments passed on to layer().
#'
#' @return A list of ggplot2 geom layers for visualizing the neuron.
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' ggplot() + geom_neuron(my_neuron)
#' ggplot() + geom_neuron(my_neuronlist)
#' ggplot() + geom_neuron(my_mesh)
#' }
#'
#' @importFrom rlang .data
#' @export
geom_neuron <-function(x, rotation_matrix = NULL, low = "turquoise", high = "navy",
                       stat = "identity", position = "identity", na.rm = FALSE, show.legend = NA,
                       inherit.aes = FALSE, ...) UseMethod('geom_neuron')

#' @rdname geom_neuron
#' @method geom_neuron neuron
#' @export
geom_neuron.neuron <- function(x = NULL, rotation_matrix = NULL, low = "turquoise", high = "navy",
                               stat = "identity", position = "identity", na.rm = FALSE, show.legend = NA,
                               inherit.aes = FALSE, ...) {
  soma <- catmaid::soma(x)
  if(!is.null(rotation_matrix)){
    soma <- as.data.frame(t(rotation_matrix[,1:3] %*% t(nat::xyzmatrix(soma))))
    soma <- soma[,-4]
    colnames(soma) <- c("X","Y","Z")
  }
  x <- ggplot2_neuron_path.neuron(x, rotation_matrix = rotation_matrix)
  list(
    ggplot2::geom_path(mapping = ggplot2::aes(x = .data$X, y = .data$Y, color = .data$Z, group = .data$group),
                       data = x,
              stat = stat, position = position, na.rm = na.rm,
              show.legend = show.legend, inherit.aes = inherit.aes, ...),
    ggplot2::geom_point(mapping = ggplot2::aes(x = .data$X, y = .data$Y), data = soma,
                         color = high, alpha = 0.5, size = 3),
    ggplot2::scale_color_gradient(low = low, high = high),
    ggnewscale::new_scale_colour()
  )
}

#' @rdname geom_neuron
#' @method geom_neuron neuronlist
#' @export
geom_neuron.neuronlist <- function(x = NULL, rotation_matrix = NULL, low = "turquoise", high = "navy",
                               stat = "identity", position = "identity", na.rm = FALSE, show.legend = NA,
                               inherit.aes = FALSE, ...) {
  x <- ggplot2_neuron_path.neuronlist(x, rotation_matrix = rotation_matrix)
  list(
    ggplot2::geom_path(mapping = ggplot2::aes(x = .data$X, y = .data$Y, color = .data$id, group = .data$group), data = x,
                       stat = stat, position = position, na.rm = na.rm,
                       show.legend = show.legend, inherit.aes = inherit.aes, ...)
  )
}

#' @rdname geom_neuron
#' @method geom_neuron mesh3d
#' @export
geom_neuron.mesh3d <- function(x = NULL, rotation_matrix = NULL, low = "turquoise", high = "navy",
                                   stat = "identity", position = "identity", na.rm = FALSE, show.legend = NA,
                                   inherit.aes = FALSE, ...) {
  x <- ggplot2_neuron_path.mesh3d(x, rotation_matrix = rotation_matrix)
  list(
    ggplot2::geom_polygon(data = x, mapping = ggplot2::aes(x = .data$X, y = .data$Y, fill = .data$Z, group = .data$group),
                    color = NA,
                    stat = stat, position = position, na.rm = na.rm,
                    show.legend = show.legend, inherit.aes = inherit.aes, ...),
    ggplot2::scale_fill_gradient(low = low, high = high),
    ggnewscale::new_scale_fill()
  )
}

#' @rdname geom_neuron
#' @method geom_neuron hxsurf
#' @export
geom_neuron.hxsurf <- function(x = NULL, rotation_matrix = NULL, low = "turquoise", high = "navy",
                               stat = "identity", position = "identity", na.rm = FALSE, show.legend = NA,
                               inherit.aes = FALSE, ...) {
  x <- rgl::as.mesh3d(x)
  geom_neuron.mesh3d(x=x, rotation_matrix=rotation_matrix, low=low, high=high,
                     stat=stat, position=position, na.rm=na.rm, show.legend=show.legend, inherit.aes=inherit.aes,
                     ...)
}

#' @rdname geom_neuron
#' @method geom_neuron NULL
#' @export
geom_neuron.NULL <- function(x = NULL, rotation_matrix = NULL, low = "turquoise", high = "navy",
                               stat = "identity", position = "identity", na.rm = FALSE, show.legend = NA,
                               inherit.aes = FALSE, ...) {
  list(
    ggplot2::geom_polygon(...)
  )
}










