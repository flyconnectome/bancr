#' Compare two neurons from the BANC connectome dataset
#'
#' @description
#' This function creates a visual comparison of two neurons from the BANC (Bilateral Antennal lobe Neuron Connectome)
#' dataset. It generates four different anatomical views to facilitate comparison.
#'
#' @param neuron1 A neuronlist object representing the first neuron to be compared.
#' @param neuron2 A neuronlist object representing the second neuron to be compared.
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
                                        banc_brain_neuropil = banc_brain_neuropil.surf,
                                        banc_vnc_neuropil = banc_vnc_neuropil.surf,
                                        banc_neuropil = banc_neuropil.surf,
                                        filename = NULL,
                                        width = 16,
                                        height = 16) {

  # Get 3D spatial points
  check_package_available('ggplot2')
  check_package_available('ggpubr')
  glist <- list()
  title.col <- "black"
  if(!"mesh3d"%in%class(banc_brain_neuropil)){
    banc_brain_neuropil <- rgl::as.mesh3d(banc_brain_neuropil)
  }
  if(!"mesh3d"%in%class(banc_vnc_neuropil)){
    banc_vnc_neuropil <- rgl::as.mesh3d(banc_vnc_neuropil)
  }
  if(!"mesh3d"%in%class(banc_neuropil)){
    banc_neuropil <- rgl::as.mesh3d(banc_neuropil)
  }
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
      geom_neuron.mesh3d(x = mesh, rotation_matrix = rotation_matrix, alpha = 0.05) +
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

#' Convert neuron objects to ggplot2-compatible data
#'
#' @description
#' This function converts 'neuron', 'mesh3d' or 'neuronlist' objects, which represent 3D
#' points linked by lines in space, into data frames that describe paths
#' compatible with ggplot2's geom_path, or geom_polygon for mesh3d objects.
#'
#' @param x A 'neuron' or 'neuronlist' object to be converted.
#' @param rotation_matrix An optional 4x4 rotation matrix to apply to the neuron coordinates.
#' @param ... Additional arguments passed to methods.
#'
#' @return A data frame with columns X, Y, Z, and group, where each group
#' represents a continuous path in the neuron.
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
    dplyr::arrange(dplyr::desc(Z))

  # return: triangles
  triangles
}


#' @rdname ggplot2_neuron_path
#' @export
geom_neuron <-function(x, rotation_matrix = NULL, low = "turquoise", high = "navy",
                       stat = "identity", position = "identity", na.rm = FALSE, show.legend = NA,
                       inherit.aes = FALSE, ...) UseMethod('geom_neuron')

#' @rdname ggplot2_neuron_path
#' @method geom_neuron neuron
#' @export
geom_neuron.neuron <- function(x = NULL, rotation_matrix = NULL, low = "turquoise", high = "navy",
                               stat = "identity", position = "identity", na.rm = FALSE, show.legend = NA,
                               inherit.aes = FALSE, ...) {

  check_package_available('ggnewscale')
  check_package_available('catmaid')
  soma <- catmaid::soma(x)
  if(!is.null(rotation_matrix)){
    soma <- as.data.frame(t(rotation_matrix[,1:3] %*% t(nat::xyzmatrix(soma))))
    soma <- soma[,-4]
    colnames(soma) <- c("X","Y","Z")
  }
  x <- ggplot2_neuron_path.neuron(x, rotation_matrix = rotation_matrix)
  list(
    ggplot2::geom_path(mapping = ggplot2::aes(x = X, y = Y, color = Z, group = group), data = x,
              stat = stat, position = position, na.rm = na.rm,
              show.legend = show.legend, inherit.aes = inherit.aes, ...),
    ggplot2::geom_point(mapping = ggplot2::aes(x = X, y = Y), data = soma,
                         color = high, alpha = 0.5, size = 3),
    ggplot2::scale_color_gradient(low = low, high = high),
    ggnewscale::new_scale_colour()
  )
}

#' @rdname ggplot2_neuron_path
#' @method geom_neuron neuronlist
#' @export
geom_neuron.neuronlist <- function(x = NULL, rotation_matrix = NULL, low = "turquoise", high = "navy",
                               stat = "identity", position = "identity", na.rm = FALSE, show.legend = NA,
                               inherit.aes = FALSE, ...) {
  x <- ggplot2_neuron_path.neuronlist(x, rotation_matrix = rotation_matrix)
  list(
    ggplot2::geom_path(mapping = ggplot2::aes(x = X, y = Y, color = id, group = group), data = x,
                       stat = stat, position = position, na.rm = na.rm,
                       show.legend = show.legend, inherit.aes = inherit.aes, ...)
  )
}

#' @rdname ggplot2_neuron_path
#' @method geom_neuron neuronlist
#' @export
geom_neuron.mesh3d <- function(x = NULL, rotation_matrix = NULL, low = "grey90", high = "grey50",
                                   stat = "identity", position = "identity", na.rm = FALSE, show.legend = NA,
                                   inherit.aes = FALSE, ...) {
  check_package_available('ggnewscale')
  x <- ggplot2_neuron_path.mesh3d(x, rotation_matrix = rotation_matrix)
  list(
    ggplot2::geom_polygon(data = x, mapping = ggplot2::aes(x = X, y = Y, fill = Z, group = group),
                    color = NA,
                    stat = stat, position = position, na.rm = na.rm,
                    show.legend = show.legend, inherit.aes = inherit.aes, ...),
    ggplot2::scale_fill_gradient(low = low, high = high),
    ggnewscale::new_scale_fill()
  )
}

#' @rdname ggplot2_neuron_path
#' @method geom_neuron neuronlist
#' @export
geom_neuron.NULL <- function(x = NULL, rotation_matrix = NULL, low = "grey90", high = "grey50",
                               stat = "identity", position = "identity", na.rm = FALSE, show.legend = NA,
                               inherit.aes = FALSE, ...) {
  list(
    ggplot2::geom_polygon(...)
  )
}










