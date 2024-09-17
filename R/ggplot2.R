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
#' @param volume A mesh3d  or hxsurf object in BANC space that you wish you co-plot.
#' @param region Character, whether to plot the brain area, VNC area or both (default).
#' @param banc_brain_neuropil A mesh object representing the brain neuropil. Default is banc_brain_neuropil.surf.
#' @param banc_vnc_neuropil A mesh object representing the ventral nerve cord (VNC) neuropil. Default is banc_vnc_neuropil.surf.
#' @param banc_neuropil A mesh object representing the entire neuropil. Default is banc_neuropil.surf.
#' @param filename Optional. If provided, the plot will be saved to this file. If NULL (default), the plot will be displayed but not saved.
#' @param width Numeric. The width of the output plot in inches. Default is 16.
#' @param height Numeric. The height of the output plot in inches. Default is 16.
#' @param alpha Vector of alpha values for neurons 1, 2 and 3. If a single value, applied to all.
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
#' banc.meta <- banctable_query()
#'
#' # Get some neurons to plot
#' banc.meta.dnao1 <- subset(banc.meta, cell_type=="DNa01")
#' dna01 <- banc_read_neuron_meshes(banc.meta.dnao1$root_id)
#' banc.meta.dnao2 <- subset(banc.meta, cell_type=="DNa02")
#' dna02 <- banc_read_neuron_meshes(banc.meta.dnao2$root_id)
#'
#' # Simplify neurons to make them easier to plot
#' dna01 <- nat::nlapply(dna01,Rvcg::vcgQEdecim,percent = 0.1)
#' dna02 <- nat::nlapply(dna02,Rvcg::vcgQEdecim,percent = 0.1)
#'
#' # Make plot!
#' banc_neuron_comparison_plot(dna01, dna02, neuron1.info = "DNa01",
#' neuron2.info = "DNa02", filename = "neuron_comparison.png")
#' }
#'
#' @export
banc_neuron_comparison_plot <- function(neuron1 = NULL,
                                        neuron2 = NULL,
                                        neuron3 = NULL,
                                        neuron1.info = NULL,
                                        neuron2.info = NULL,
                                        neuron3.info = NULL,
                                        volume = NULL,
                                        region = c("both","brain","vnc"),
                                        banc_brain_neuropil = NULL,
                                        banc_vnc_neuropil = NULL,
                                        banc_neuropil = NULL,
                                        alpha = 0.5,
                                        filename = NULL,
                                        width = 16,
                                        height = 16) {

  # Get 3D spatial points
  region <- match.arg(region)
  check_package_available('ggplot2')
  check_package_available('ggpubr')
  glist <- list()
  title.col <- "black"
  if(is.null(banc_brain_neuropil)) {banc_brain_neuropil <- bancr::banc_brain_neuropil.surf}
  if(is.null(banc_vnc_neuropil)) {banc_vnc_neuropil <- bancr::banc_vnc_neuropil.surf}
  if(is.null(banc_neuropil)) {banc_neuropil <- bancr::banc_neuropil.surf}
  views <- names(banc_rotation_matrices)
  if(region=="brain"){
    views <- c("front", "brain_side")
    banc_neuropil <- NULL
    banc_vnc_neuropil <- NULL
  }else if(region == "vnc"){
    views <- c("vnc", "vnc_side")
    banc_neuropil <- NULL
    banc_brain_neuropil <- NULL
  }else{
    views <- c("main", "side", "front", "vnc")
  }
  for(view in views){

    # Choose mesh
    rotation_matrix <- banc_rotation_matrices[[view]]
    if(view%in%c("front")){
      mesh <- banc_brain_neuropil
      decaptitate <- "brain"
      title.col <- "blue"
    }else if(view%in%c("brain_side")){
      mesh <- banc_brain_neuropil
      decaptitate <- "brain"
      title.col <- "firebrick"
    }else if(view%in%c("vnc","vnc_side")){
      mesh <- banc_vnc_neuropil
      decaptitate <- "vnc"
      title.col <- "firebrick"
    }else if (view=="main"){
      mesh <- banc_neuropil
      decaptitate <- "none"
      title.col <- "blue"
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

    # Work out colours
    if(length(neuron1)==1){
      cols1 <- c("blue", "navy")
    }else{
      cols1 <- grDevices::colorRampPalette(c("#0000CD", "#4169E1"))(length(neuron1))
    }
    if(length(neuron2)==1){
      cols2 <- c("darkred", "red")
    }else{
      cols2 <- grDevices::colorRampPalette(c("#8f0723","#DC143C", "#FF4500"))(length(neuron2))
    }
    if(length(neuron3)==1){
      cols3 <- c("darkgreen", "green")
    }else{
      cols3 <- grDevices::colorRampPalette(c("#076b3e", "#32c080"))(length(neuron3))
    }
    if(length(alpha)<3){
      alpha<-rep(alpha[1],3)
    }

    # Create the plot
    p <- ggplot2::ggplot() +
      geom_neuron(x = mesh, rotation_matrix = rotation_matrix, alpha = 0.05, cols = c("grey90", "grey50")) +
      geom_neuron(x = volume, rotation_matrix = rotation_matrix, alpha = 0.05, cols = c("grey30")) +
      geom_neuron(x=neuron_pruned1, rotation_matrix = rotation_matrix, cols = cols1, alpha = alpha[1], linewidth = 0.3) +
      geom_neuron(x=neuron_pruned2, rotation_matrix = rotation_matrix, cols = cols2, alpha = alpha[2], linewidth = 0.3) +
      geom_neuron(x=neuron_pruned3, rotation_matrix = rotation_matrix, cols = cols3, alpha = alpha[3], linewidth = 0.3) +
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
      if(view%in%c("front")){
        ggplot2::labs(title = neuron1.info)
      }else if (view%in%c("vnc","brain_side","vnc_side")){
        ggplot2::labs(title = neuron2.info)
      }else{
        ggplot2::labs(title = "")
      }

    # Add to list
    glist[[view]] <- p
  }

  # Arrange
  if(region=="both"){
    ga <- ggpubr::ggarrange(glist[["front"]], glist[["vnc"]], glist[["main"]], glist[["side"]],
                            heights = c(1, 1, 1, 1), widths = c(1, 1, 1, 1),
                            ncol = 2, nrow = 2) +
      ggplot2::theme(plot.margin = ggplot2::margin(0,0,0,0, "cm"))
  }else{
    ga <- ggpubr::ggarrange(glist[[1]], glist[[2]],
                            heights = c(1, 1), widths = c(1, 1),
                            ncol = 1, nrow = 2) +
      ggplot2::theme(plot.margin = ggplot2::margin(0,0,0,0, "cm"))
  }

  # # Annotate with some more information
  # if(!is.null(neuron1.info)&!is.null(neuron2.info)){
  #   ggpubr::annotate_figure(ga,
  #                           top = ggpubr::text_grob(neuron1.info, color = "blue", face = "bold", size = 12),
  #                           bottom = ggpubr::text_grob(neuron2.info, color = "coral", size = 12)
  #   )
  # }

  # Save
  if(!is.null(filename)){
    ggplot2::ggsave(filename, plot = ga, width = width, height = height,
                    bg = "#fcfcfa", dpi = 300, limitsize = FALSE)
    message("saved as:  ", filename)
    invisible()
  }else{
    plot(ga)
  }
}

#' Plot a neuron in the BANC connectomic dataset using ggplot2
#'
#' This function visualizes a neuron or neuron-related object from the BANC connectomic dataset using ggplot2.
#' The only thing specific to the BANC data set is are the prreset 'view' angles.
#'
#' @param x A 'neuron', 'neuronlist', 'mesh3d', or 'hxsurf' object to be visualized.
#' @param volume A brain/neuropil volume to be plotted in grey, for context.
#'   Defaults to NULL, no volume plotted.
#' @param info Optional. A string to be used as the plot title.
#' @param view A character string specifying the view orientation.
#'   Options are "main", "side", "front", "vnc", "vnc_side", "brain_side".
#' @param cols1 A vector of two colors for the lowest Z values. Default is c("turquoise", "navy").
#' @param cols2 A vector of two colors for the highest Z values. Default is c("grey75", "grey50").
#' @param alpha Transparency of the neuron visualization. Default is 0.5.
#' @param title.col Color of the plot title. Default is "darkgrey".
#' @param ... Additional arguments passed to geom_neuron().
#'
#' @return A ggplot object representing the visualized neuron.
#'
#' @details This function is a wrapper around the ggneuron function, specifically tailored for the BANC dataset.
#'   It applies a rotation matrix based on the specified view and uses predefined color schemes.
#'
#' @export
banc_ggneuron <-function(x,
                         volume = NULL,
                         info = NULL,
                         view = c("main", "side", "front", "vnc", "vnc_side", "brain_side"),
                         cols1 = c("turquoise","navy"),
                         cols2 =  c("grey75", "grey50"),
                         alpha = 0.5,
                         title.col = "darkgrey",
                         ...){
  view <- match.arg(view)
  rotation_matrix <- banc_rotation_matrices[[view]]
  ggneuron(x,
           volume = volume,
           info = info,
           rotation_matrix = rotation_matrix,
           cols1 = cols1,
           cols2 =  cols2,
           alpha = alpha,
           title.col = title.col,
           ...)
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

#' @rdname ggplot2_neuron_path
#' @method ggplot2_neuron_path NULL
#' @export
ggplot2_neuron_path.NULL <- function(x, rotation_matrix = NULL, ...) {
  NULL
}

#' Create ggplot2 Geom Layer for Neuron Visualization
#'
#' @description
#' This function creates a ggplot2 geom layer for visualizing neuron objects.
#' It supports 'neuron', 'neuronlist', and 'mesh3d' objects.
#'
#' @param x A 'neuron', 'neuronlist', or 'mesh3d' object to be visualized.
#' @param rotation_matrix An optional 4x4 rotation matrix to apply to the neuron coordinates.
#' @param root Numeric, if >0 and x is or contains `neuron` objects,
#' then the root node is plotted as a dot of size `root`. If `FALSE` or `0` no root node is plotted.
#' @param cols The color to plot the neurons in. If \code{length(cols)==length(x)} each neuron will be coloured
#' by its index in `x` applied to `cols`.
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
#' ggplot() + geom_neuron(nat::Cell07PNs[[1]], root = 10)
#' ggplot() + geom_neuron(nat::Cell07PNs)
#' ggplot() + geom_neuron(nat::MBL.surf)
#' }
#'
#' @importFrom rlang .data
#' @export
geom_neuron <-function(x, rotation_matrix = NULL, root = 3, cols = c("navy", "turquoise"),
                       stat = "identity", position = "identity", na.rm = FALSE, show.legend = NA,
                       inherit.aes = FALSE, ...) UseMethod('geom_neuron')

#' @rdname geom_neuron
#' @method geom_neuron neuron
#' @export
geom_neuron.neuron <- function(x = NULL, rotation_matrix = NULL, root = 3, cols = c("navy", "turquoise"),
                               stat = "identity", position = "identity", na.rm = FALSE, show.legend = NA,
                               inherit.aes = FALSE, ...) {
  check_package_available('ggnewscale')
  check_package_available('catmaid')
  if(root){
    x$tags$soma <- nat::rootpoints(x)
  }
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
                         color = cols[1], alpha = 0.5, size = root),
    ggplot2::scale_color_gradient(low = cols[1], high = cols[length(cols)]),
    ggnewscale::new_scale_colour()
  )
}

#' @rdname geom_neuron
#' @method geom_neuron neuronlist
#' @export
geom_neuron.neuronlist <- function(x = NULL, rotation_matrix = NULL, root = 3, cols = c("navy", "turquoise"),
                               stat = "identity", position = "identity", na.rm = FALSE, show.legend = NA,
                               inherit.aes = FALSE, ...) {
  glist <- list()
  if(length(x)!=1){
    if(cols[1]=="rainbow"){
      cols <- grDevices::rainbow(length(x))
    }else if(length(cols)==1){
      cols <- rep(cols, length(x))
    }else{
      cols <- grDevices::colorRampPalette(c(cols[1], cols[length(cols)]))(length(x))
    }
    for(i in 1:length(x)){
      glist[[i]] <- geom_neuron(x = x[[i]], rotation_matrix = rotation_matrix, cols = cols[i],
                                stat = stat, position = position, na.rm = na.rm, show.legend = show.legend,
                                inherit.aes = FALSE, ...)
    }
  }else{
    if(cols[1]=="rainbow"){
      cols <-c("purple","pink")
    }
    for(i in 1:length(x)){
      glist[[i]] <- geom_neuron(x = x[[i]], rotation_matrix = rotation_matrix, cols = cols,
                                stat = stat, position = position, na.rm = na.rm, show.legend = show.legend,
                                inherit.aes = FALSE, ...)
    }
  }
  glist
}

#' @rdname geom_neuron
#' @method geom_neuron mesh3d
#' @export
geom_neuron.mesh3d <- function(x = NULL, rotation_matrix = NULL, root = 3, cols = c("navy", "turquoise"),
                                   stat = "identity", position = "identity", na.rm = FALSE, show.legend = NA,
                                   inherit.aes = FALSE, ...) {
  check_package_available('ggnewscale')
  x <- ggplot2_neuron_path.mesh3d(x, rotation_matrix = rotation_matrix)
  list(
    ggplot2::geom_polygon(data = x, mapping = ggplot2::aes(x = .data$X, y = .data$Y, fill = .data$Z, group = .data$group),
                    color = NA,
                    stat = stat, position = position, na.rm = na.rm,
                    show.legend = show.legend, inherit.aes = inherit.aes, ...),
    ggplot2::scale_fill_gradient(low = cols[1], high = cols[length(cols)]),
    ggnewscale::new_scale_fill()
  )
}

#' @rdname geom_neuron
#' @method geom_neuron hxsurf
#' @export
geom_neuron.hxsurf <- function(x = NULL, rotation_matrix = NULL, root = 3, cols = c("navy", "turquoise"),
                               stat = "identity", position = "identity", na.rm = FALSE, show.legend = NA,
                               inherit.aes = FALSE, ...) {
  x <- rgl::as.mesh3d(x)
  geom_neuron.mesh3d(x=x, rotation_matrix=rotation_matrix, cols=cols,
                     stat=stat, position=position, na.rm=na.rm, show.legend=show.legend, inherit.aes=inherit.aes,
                     ...)
}

#' @rdname geom_neuron
#' @method geom_neuron NULL
#' @export
geom_neuron.NULL <- function(x = NULL, rotation_matrix = NULL, root = 3, cols = c("navy", "turquoise"),
                               stat = "identity", position = "identity", na.rm = FALSE, show.legend = NA,
                               inherit.aes = FALSE, ...) {
  list(
    ggplot2::geom_polygon(...)
  )
}

#' @rdname geom_neuron
#' @method geom_neuron list
#' @export
geom_neuron.list <- function(x = NULL, rotation_matrix = NULL, root = 3, cols = c("navy", "turquoise"),
                             stat = "identity", position = "identity", na.rm = FALSE, show.legend = NA,
                             inherit.aes = FALSE, ...) {

  if(length(x)){
    geom_neuron.neuronlist(x=x, rotation_matrix=rotation_matrix, cols=cols,
                       stat=stat, position=position, na.rm=na.rm, show.legend=show.legend, inherit.aes=inherit.aes,
                       ...)
  }else{
    geom_neuron.NULL(x = x, ...)
  }
}

#' @rdname geom_neuron
#' @method geom_neuron matrix
#' @export
geom_neuron.matrix <- function(x = NULL, rotation_matrix = NULL, root = 3, cols = c("navy", "turquoise"),
                             stat = "identity", position = "identity", na.rm = FALSE, show.legend = NA,
                             inherit.aes = FALSE, ...) {
  x<-as.data.frame(nat::xyzmatrix(x))
  if(!is.null(rotation_matrix)){
    x <- as.data.frame(t(rotation_matrix[,1:3] %*% t(nat::xyzmatrix(x))))
    x <- x[,-4]
    colnames(x) <- c("X","Y","Z")
  }
  list(
    ggplot2::geom_point(data = x, mapping = ggplot2::aes(x = .data$X, y = .data$Y, color = .data$Z),
                        size = root,  ...),
    ggplot2::scale_color_gradient(low = cols[1], high = cols[length(cols)]),
    ggnewscale::new_scale_colour()
  )
}

#' @rdname geom_neuron
#' @method geom_neuron data.frame
#' @export
geom_neuron.data.frame <- function(x = NULL, rotation_matrix = NULL, root = 3, cols = c("navy", "turquoise"),
                               stat = "identity", position = "identity", na.rm = FALSE, show.legend = NA,
                               inherit.aes = FALSE, ...) {
  geom_neuron.matrix(x, rotation_matrix = rotation_matrix, root = root,
                     low = cols[1], high = cols[length(cols)],
                     stat = stat, position = position,
                     na.rm = FALSE, show.legend = NA,
                     inherit.aes = FALSE,
                     ...)
}

#' @rdname geom_neuron
#' @method geom_neuron dotprops
#' @export
geom_neuron.dotprops <- function(x = NULL, rotation_matrix = NULL, root = 3, cols = c("navy", "turquoise"),
                                   stat = "identity", position = "identity", na.rm = FALSE, show.legend = NA,
                                   inherit.aes = FALSE, ...) {
  x<-as.data.frame(nat::xyzmatrix(x))
  geom_neuron.data.frame(x, rotation_matrix = rotation_matrix, root = root,
                     low = cols[1], high = cols[length(cols)],
                     stat = stat, position = position,
                     na.rm = FALSE, show.legend = NA,
                     inherit.aes = FALSE,
                     ...)
}

#' @rdname geom_neuron
#' @method geom_neuron synapticneuron
#' @export
geom_neuron.synapticneuron <- function(x = NULL,
                                       rotation_matrix = NULL,
                                       root = 3,
                                       cols = c("navy", "turquoise"),
                                       stat = "identity", position = "identity",
                                       na.rm = FALSE,
                                       show.legend = NA,
                                       inherit.aes = FALSE,
                                       ...) {
  geomneuron<-geom_neuron.neuron(x = x,
                                 rotation_matrix = rotation_matrix,
                                 root = root,
                                 cols = cols,
                                 stat = stat,
                                 position = position,
                                 na.rm = na.rm,
                                 show.legend = show.legend,
                                 inherit.aes = inherit.aes,
                                 ...)
  if(!is.null(x$connectors)){
    syns.in <- nat::xyzmatrix(subset(x$connectors, x$connectors$prepost==1))
    syns.out <- nat::xyzmatrix(subset(x$connectors, x$connectors$prepost==0))
    if(!is.null(rotation_matrix)){
      syns.in <- as.data.frame(t(rotation_matrix[,1:3] %*% t(syns.in)))
      syns.in <- syns.in[,-4]
      colnames(syns.in) <- c("X","Y","Z")
      syns.out <- as.data.frame(t(rotation_matrix[,1:3] %*% t(syns.out)))
      syns.out <- syns.out[,-4]
      colnames(syns.out) <- c("X","Y","Z")
    }
    glist <- list(
      ggplot2::geom_point(data = syns.in,
                          mapping = ggplot2::aes(x = .data$X,
                                                 y = .data$Y),
                          color = "#132157",
                          size = root/100,
                          alpha = 0.25),
      ggplot2::geom_point(data = syns.out,
                          mapping = ggplot2::aes(x = .data$X,
                                                 y = .data$Y),
                          color = "#D72000",
                          size = root/100,
                          alpha = 0.25)
    )
    c(geomneuron,glist)
  }else{
    geomneuron
  }
}

#' @rdname geom_neuron
#' @method geom_neuron splitneuron
#' @export
geom_neuron.splitneuron <- function(x = NULL,
                                   rotation_matrix = NULL,
                                   root = 3,
                                   cols = c("navy", "turquoise"),
                                   stat = "identity", position = "identity",
                                   na.rm = FALSE,
                                   show.legend = NA,
                                   inherit.aes = FALSE,
                                   ...) {

  # Get parts
  if(root){
    x$tags$soma <- nat::rootpoints(x)
  }
  soma <- nat::xyzmatrix(x)[rootpoints(x),] # catmaid::soma(x)
  soma = t(as.data.frame(soma))
  if(!is.null(rotation_matrix)){
    soma <- as.data.frame(t(rotation_matrix[,1:3] %*% t(nat::xyzmatrix(soma))))
    soma <- soma[,-4]
    colnames(soma) <- c("X","Y","Z")
  }
  rownames(x$d) <- 1:nrow(x$d)
  dendrites.v = subset(rownames(x$d), x$d$Label == 3)
  axon.v = subset(rownames(x$d), x$d$Label == 2)
  p.d.v = subset(rownames(x$d), x$d$Label == 4)
  p.n.v = subset(rownames(x$d), x$d$Label == 7)
  null.v = subset(rownames(x$d), x$d$Label ==  0 | is.na(x$d$Label))

  # Get cable
  dendrites <- tryCatch(nat::prune_vertices(x,
                                           verticestoprune = as.integer(c(axon.v,p.d.v, p.n.v, null.v))),
                       error = function(e) NULL)
  axon <- tryCatch(nat::prune_vertices(x,
                                      verticestoprune = as.integer(c(dendrites.v, p.d.v, p.n.v, null.v))),
                  error = function(e) NULL)
  p.d <- tryCatch(nat::prune_vertices(x, verticestoprune = as.integer(c(axon.v, dendrites.v, p.n.v, null.v))),
                 error = function(e) NULL)
  p.n <- tryCatch(nat::prune_vertices(x, verticestoprune = as.integer(c(axon.v, dendrites.v, p.d.v, null.v))),
                 error = function(e) NULL)
  nulls <- tryCatch(nat::prune_vertices(x, verticestoprune = as.integer(c(axon.v, dendrites.v, p.d.v, p.n.v))),
                 error = function(e) NULL)

  # Stitch subtree
  dendrites <- nat::stitch_neurons_mst(dendrites)
  axon <- nat::stitch_neurons_mst(axon)
  p.d <- nat::stitch_neurons_mst(p.d)
  p.n <- nat::stitch_neurons_mst(p.n)

  # Make into a multi-segmen neuroblist
  nulls <- nat::nlapply(1:length(nulls$SubTrees), function(subt) tryCatch(nat::prune_vertices(nulls,
                                                                                verticestoprune = unlist(nulls$SubTrees[[subt]]),
                                                                                invert = TRUE),
                                                                          error = function(e) NULL),
                                                                          .progress = FALSE)
  nulls <- nulls[unlist(lapply(nulls, length))>0]

  # Make ggplot2 objects
  g.dendrites <- ggplot2_neuron_path(dendrites, rotation_matrix = rotation_matrix)
  g.axon <- ggplot2_neuron_path(axon, rotation_matrix = rotation_matrix)
  g.p.d <- ggplot2_neuron_path(p.d, rotation_matrix = rotation_matrix)
  g.p.n <- ggplot2_neuron_path(p.n, rotation_matrix = rotation_matrix)
  g.nulls <- ggplot2_neuron_path(nulls, rotation_matrix = rotation_matrix)

  # Make geom objects
  glist <- list(
    if(length(g.dendrites)){
      ggplot2::geom_path(mapping = ggplot2::aes(x = .data$X, y = .data$Y, group = .data$group),
                         data = g.dendrites, col = "#54BCD1", na.rm = TRUE,
                         stat = stat, position = position,
                         show.legend = show.legend, inherit.aes = inherit.aes, alpha = 1)
    }else{
      NULL
    },
    if(length(g.axon)){
      ggplot2::geom_path(mapping = ggplot2::aes(x = .data$X, y = .data$Y, group = .data$group),
                         data = g.axon, col = "#EF7C12", na.rm = TRUE,
                         stat = stat, position = position,
                         show.legend = show.legend, inherit.aes = inherit.aes,  alpha = 1)
    }else{
      NULL
    },
    if(length(g.p.d)){
      ggplot2::geom_path(mapping = ggplot2::aes(x = .data$X, y = .data$Y, group = .data$group),
                         data = g.p.d, col = "#8FDA04", na.rm = TRUE,
                         stat = stat, position = position,
                         show.legend = show.legend, inherit.aes = inherit.aes,  alpha = 1)
    }else{
      NULL
    },
    if(length(g.p.n)){
      ggplot2::geom_path(mapping = ggplot2::aes(x = .data$X, y = .data$Y, group = .data$group),
                         data = g.p.n, col = "#C70E7B", na.rm = TRUE,
                         stat = stat, position = position,
                         show.legend = show.legend, inherit.aes = inherit.aes,  alpha = 1)
    }else{
      NULL
    },
    if(length(g.nulls)){
      ggplot2::geom_path(mapping = ggplot2::aes(x = .data$X, y = .data$Y, group = .data$group),
                         data = g.nulls, col = "#B3B3B3", na.rm = TRUE,
                         stat = stat, position = position,
                         show.legend = show.legend, inherit.aes = inherit.aes,  alpha = 1)
    }else{
      NULL
    },
    ggplot2::geom_point(mapping = ggplot2::aes(x = .data$X, y = .data$Y),
                        data = soma, col = "black",
                        color = cols[1], alpha = 0.75, size = root)
    )

  # And synapses?
  if(!is.null(x$connectors)){
    syns.in <- nat::xyzmatrix(subset(x$connectors, x$connectors$prepost==1))
    syns.out <- nat::xyzmatrix(subset(x$connectors, x$connectors$prepost==0))
    if(!is.null(rotation_matrix)){
      syns.in <- as.data.frame(t(rotation_matrix[,1:3] %*% t(syns.in)))
      syns.in <- syns.in[,-4]
      colnames(syns.in) <- c("X","Y","Z")
      syns.out <- as.data.frame(t(rotation_matrix[,1:3] %*% t(syns.out)))
      syns.out <- syns.out[,-4]
      colnames(syns.out) <- c("X","Y","Z")
    }
    syn.glist <- list(
      ggplot2::geom_point(data = syns.in,
                          mapping = ggplot2::aes(x = .data$X,
                                                 y = .data$Y),
                          color = "#132157",
                          size = root/100,
                          alpha = 0.25),
      ggplot2::geom_point(data = syns.out,
                          mapping = ggplot2::aes(x = .data$X,
                                                 y = .data$Y),
                          color = "#D72000",
                          size = root/100,
                          alpha = 0.25)
    )
    c(glist,syn.glist)
  }else{
    glist
  }
}



#' Create a ggplot2 Visualisation of Neuron Objects
#'
#' @description
#' This function creates a complete ggplot2 visualization for neuron objects,
#' including 'neuron', 'neuronlist', 'mesh3d', and 'hxsurf' objects. It sets up
#' a minimal theme and applies consistent styling to the plot.
#'
#' @param x A 'neuron', 'neuronlist', 'mesh3d', or 'hxsurf' object to be visualized.
#' @param volume a brain/neuropil volume to be plotted in grey, for context.
#' Defaults to NULL, no volume plotted.
#' @param info Optional. A string to be used as the plot title.
#' @param rotation_matrix An optional 4x4 rotation matrix to apply to the neuron coordinates.
#' @param cols1 Color for the lowest Z values. Default is "turquoise".
#' @param cols2 Color for the highest Z values. Default is "navy".
#' @param alpha Transparency of the neuron visualization. Default is 0.5.
#' @param title.col Color of the plot title. Default is "darkgrey".
#' @param ... Additional arguments passed to geom_neuron().
#'
#' @return A ggplot2 object representing the neuron visualization.
#'
#' @details
#' This function wraps around geom_neuron() to create a complete plot with a
#' consistent, minimal theme. It removes axes, legends, and other extraneous
#' elements to focus on the neuron visualization itself.
#'
#' @examples
#' \dontrun{
#' # Visualize the banc volume
#' ggneuron(banc_neuropil.surf, banc.surf)
#'
#' # Visualize the banc brain neuropil
#' ggneuron(banc_neuropil.surf,
#' rotation_matrix = bancr:::banc_rotation_matrices[["front"]])
#'
#' # See constituents of a neuronlist
#' ggneuron(nat::Cell07PNs, volume = nat::MBL.surf,
#' root = 10, high = "red", low = "black")
#' }
#'
#' @seealso
#' \code{\link{geom_neuron}} for the underlying geom used in this function.
#'
#' @export
ggneuron <- function(x,
                     volume = NULL,
                     info = NULL,
                     rotation_matrix = NULL,
                     cols1 = c("turquoise","navy"),
                     cols2 =  c("grey75", "grey50"),
                     alpha = 0.5,
                     title.col = "darkgrey",
                     ...){
  ggplot2::ggplot() +
    {if(!is.null(volume)){
      geom_neuron(x = volume, rotation_matrix = rotation_matrix, alpha = max(alpha-0.25,0.01), cols = cols2)
    }} +
    geom_neuron(x = x, rotation_matrix = rotation_matrix, cols = cols1, alpha = alpha, ...) +
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
    ggplot2::labs(title = info)
}

prune_vertices.synapticneuron <- function (x, verticestoprune, invert = FALSE, ...){
  if(length(verticestoprune)==nrow(x$d)){
    warning('no points left after pruning')
    return(NULL)
  }
  soma <- catmaid::somaid(x)
  if(!is.null(soma)&&!is.na(soma)){
    x$d[catmaid::somaindex(x),"Label"] <- 1
  }
  pruned <- nat::prune_vertices(x, verticestoprune, invert = invert, ...)
  root <- nat::rootpoints(x)
  if(!is.null(root)){
    root <- nat::xyzmatrix(x)[root,]
    pruned <- nat::reroot(x = pruned, point = c(root))
  }
  pruned$connectors <- x$connectors[x$connectors$treenode_id %in%
                                     pruned$d$PointNo, ]
  relevant.points <- subset(x$d, x$d$PointNo %in% pruned$d$PointNo)
  y <- pruned
  y$d <- relevant.points[match(pruned$d$PointNo, relevant.points$PointNo),]
  y$d$Parent <-  pruned$d$Parent
  y$tags <- lapply(x$tags, function(t) t[t %in% pruned$d$PointNo])
  y$url <- x$url
  y$headers <- x$headers
  y$AD.segregation.index = x$AD.segregation.index
  smid <- which(y$d$Label==1)[1]
  if(length(smid)){
    y$tags$soma <- y$d$PointNo[smid]
  }
  class(y) <- class(x)
  y
}

