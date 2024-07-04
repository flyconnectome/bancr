banc_neuron_comparison_plot <- function(neuron1, neuron2,
                                        filename =  NULL,
                                        width = 8, height = 6) {

  # Get 3D spatial points
  glist <- list()
  for(view in names(banc_rotation_matrices)){

    # Choose mesh
    if(view=="front"){
      mesh <- banc_brain_neuropil.surf
      decaptitate <- "brain"
    }else if(view=="vnc"){
      mesh <- banc_vnc_neuropil.surf
      decaptitate <- "vnc"
    }else{
      mesh <- banc_neuropil.surf
      decaptitate <- "none"
    }

    # Get 3D points
    vertices <- nat::xyzmatrix(mesh)
    if(decaptitate=="brain"){
      neuron_pruned1 <- banc_decapitate(neuron1, invert = TRUE)
      neuron_pruned2 <- banc_decapitate(neuron2, invert = TRUE)
    }else if(decaptitate=="vnc"){
      neuron_pruned1 <- banc_decapitate(neuron1, invert = FALSE)
      neuron_pruned2 <- banc_decapitate(neuron2, invert = FALSE)
    }else{
      neuron_pruned1=neuron1
      neuron_pruned2=neuron2
    }

    # Rotate mesh
    rotation_matrix <- banc_rotation_matrices[[view]]
    rotated_vertices <- as.data.frame(t(rotation_matrix[,1:3] %*% t(as.matrix(vertices))))
    rotated_vertices <- rotated_vertices[,-4]
    colnames(rotated_vertices) <- c("X","Y","Z")

    # Process neuron 1
    edges_df1 <- ggplot2_neuron_path(neuron_pruned1, rotation_matrix=rotation_matrix)

    # Process neuron 2
    edges_df2 <- ggplot2_neuron_path(neuron_pruned2, rotation_matrix=rotation_matrix)

    # Create the plot
    p <- ggplot() +
      geom_point(data = rotated_vertices, aes(x = X, y = Y), color = "grey", size = 1, alpha = 0.05) +
      #scale_color_gradient(low = "grey100", high = "grey50") +
      ggnewscale::new_scale_colour() +
      geom_path(data = edges_df1, aes(x = X, y = Y, color = Z, group = group), alpha = 0.5, linewidth = 0.5) +
      scale_color_gradient(low = "navy", high = "cyan") +
      ggnewscale::new_scale_colour() +
      geom_path(data = edges_df2, aes(x = X, y = Y, color = Z, group = group), alpha = 0.5, linewidth = 0.5) +
      scale_color_gradient(low = "coral", high = "darkred") +
      theme_minimal() +
      coord_equal() +
      labs(title = view) +
      theme(legend.position = "none") +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank())

    # Add to list
    glist[[view]] <- p
  }

  # Arrange plot
  g <- gridExtra::grid.arrange(glist[["main"]], glist[["side"]], glist[["front"]], glist[["vnc"]], ncol = 2,
                               layout_matrix = rbind(c(1,1,1,2,2,3,3,3,3),
                                                     c(1,1,1,2,2,3,3,3,3),
                                                     c(1,1,1,2,2,4,4,4,4),
                                                     c(1,1,1,2,2,4,4,4,4)
                               ))
  if(filename){
    ggsave(filename, plot = p, width = width, height = height)
    invisible()
  }else{
    plot(g)
  }
}

#' Convert Neuron Objects to ggplot2-Compatible Path Data
#'
#' @description
#' This function converts 'neuron' or 'neuronlist' objects, which represent 3D
#' points linked by lines in space, into data frames that describe paths
#' compatible with ggplot2's geom_path.
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
ggplot2_neuron_path.neuron <- function(x, rotation_matrix = NULL){
  x$d <- x$d[order(x$d$Parent),]
  x$d <- x$d[order(x$d$PointNo),]
  npoints <- nat::xyzmatrix(x)
  edges_df <- data.frame()
  if(!is.null(rotation_matrix)){
    npoints <- as.data.frame(t(rotation_matrix[,1:3] %*% t(xyzmatrix(npoints))))
    npoints <- npoints[,-4]
    colnames(npoints) <- c("X","Y","Z")
  }
  ss <-nat::seglist(x)
  seglist <- ss[[1]]$SegList
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
ggplot2_neuron_path.neuronlist <- function(x, rotation_matrix = NULL){
  ll <- lapply(x, ggplot2_neuron_path, rotation_matrix = rotation_matrix, ...)
  max.group <- 0
  for(i in 1:length(ll)){
    ll[[i]]$group <- ll[[i]]$group + max.group
    max.group <- max(ll[[i]]$group, na.rm = TRUE)
  }
  do.call(rbind, ll)
}



