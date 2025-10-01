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
#' dna01 <- banc_read_neuron_meshes(unique(banc.meta.dnao1$root_id))
#' banc.meta.dnao2 <- subset(banc.meta, cell_type=="DNa02")
#' dna02 <- banc_read_neuron_meshes(unique(banc.meta.dnao2$root_id))
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
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for this function. Please install it with: install.packages('ggplot2')")
  }
  if (!requireNamespace("ggpubr", quietly = TRUE)) {
    stop("Package 'ggpubr' is required for this function. Please install it with: install.packages('ggpubr')")
  }
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
      cols1 <- c("darkblue", "navy")
    }else{
      cols1 <- grDevices::colorRampPalette(c("#00008B","#0000CD","#4169E1","#1E90FF","#87CEEB", "#B0E0E6"))(length(neuron1))
    }
    if(length(neuron2)==1){
      cols2 <- c("darkred", "#F88379")
    }else{
      cols2 <- grDevices::colorRampPalette(c("#8f0723","#DC143C","#FF4500","#FF7F50","#F88379","#FFB6C1","#FF69B4"))(length(neuron2))
    }
    if(length(neuron3)==1){
      cols3 <- c("darkgreen", "green")
    }else{
      cols3 <- grDevices::colorRampPalette(c("#006400","#076b3e","#228B22", "#32c080","#90EE90"))(length(neuron3))
    }
    if(length(alpha)<3){
      alpha<-rep(alpha[1],3)
    }

    # Create the plot
    p <- ggplot2::ggplot() +
      nat.ggplot::geom_neuron(x = mesh, rotation_matrix = rotation_matrix, alpha = 0.05, cols = c("grey90", "grey60")) +
      nat.ggplot::geom_neuron(x = volume, rotation_matrix = rotation_matrix, alpha = 0.05, cols = c("grey50")) +
      nat.ggplot::geom_neuron(x=neuron_pruned1, rotation_matrix = rotation_matrix, cols = cols1, alpha = alpha[1], linewidth = 0.3) +
      nat.ggplot::geom_neuron(x=neuron_pruned2, rotation_matrix = rotation_matrix, cols = cols2, alpha = alpha[2], linewidth = 0.3) +
      nat.ggplot::geom_neuron(x=neuron_pruned3, rotation_matrix = rotation_matrix, cols = cols3, alpha = alpha[3], linewidth = 0.3) +
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
    ggplot2::ggsave(filename,
                    plot = ga,
                    width = width,
                    height = height,
                    bg = "#fcfcfa",
                    dpi = 300,
                    limitsize = FALSE,
                    create.dir = TRUE)
    message("saved as:  ", filename)
    invisible()
  }else{
    plot(ga)
  }
}


