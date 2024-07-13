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

#' Perform Elastix Transform on 3D Points
#'
#' This function applies an Elastix spatial transform to a set of 3D points.
#'
#' @param points A matrix with 3 columns or a data frame with x, y, z columns representing 3D points.
#' @param transform_file Path to the Elastix transform file, usually a `.txt` file, usually a `.txt` file.
#' @param copy_files Vector of additional file paths to copy to the temporary directory.
#' @param return_logs Logical, if TRUE, returns the Elastix log instead of transformed points.
#'
#' @return A matrix of transformed 3D points, or Elastix logs if return_logs is TRUE.
#'
#' @details
#' This function requires Elastix to be installed and added to the system PATH.
#' It creates a temporary directory for processing, applies the Elastix transform,
#' and cleans up afterwards.
#'
#' @examples
#' \dontrun{
#' points <- matrix(rnorm(30), ncol = 3)
#' transformed_points <- elastix_xform(points, "path/to/transform.txt")
#' }
#'
#' @export
elastix_xform <- function(points,
                          transform_file,
                          copy_files = c(),
                          return_logs = FALSE) {

  # Does transform file exist
  check_if_possible(transform_file)

  # Do we have valid points
  points <- nat::xyzmatrix(points)
  if (ncol(points) != 3) {
    stop("points must be a matrix with 3 columns or a data frame with x/y/z columns")
  }

  # Create a temporary directory
  temp_dir <- tempdir(check = TRUE)

  # Copy additional files if required
  if (length(copy_files) > 0) {
    file.copy(copy_files, temp_dir)
  }

  # Write points to file
  in_file <- file.path(temp_dir, "inputpoints.txt")
  write_elastix_input_file(points, in_file)
  out_file <- file.path(temp_dir, "outputpoints.txt")

  # Prepare the command
  elastix_path <- Sys.which("transformix")
  if (elastix_path == "") {
    elastix_path <- "/opt/elastix-5.1.0-mac/bin/transformix"
    #stop("Could not find elastix binary. Make sure it's in your PATH.")
  }

  # Construct system command
  command <- paste(
    elastix_path,
    "-out", temp_dir,
    "-tp", transform_file,
    "-def", in_file
  )

  # Run the transform
  sys.run <- system(command, intern = TRUE)
  if (return_logs) {
    log_file <- file.path(temp_dir, "transformix.log")
    if (!file.exists(log_file)) {
      warning("No log file found.")
      stop(sys.run)
    }
    return(readLines(log_file))
  }
  if (!file.exists(out_file)) {
    warning("Elastix transform did not produce any output.")
    stop(sys.run)
  }

  # Parse points
  points_xf <- read_elastix_output_file(out_file)

  # Clean up
  unlink(temp_dir, recursive = TRUE)

  # Return
  return(points_xf)
}

#' Check if Elastix Transform is possible
#'
#' Verifies if the specified transform file exists.
#'
#' @param file Path to the Elastix transform file.
#' @param on_error Action to take on error: "raise" to stop execution, or any other value to return an error message.
#'
#' @return NULL if file exists, or an error message if the file doesn't exist and on_error is not "raise".
#'
#' @keywords internal
check_if_possible <- function(file, on_error = "raise") {
  if (!file.exists(file)) {
    msg <- paste("Transformation file", file, "not found.")
    if (on_error == "raise") {
      stop(msg)
    }
    return(msg)
  }
}

#' Write Elastix Input File
#'
#' Writes 3D points to a file in the format required by Elastix.
#'
#' @param points Matrix of 3D points.
#' @param filepath Path where the input file should be written.
#'
#' @keywords internal
#' @importFrom utils write.table
write_elastix_input_file <- function(points, filepath) {
  cat("point\n", nrow(points), "\n", file = filepath)
  utils::write.table(points, filepath, append = TRUE, col.names = FALSE,
              row.names = FALSE, sep = " ")
}

#' Read Elastix Output File
#'
#' Reads and parses the output file produced by Elastix transform.
#'
#' @param filepath Path to the Elastix output file.
#'
#' @return A matrix of transformed 3D points.
#'
#' @keywords internal
read_elastix_output_file <- function(filepath) {

  # Read all lines from the file
  lines <- readLines(filepath)

  # Process each line
  points <- lapply(lines, function(line) {
    # Extract the part between 'OutputPoint = [' and ']'
    output <- strsplit(strsplit(line, "OutputPoint = \\[ ")[[1]][2], " \\]")[[1]][1]

    # Split the string into numeric values
    as.numeric(strsplit(output, " ")[[1]])
  })

  # Convert the list of points to a matrix
  points_matrix <- do.call(rbind, points)

  # Return
  return(points_matrix)
}

# hidden
update_elastix_transforms_locations <- function(transform_file,
                                                search = "1_elastix_affine",
                                                file_path = NULL){
   # Read the file content
  lines <- readLines(transform_file)

  # Define the target pattern and replacement
  target_pattern <- sprintf('(InitialTransformParametersFileName\\s*")(.+%s\\.txt)(")', search)

  # Function to replace the target pattern
  replace_path <- function(line) {
    gsub(target_pattern, paste0("\\1", file_path, "\\3"), line)
  }

  # Apply the replacement to each line
  updated_lines <- sapply(lines, replace_path)

  # Write the updated content back to the file
  writeLines(updated_lines, transform_file)
}

#' Apply Elastix Transform using Navis
#'
#' Applies an Elastix transform to 3D points using the Navis Python library.
#'
#' @param x Matrix or data frame of 3D points.
#' @param transform_file Path to the Elastix transform file, usually a `.txt` file.
#'
#' @return A matrix of transformed 3D points.
#'
#' @details
#' This function requires the reticulate R package and the Navis Python library.
#'
#' @examples
#' \dontrun{
#' neuron.mesh <- banc_read_neuron_meshes("720575941478275714")
#' points <- nat::xyzmatrix(neuron.mesh)
#' transformed_points <- navis_elastix_xform(points,
#' transform_file = "brain_240707/BANC_to_template.txt")
#' points3d(points)
#' plot3d(nat.flybrains::JRC2018F)
#' }
#'
#' @export
navis_elastix_xform <- function(x, transform_file){
  xyz <- nat::xyzmatrix(x)
  if (ncol(xyz) != 3) {
    stop("Input 'x' must have exactly 3 columns representing x, y and z coordinates.")
  }
  reticulate::py_run_string("from navis import transforms")
  reticulate::py_run_string(sprintf("tr = transforms.ElastixTransform('%s')",
                                    transform_file))
  reticulate::py_run_string("xform = tr.xform")
  result <- reticulate::py$xform(x)
  colnames(result) <- colnames(xyz)
  nat::xyzmatrix(x) <- result
  x
}

#' Transform Points between BANC Connectome and JRC2018F Template Brain
#'
#' @description
#' This function transforms 3D points between the BANC (Buhmann et al. Adult Neural Connectome)
#' coordinate system and the D. melanogaster template brain JRC2018F coordinate system.
#'
#' @param x An object containing 3D points (must be compatible with nat::xyzmatrix).
#' @param banc.units Character string specifying the units of the BANC space data (input or output, depending on the inverse argument).
#'   Must be one of "nm" (nanometers), "um", or "raw" (BANC raw banc.units). Default is "nm".
#' @param subset Optional. A logical vector or expression to subset the input object.
#' @param inverse Logical. If TRUE, performs the inverse transformation (JRC2018F to BANC).
#'   Default is FALSE.
#' @param transform_file Optional. Path to a custom transform file. If NULL, uses default files.
#' @param method Character string specifying the transformation method.
#'   Must be either "elastix" or "tpsreg". Default is "elastix".
#'
#' @return The input object with transformed 3D points.
#'
#' @details
#' This function applies either an Elastix transform or a thin-plate spline registration to convert
#' points between the BANC and JRC2018F coordinate systems. It handles unit conversions as necessary.
#'
#' The default transformation files are included with the package and are located in the
#' 'inst/extdata/brain_240707' directory.
#'
#' @examples
#' \dontrun{
#' # Transform points from BANC to JRC2018F
#' transformed_points <- banc_to_JRC2018F(points, banc.units = "nm")
#'
#' # Use a custom transform file
#' custom_transformed <- banc_to_JRC2018F(points, transform_file = "path/to/custom/transform.txt")
#'
#' # Where the default transform files are located:
#' banc_to_JRC2018F_file <- system.file(file.path("extdata","brain_240707"),
#' "BANC_to_template.txt", package="bancr")
#' JRC2018F_to_banc_file <- system.file(file.path("extdata","brain_240707"),
#' "template_to_BANC.txt", package="bancr")
#' }
#'
#' @seealso
#' \code{\link{elastix_xform}} for the underlying Elastix transformation function.
#' \code{\link{banc_raw2nm}} and \code{\link{banc_nm2raw}} for unit conversion functions.
#'
#' @export
banc_to_JRC2018F <- function(x,
                             banc.units = c("nm", "um", "raw"),
                             subset = NULL,
                             inverse = FALSE,
                             transform_file = NULL,
                             method = c("tpsreg","elastix","navis_elastix_xform")){

  # manage arguments
  banc.units <- match.arg(banc.units)
  method <- match.arg(method)

  # find transform
  if(is.null(transform_file)){
    if(inverse){
      transform_file <- system.file(file.path("extdata","brain_240707"), "3_elastix_Bspline_fine.txt", package="bancr")
      transform_file2 <- system.file(file.path("extdata","brain_240707"), "2_elastix_Bspline_coarse.txt", package="bancr")
      transform_file1 <- system.file(file.path("extdata","brain_240707"), "1_elastix_affine.txt", package="bancr")
      update_elastix_transforms_locations(transform_file, search = "2_elastix_Bspline_coarse", file_path = system.file(file.path("extdata","brain_240707"), package="bancr"))
      update_elastix_transforms_locations(transform_file2, search = "1_elastix_affine", file_path = system.file(file.path("extdata","brain_240707"), package="bancr"))
      update_elastix_transforms_locations(transform_file1, search = "0_manual_affine", file_path = system.file(file.path("extdata","brain_240707"), package="bancr"))
    }else{
      transform_file <- system.file(file.path("extdata","brain_240707"), "BANC_to_template.txt", package="bancr")
    }
  }

  # apply subset
  if(!is.null(subset)) {
    xs <- x[subset]
    xst <- banc_to_JRC2018F(xs, banc.units = banc.units, inverse=inverse, transform_file=transform_file, method=method)
    x[subset] <- xst
    return(x)
  }

  # convert to um if necessary
  xyz <- nat::xyzmatrix(x)
  if(isFALSE(inverse) && grepl("elastix",method)){
    if(banc.units=='nm'){
      xyz <- xyz/1e3
    }else if(banc.units=='raw'){
      xyz <- banc_raw2nm(xyz)/1e3
    }
  }

  ## Decapitate
  # if(isFALSE(inverse)){
  #   xyz <- banc_decapitate(xyz*100)/1000
  # }

  # do transformation
  if(method=="elastix"){
    # Result is in um
    xyz2 <- elastix_xform(xyz, transform_file = transform_file)
  }else if(method=="navis_elastix_xform"){
    # Result is in um
    xyz2 <- navis_elastix_xform(xyz, transform_file = transform_file)
  }else{
    if(inverse){
      # Result is in nm
      # utils::data("jrc2018f_to_banc_tpsreg", envir = environment())
      xyz2 <- Morpho::applyTransform(xyz,
                                     trafo = bancr::jrc2018f_to_banc_tpsreg,
                                     inverse = FALSE)
    }else{
      # Result is in um
      # utils::data("banc_to_jrc2018f_tpsreg", envir = environment())
      xyz2 <- Morpho::applyTransform(xyz,
                                     trafo = bancr::banc_to_jrc2018f_tpsreg,
                                     inverse = FALSE)
    }
  }

  ## Decapitate
  # if(isTRUE(inverse)){
  #   xyz2 <- banc_decapitate(xyz2*100)/1000
  # }

  # convert from um to original banc.units if necessary
  if(isTRUE(inverse) && grepl("elastix",method)){
    if(banc.units=='nm'){
      xyz2 <- xyz2*1e3
    }else if(banc.units=='raw'){
      xyz2 <- banc_nm2raw(xyz2*1e3)
    }
  }else if(isTRUE(inverse)){
    if(banc.units=='um'){
      xyz2 <- xyz2/1e3
    }else if(banc.units=='raw'){
      xyz2 <- banc_nm2raw(xyz2)
    }
  }

  # put points back
  nat::xyzmatrix(x) <- xyz2

  # return object
  return(x)
}
# Jasper's Elastix transform
# transform_file <- "/Users/GD/LMBD/Papers/banc/the-BANC-fly-connectome/fanc/transforms/transform_parameters/brain_240707/BANC_to_template.txt"

#' Mirror BANC Connectome Points
#'
#' @description
#' This function mirrors 3D points in the BANC (Buhmann et al. Adult Neural Connectome)
#' coordinate system by transforming to JRC2018F, mirroring, and transforming back.
#'
#' @param x An object containing 3D points (must be compatible with nat::xyzmatrix).
#' @param banc.units Character string specifying the banc.units of the input points.
#'   Must be one of "nm" (nanometers), "um", or "raw" (BANC raw banc.units). Default is "nm".
#' @param subset Optional. A logical vector or expression to subset the input object.
#' @param inverse Logical. Not used in this function, kept for compatibility with banc_to_JRC2018F.
#' @param transform_files Optional. A vector of two file paths for custom transform files.
#'   If NULL, uses default files.
#' @param method Character string specifying the transformation method.
#'   Must be either "elastix" or "tpsreg". Default is "elastix".
#' @param ... Additional arguments passed to \code{nat.templatebrains::\link{mirror_brain}}.
#'
#' @return The input object with mirrored 3D points.
#'
#' @details
#' This function performs mirroring of BANC points by first transforming them to the JRC2018F
#' coordinate system, applying the mirroring operation, and then transforming them back to BANC.
#' It can use either Elastix transforms or thin-plate spline registration for the coordinate
#' system transformations.
#'
#' @examples
#' \dontrun{
#' # Example using saved tpsreg
#' banc_neuropil.surf.m <- banc_mirror(banc_neuropil.surf, method = "tpsreg")
#' clear3d()
#' banc_view()
#' plot3d(banc_neuropil.surf, alpha = 0.5, col = "lightgrey")
#' plot3d(banc_neuropil.surf.m, alpha = 0.5, col = "green")
#'
#' # Example using custom Elastix transforms
#' choose_banc()
#' rootid <- "720575941626035769"
#' neuron.mesh <- banc_read_neuron_meshes(rootid)
#'
#' # Show neuron in BANC neuropil
#' banc_view()
#' plot3d(neuron.mesh, col = "red")
#' plot3d(banc_neuropil.surf, alpha = 0.1, col = "lightgrey")
#'
#' # Show only the portion in the brain
#' neuron.mesh.brain <- banc_decapitate(neuron.mesh, invert = TRUE)
#'
#' # Mirror in BANC space
#' neuron.mesh.mirror <- banc_mirror(neuron.mesh.brain,
#' transform_files = c("brain_240707/BANC_to_template.txt",
#'  "brain_240707/template_to_BANC.txt"))
#' plot3d(neuron.mesh.mirror, col = "cyan")
#' }
#'
#' @seealso
#' \code{\link{banc_to_JRC2018F}} for the underlying transformation function.
#' \code{nat.templatebrains::\link{mirror_brain}} for the mirroring operation in JRC2018F space.
#'
#' @export
banc_mirror <- function(x,
                        banc.units = c("nm", "um", "raw"),
                        subset = NULL,
                        inverse = FALSE,
                        transform_files = NULL,
                        method = c("tpsreg","elastix","navis_elastix_xform"),
                        ...){

  # Manage arguments
  banc.units <- match.arg(banc.units)
  method <- match.arg(method)

  #Get 3D points
  xyz <- nat::xyzmatrix(x)

  # Use elastix transform
  if(method=="elastix"){

    # Convert to JRC2018F
    x.jrc2018f <- banc_to_JRC2018F(x=xyz, banc.units=banc.units, subset=NULL, inverse=FALSE, transform_file = transform_files[1], method = method)

    # Mirror
    x.jrc2018f.m <-  nat.templatebrains::mirror_brain(x.jrc2018f, brain = nat.flybrains::JRC2018F, transform = "flip")

    # Back to BANC
    x.banc.m <- banc_to_JRC2018F(x=x.jrc2018f.m, banc.units=banc.units, subset=NULL, inverse=TRUE , transform_file = transform_files[2], method = method)

  }else{

    # convert to um if necessary
    if(banc.units=='um'){
      xyz <- xyz*1e3
    }else if(banc.units=='raw'){
      xyz <- banc_raw2nm(xyz)
    }

    # use pre-calculated tps reg
    # utils::data("banc_mirror_tpsreg", envir = environment())
    x.banc.m <- Morpho::applyTransform(xyz,
                                       trafo = bancr::banc_mirror_tpsreg)

    # convert from um to original banc.units if necessary
    if(banc.units=='um'){
      x.banc.m <- x.banc.m/1e3
    }else if(banc.units=='raw'){
      x.banc.m <- banc_nm2raw(x.banc.m)
    }

  }

  # return
  nat::xyzmatrix(x) <- x.banc.m
  return(x)

}
