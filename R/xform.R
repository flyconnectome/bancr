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
#' @param transform_file Path to the Elastix transform file.
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
elastix_xform <- function(points, transform_file, copy_files = c(), return_logs = FALSE) {
  check_if_possible(transform_file)

  # Do we have valid points
  points <- nat::xyzmatrix(points)
  if (ncol(points) != 3) {
    stop("points must be a matrix with 3 columns or a data frame with x/y/z columns")
  }

  # Create a temporary directory
  temp_dir <- tempdir()

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

#' Check if Elastix Transform is Possible
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
write_elastix_input_file <- function(points, filepath) {
  cat("point\n", nrow(points), "\n", file = filepath)
  write.table(points, filepath, append = TRUE, col.names = FALSE,
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

  return(points_matrix)
}

#' Apply Elastix Transform using Navis
#'
#' Applies an Elastix transform to 3D points using the Navis Python library.
#'
#' @param x Matrix or data frame of 3D points.
#' @param transform_path Path to the Elastix transform file.
#'
#' @return A matrix of transformed 3D points.
#'
#' @details
#' This function requires the reticulate R package and the Navis Python library.
#'
#' @examples
#' \dontrun{
#' points <- matrix(rnorm(30), ncol = 3)
#' transformed_points <- navis_elastix_xform(points, "path/to/transform.txt")
#' }
#'
#' @export
navis_elastix_xform <- function(x, transform_path){
  x <- nat::xyzmatrix(x)
  if (ncol(x) != 3) {
    stop("Input 'x' must have exactly 3 columns representing x, y, and z coordinates.")
  }
  reticulate::py_run_string("from navis import transforms")
  reticulate::py_run_string(sprintf("tr = transforms.ElastixTransform('%s')",
                                    transform_path))
  reticulate::py_run_string("xform = tr.xform")
  result <- reticulate::py$xform(r_matrix)
  colnames(result) <- colnmes(x)
  result
}

# Jasper's elastix transform
# transform_file <- "/Users/GD/LMBD/Papers/banc/the-BANC-fly-connectome/fanc/transforms/transform_parameters/brain_240707/BANC_to_template.txt"
# Transforms, note: elastix transforms cannot be inverted
banc_to_JRC2018F <- system.file(file.path("extdata","brain_240707"), "BANC_to_template.txt", package="bancr")
JRC2018F_to_banc <- system.file(file.path("extdata","brain_240707"), "template_to_BANC.txt", package="bancr")
