#' Upload Mesh to Google Cloud Storage for BANC Neuroglancer
#'
#' @description This function uploads a mesh.obj file to a Google Cloud Storage bucket
#' that can be accessed by BANC neuroglancer. It uses reticulate to call
#' Python code and manages a conda environment.
#'
#' The underlying Python code is based on the bikinibottom package:
#' \url{https://github.com/jasper-tms/bikini-bottom/tree/main}
#'
#' @param mesh Character or mesh3d object. Either a path to the .obj file or an R mesh object.
#' @param mesh_id Integer. ID for the mesh.
#' @param vol Character. URL for the precomputed volume.
#' @param compress Logical. Whether to compress the mesh. Default is `TRUE`.
#' @param overwrite Logical. Whether or not to overwrite an extant file of the same name and path at
#' the location specified by `vol`+`mesh_id`. Default is `FALSE`.
#  @param conda_env Character. Name of the conda environment. Default is "banc".
#  @param python_version Character. Python version to use. Default is "3.8".
#  @param force_env_creation Logical. If TRUE, forces creation of a new environment even if it exists. Default is FALSE.
#'
#' @return Invisible `NULL`. The function is called for it upload effect.
#'
#' @details For this function to work the global environment variable `GOOGLE_APPLICATION_CREDENTIALS` to
#' something like `"$HOME/.config/gcloud/legacy_credentials/firstname.secondname@gmail.com/adc.json"`.
#' This file can be obtains by installing `gsutil` and running `gsutil init` in the terminal.
#'
#' @importFrom reticulate use_condaenv conda_create conda_install py_install import conda_list
#' @importFrom utils packageVersion
#' @export
#'
#' @examples
#' \dontrun{
#' banc_upload_mesh(
#' mesh = "/Users/abates/projects/flyconnectome/bancpipeline/deformetrica/obj/banc_brain_neuropil.obj",
#' mesh_id = 3,
#' vol = "precomputed://gs://lee-lab_brain-and-nerve-cord-fly-connectome/volume_meshes/meshes"
#' )
#' }
#' @references
#' bikinibottom package: \url{https://github.com/jasper-tms/bikini-bottom/tree/main}
banc_upload_mesh <- function(mesh,
                             mesh_id,
                             vol = 'precomputed://gs://lee-lab_brain-and-nerve-cord-fly-connectome/volume_meshes/',
                             compress = TRUE,
                             overwrite = FALSE) {
    # Had tried to use: conda_env = "banc", python_version = "3.9", force_env_creation = FALSE

    # Input validation
    if (!is.character(mesh) && !inherits(mesh, "mesh3d")) {
      stop("mesh must be either a character string (file path) or a mesh3d object.")
    }
    if (is.character(mesh) && !file.exists(mesh)) {
      stop("Invalid mesh path. Please provide a valid path to an existing .obj file.")
    }
    if (!is.character(vol) || !grepl("^precomputed://", vol)) {
      stop("Invalid vol. It should be a character string starting with 'precomputed://'.")
    }
    if (!is.logical(compress)) {
      stop("compress must be a logical value (TRUE or FALSE).")
    }
    # if (!is.character(conda_env)) {
    #   stop("conda_env must be a character string.")
    # }
    # if (!is.character(python_version)) {
    #   stop("python_version must be a character string.")
    # }
    # if (!is.logical(force_env_creation)) {
    #   stop("force_env_creation must be a logical value (TRUE or FALSE).")
    # }

    # Check if reticulate is installed
    if (!requireNamespace("reticulate", quietly = TRUE)) {
      stop("Please install the 'reticulate' package to use this function.")
    }

    # Check reticulate version
    if (packageVersion("reticulate") < "1.16") {
      warning("Your reticulate version might be outdated. Consider updating for better compatibility.")
    }

    # # Function to create and setup the conda environment
    # setup_conda_env <- function() {
    #   tryCatch({
    #     reticulate::conda_create(conda_env, python_version = python_version)
    #     message("Created conda environment: ", conda_env)
    #
    #     reticulate::conda_install(conda_env, packages = c("numpy", "trimesh"))
    #     reticulate::conda_install(conda_env, packages = c("pip"))
    #     reticulate::conda_install(conda_env, pip = TRUE, packages = c(
    #       "cloud-volume==8.32.1",
    #       "git+https://github.com/jasper-tms/npimage.git",
    #       "bikinibottom"
    #     ))
    #
    #     message("Installed required packages in the conda environment.")
    #   }, error = function(e) {
    #     stop("Failed to create or setup conda environment: ", e$message)
    #   })
    # }

    # # Check if conda environment exists and create if necessary
    # if (!conda_env %in% reticulate::conda_list()$name || force_env_creation) {
    #   setup_conda_env()
    # } else {
    #   message("Using existing conda environment: ", conda_env)
    # }

    # # Try to use the conda environment, if it fails, attempt to reset Python configuration
    # tryCatch({
    #   reticulate::use_condaenv(conda_env, required = TRUE)
    # }, error = function(e) {
    #   message("Attempting to reset Python configuration...")
    #   current_config <- reticulate::py_config()
    #   reticulate::py_discover_config()
    #   new_config <- reticulate::py_config()
    #   if (identical(current_config, new_config)) {
    #     stop("Unable to reset Python configuration. Please restart your R session.")
    #   } else {
    #     reticulate::use_condaenv(conda_env, required = TRUE)
    #   }
    # })

    # Check if npimage and bikinibottom are installed, install if not
    check_and_install_packages <- function() {
      installed_packages <- reticulate::py_list_packages()

      if (!"numpyimage" %in% installed_packages$package) {
        message("Installing npimage...")
        reticulate::py_install("git+https://github.com/jasper-tms/npimage.git", pip = TRUE)
      }

      if (!"bikinibottom" %in% installed_packages$package) {
        message("Installing bikinibottom...")
        reticulate::py_install("git+https://github.com/jasper-tms/bikini-bottom.git", pip = TRUE)
      }
    }

    # Actually check and install
    check_and_install_packages()

    # Import bikinibottom
    tryCatch({
      bb <- reticulate::import("bikinibottom")
    }, error = function(e) {
      stop("Failed to import bikinibottom. Error: ", e$message)
    })

    # Convert R mesh to Python trimesh if necessary
    if (inherits(mesh, "mesh3d")) {
      message("Converting R mesh to Python trimesh...")
      clean_mesh <- Rvcg::vcgClean(mesh, sel = 1)  # This removes unreferenced vertices
      vertices <- t(clean_mesh$vb[-4,])  # Remove homogeneous coordinate and transpose
      faces <- t(clean_mesh$it)  # Transpose to get correct orientation
      faces <- t(clean_mesh$it - 1)  # Subtract 1 to convert to 0-based indexing

      # check
      message("Vertices shape: ", paste(dim(vertices), collapse = " x "))
      message("Faces shape: ", paste(dim(faces), collapse = " x "))
      message("Max face index: ", max(faces))
      message("Number of vertices: ", nrow(vertices))

      # Create Python trimesh object
      trimesh <- reticulate::import("trimesh")
      py_mesh <- trimesh$Trimesh(vertices = vertices, faces = faces)
    } else {
      py_mesh <- mesh  # If it's a file path, pass it directly
    }

    # Define a Python function that will receive and process the large integer
    if(bit64::is.integer64(mesh_id)){
      reticulate::py_run_string("
import numpy as np
def process_large_integer(int64_str):
    # Convert the string to a 64-bit integer using numpy
    int64 = np.int64(int64_str)
    print(f'Received 64-bit integer: {int64}')
    print(f'Type in Python: {type(int64)}')
    return int64
")
    }else if(is.character(mesh_id)){
      reticulate::py_run_string("
def process_large_integer(large_int_str):
    # Convert the string back to an integer in Python
    large_int = int(large_int_str)
    print(f'Received integer: {large_int}')
    print(f'Type in Python: {type(large_int)}')
    # Your processing logic here
    return large_int  # or whatever you want to return
")
    }

    # Pass integer to bb
    reticulate::py_run_string("
import bikinibottom as bb
def send_mesh(mesh, mesh_id, vol, compress, overwrite):
    mesh_id = process_large_integer(mesh_id)
    bb.push_mesh(mesh=mesh, mesh_id=mesh_id, vol=vol, compress=compress, overwrite=overwrite)
")

    # Upload mesh
    tryCatch({
      message("Uploading mesh_id: ", mesh_id)
      reticulate::py$send_mesh(mesh = py_mesh, mesh_id = mesh_id,
                   vol = vol, compress = compress, overwrite = overwrite)
      message("Mesh uploaded successfully.")
    }, error = function(e) {
      if(overwrite){
        message("Error uploading mesh: ", e$message)
      }else{
        warning("Error uploading mesh: ", e$message)
      }
    })

    # Return
    invisible(NULL)
  }


# import bikinibottom
# bikinibottom.push_mesh(‘banc_neuropil.obj’,
#                         mesh_id=123,
#                         vol=‘gs://lee-lab_brain-and-nerve-cord-fly-connectome/volume_meshes’)

# Write BANC neuropil as neuroglancer mesh
# ml <- neuronlist(as.mesh3d(banc_neuropil.surf))
# names(ml) <- "banc_neuropil"
# write_neuroglancer_mesh_dir(ml,"data-raw/ngmesh")


# # # googleCloudStorageR
# banc_gcs_read <- function(path = "gs://zetta_lee_fly_cns_001_synapse/240623_run/assignment/final_edgelist.df"){
#
# # Not sure hot to get this working in R, but looks useful
# googleCloudStorageR::gcs_setup(token="none")
# googleCloudStorageR::gcs_auth(token=NULL)
# #
# # # OR?
# scope <-c("https://www.googleapis.com/auth/cloud-platform")
# token <- gargle::token_fetch(scopes = scope)
# googleCloudStorageR::gcs_auth(token = token)
#
#
# # Then get data?
# gcs_auth(json_file="/Users/abates/.config/gcloud/legacy_credentials/alexander.shakeel.bates@gmail.com/adc.json")
# path <- 'gs://zetta_lee_fly_cns_001_synapse/240529_run/240604_assignment/final_edgelist.df'
# path <- 'gs://lee-lab_brain-and-nerve-cord-fly-connectome/volume_meshes/meshes/1:0'
# df <- googleCloudStorageR::gcs_get_object(path,
#                                           bucket = "lee-lab_brain-and-nerve-cord-fly-connectome",
#                                           )



# import bikinibottom as bb
# bb.push_mesh(mesh='/Users/abates/projects/flyconnectome/bancpipeline/deformetrica/obj/banc_neuropil.obj', mesh_id = 2,vol="precomputed://gs://lee-lab_brain-and-nerve-cord-fly-connectome/volume_meshes/meshes")























