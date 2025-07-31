# And the same for the mesh
write_mesh3d_to_vtk <- function(mesh, filename, simplify = TRUE, percent = 0.1) {
  if (!requireNamespace("rgl", quietly = TRUE)) {
    stop("Package 'rgl' is required but not installed.")
  }

  if (!inherits(mesh, "mesh3d")) {
    stop("Input must be a mesh3d object")
  }

  # Clean and simplify meesh
  if(simplify){
    mesh <- Rvcg::vcgQEdecim(mesh, percent = percent)
    mesh <- Rvcg::vcgClean(mesh, sel=c(0,1,2,3,4,6,7))
  }

  vertices <- t(mesh$vb[1:3,])
  faces <- t(mesh$it)

  cat("Vertex count:", nrow(vertices), "
")
  cat("Face count:", nrow(faces), "
")

  con <- file(filename, "w")
  on.exit(close(con))

  writeLines("# vtk DataFile Version 2.0", con)
  writeLines("Mesh exported from R", con)
  writeLines("ASCII", con)
  writeLines("DATASET POLYDATA", con)

  writeLines(sprintf("POINTS %d float", nrow(vertices)), con)
  utils::write.table(format(vertices, scientific = FALSE), con, row.names = FALSE, col.names = FALSE, quote = FALSE)

  writeLines(sprintf("POLYGONS %d %d", nrow(faces), nrow(faces) * 4), con)
  face_data <- cbind(3, faces - 1)
  utils::write.table(face_data, con, row.names = FALSE, col.names = FALSE, quote = FALSE)

  cat("Mesh successfully written to", filename, "
")

  # Check file content
  cat("First few lines of the VTK file:
")
  system(sprintf("head -n 10 %s", filename))
}

# Write a OBJ file from mesh3d
write_mesh3d_to_obj <- function(mesh, filename) {
  if (!inherits(mesh, "mesh3d")) {
    stop("Input must be a mesh3d object")
  }

  # Open the file for writing
  con <- file(filename, "w")
  on.exit(close(con))  # Ensure the file is closed when the function exits

  # Write vertices
  vertices <- t(mesh$vb[1:3, ] / mesh$vb[4, ])
  write(paste("v", vertices[,1], vertices[,2], vertices[,3]), con, sep = "
")

  # Write texture coordinates if present
  if (!is.null(mesh$texcoords)) {
    texcoords <- t(mesh$texcoords)
    write(paste("vt", texcoords[,1], texcoords[,2]), con, sep = "
")
  }

  # Write normals if present
  if (!is.null(mesh$normals)) {
    normals <- t(mesh$normals[1:3, ] / mesh$normals[4, ])
    write(paste("vn", normals[,1], normals[,2], normals[,3]), con, sep = "
")
  }

  # Write faces
  faces <- t(mesh$it)
  if (!is.null(mesh$texcoords) && !is.null(mesh$normals)) {
    # Faces with vertex/texture/normal indices
    write(paste("f",
                paste(faces[,1], faces[,1], faces[,1], sep="/"),
                paste(faces[,2], faces[,2], faces[,2], sep="/"),
                paste(faces[,3], faces[,3], faces[,3], sep="/")), con, sep = "
")
  } else if (!is.null(mesh$texcoords)) {
    # Faces with vertex/texture indices
    write(paste("f",
                paste(faces[,1], faces[,1], sep="/"),
                paste(faces[,2], faces[,2], sep="/"),
                paste(faces[,3], faces[,3], sep="/")), con, sep = "
")
  } else if (!is.null(mesh$normals)) {
    # Faces with vertex//normal indices
    write(paste("f",
                paste(faces[,1], "", faces[,1], sep="/"),
                paste(faces[,2], "", faces[,2], sep="/"),
                paste(faces[,3], "", faces[,3], sep="/")), con, sep = "
")
  } else {
    # Faces with only vertex indices
    write(paste("f", faces[,1], faces[,2], faces[,3]), con, sep = "
")
  }

  cat("OBJ file written successfully:", filename, "
")
}

# hidden
write_neuron_to_vtk_paired <- function(neuron, file) {

  # Extract points from the neuron
  points <- nat::xyzmatrix(neuron)

  # Open the file for writing
  con <- file(file, "w")

  # Write VTK header
  writeLines("# vtk DataFile Version 3.0", con)
  writeLines("Neuron VTK file", con)
  writeLines("ASCII", con)
  writeLines("DATASET POLYDATA", con)

  # Write POINTS section
  writeLines(sprintf("POINTS %d float", nrow(points)), con)
  write.table(points, con, row.names = FALSE, col.names = FALSE, quote = FALSE)

  # Prepare LINES section
  line_data <- neuron$d %>%
    dplyr::mutate(Parent = Parent,
                  PointNo  = PointNo) %>%
    dplyr::mutate(from = (1:dplyr::n())-1,
                  to = from[match(Parent,PointNo)],
                  pair = 2) %>%
    dplyr::filter(!is.na(to)) %>%
    dplyr::select(pair, to, from)
  num_pairs <- nrow(line_data)

  # Write LINES section
  writeLines(sprintf("LINES %d %d", num_pairs, num_pairs * 3), con)
  write.table(line_data, con, row.names = FALSE, col.names = FALSE, quote = FALSE)

  # Close the file
  close(con)
  cat(sprintf("VTK file written to: %s
", file))
}
