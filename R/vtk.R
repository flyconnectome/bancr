# And the same for the mesh
write_mesh3d_to_vtk <- function(mesh, filename) {
  if (!requireNamespace("rgl", quietly = TRUE)) {
    stop("Package 'rgl' is required but not installed.")
  }

  if (!inherits(mesh, "mesh3d")) {
    stop("Input must be a mesh3d object")
  }

  vertices <- t(mesh$vb[1:3,])
  faces <- t(mesh$it)

  cat("Vertex count:", nrow(vertices), "\n")
  cat("Face count:", nrow(faces), "\n")

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

  cat("Mesh successfully written to", filename, "\n")

  # Check file content
  cat("First few lines of the VTK file:\n")
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
  write(paste("v", vertices[,1], vertices[,2], vertices[,3]), con, sep = "\n")

  # Write texture coordinates if present
  if (!is.null(mesh$texcoords)) {
    texcoords <- t(mesh$texcoords)
    write(paste("vt", texcoords[,1], texcoords[,2]), con, sep = "\n")
  }

  # Write normals if present
  if (!is.null(mesh$normals)) {
    normals <- t(mesh$normals[1:3, ] / mesh$normals[4, ])
    write(paste("vn", normals[,1], normals[,2], normals[,3]), con, sep = "\n")
  }

  # Write faces
  faces <- t(mesh$it)
  if (!is.null(mesh$texcoords) && !is.null(mesh$normals)) {
    # Faces with vertex/texture/normal indices
    write(paste("f",
                paste(faces[,1], faces[,1], faces[,1], sep="/"),
                paste(faces[,2], faces[,2], faces[,2], sep="/"),
                paste(faces[,3], faces[,3], faces[,3], sep="/")), con, sep = "\n")
  } else if (!is.null(mesh$texcoords)) {
    # Faces with vertex/texture indices
    write(paste("f",
                paste(faces[,1], faces[,1], sep="/"),
                paste(faces[,2], faces[,2], sep="/"),
                paste(faces[,3], faces[,3], sep="/")), con, sep = "\n")
  } else if (!is.null(mesh$normals)) {
    # Faces with vertex//normal indices
    write(paste("f",
                paste(faces[,1], "", faces[,1], sep="/"),
                paste(faces[,2], "", faces[,2], sep="/"),
                paste(faces[,3], "", faces[,3], sep="/")), con, sep = "\n")
  } else {
    # Faces with only vertex indices
    write(paste("f", faces[,1], faces[,2], faces[,3]), con, sep = "\n")
  }

  cat("OBJ file written successfully:", filename, "\n")
}
