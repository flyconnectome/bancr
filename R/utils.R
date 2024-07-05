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
  write.table(format(vertices, scientific = FALSE), con, row.names = FALSE, col.names = FALSE, quote = FALSE)

  writeLines(sprintf("POLYGONS %d %d", nrow(faces), nrow(faces) * 4), con)
  face_data <- cbind(3, faces - 1)
  write.table(face_data, con, row.names = FALSE, col.names = FALSE, quote = FALSE)

  cat("Mesh successfully written to", filename, "\n")

  # Check file content
  cat("First few lines of the VTK file:\n")
  system(sprintf("head -n 10 %s", filename))
}

