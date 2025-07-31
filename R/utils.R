# hidden
check_package_available <- function(pkg, repo=c("CRAN", "Bioconductor")) {
  if(!requireNamespace(pkg, quietly = TRUE)) {
    repo=match.arg(repo)
    installmsg=switch(repo,
                      CRAN=paste0("install.packages('",pkg,"')"),
                      Bioconductor=paste0('if (!require("BiocManager", quietly = TRUE))',
                                          '
  install.packages("BiocManager")',
                                          '
BiocManager::install("',pkg,'")'))
    stop("Please install suggested package: ", pkg, " by doing
",
         installmsg, call. = F)
  }
}

# hidden
euclidean_distances <- function(A, B) {
  sqrt(rowSums((A - B)^2))
}

# Helper
express_lane <- function(base_dir, search = "^1_", link = c("move","symlink","copy")) {
  link <- match.arg(link)
  todo_dir <- base_dir #fs::path(base_dir, "todo")
  if (fs::dir_exists(todo_dir)) {
    express_dir <- fs::path(base_dir, "express")

    # Create express directory if it doesn't exist
    fs::dir_create(express_dir)
    express.files <- list.files(express_dir, full.names = TRUE)
    remove.old <- file.remove(express.files)

    # Get all PNG files in todo directory that start with '1'
    png_files <- fs::dir_ls(todo_dir, recurse = TRUE, glob = "*.png")
    png_files<- png_files[grepl(search, fs::path_file(png_files))]

    # Create symlinks
    purrr::walk(png_files, function(file) {
      # Use only the filename for the symlink, not the full path
      file_name <- fs::path_file(file)
      symlink_path <- fs::path(express_dir, file_name)

      # If a symlink with this name already exists, add a numeric suffix
      if (fs::file_exists(symlink_path)) {
        i <- 1
        while (fs::file_exists(symlink_path)) {
          symlink_path <- fs::path(express_dir, paste0(fs::path_ext_remove(file_name), "_", i, ".png"))
          i <- i + 1
        }
      }

      # Create symlink
      if(link=="symlink"){
        fs::link_create(file, symlink_path)
      }else if(link=="move"){
        fs::file_move(file,symlink_path)
      }else{
        fs::file_copy(file,symlink_path)
      }
      cat("Created symlink:", symlink_path, "->", file, "
")
    })
    remove_empty_dirs(base_dir)
    cat("Symlink creation completed.
")
  }else{
    cat("No 'todo' folder from which to connect
")
  }
  invisible()
}

# helper
remove_empty_dirs <- function(base_dir) {
  dirs <- fs::dir_ls(base_dir, recurse = TRUE, type = "directory") %>%
    purrr::discard(~grepl("/express$|/done$", .x))
  dirs <- rev(dirs)  # Start from the deepest directories
  purrr::walk(dirs, function(dir) {
    if (length(fs::dir_ls(dir)) == 0) {
      fs::dir_delete(dir)
      cat("Removed empty directory:", dir, "
")
    }
  })
}


