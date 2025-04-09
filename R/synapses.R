#' Download all of the BANC synapses as a .sqlite file that you can read lazily from later
#'
#' @details Downloads all automatic Zetta.ai synapse detections for the BANC and saves them as
#' a \code{banc_data.sqlite} file. Once this is done, in the future the function will read from this file
#' lazily so as not to throw the whole thing into system memory.
#'
#' @param path The google storage path to the desired synapses file. Read using \code{readr::read_csv}.
#' @param overwrite Logical, whether or not to overwrite an extant \code{banc_data.sqlite} file.
#' @param n_max Numeric, the maximum number of rows to read from \code{path} if you just want to see
#' a taster of the file.
#' @param details Logical Whether or not to read all data columns in the target synapse \code{.csv}. Defaults to
#' \code{FALSE} in order to read only the essential presynapse position data.
#' @param min_size Numeric, filter parameter, the minimum size (in nm) of the detected synapse.
#' @param rawcoords Logical, whether or not to convert from raw coordinates into nanometers. Default is `FALSE`.
#'
#' @return a data.frame
#'
#' @seealso \code{\link{banc_partner_summary}}, \code{\link{banc_partners}}
#' @export
#'
#' @examples
#' \dontrun{
#' syns <- banc_all_synapses()
#' }
banc_all_synapses <- function(path = c("gs://lee-lab_brain-and-nerve-cord-fly-connectome/synapses/v2.0/final_edgelist.csv",
                                       "gs://zetta_lee_fly_cns_001_synapse/240623_run/assignment/final_edgelist.df"),
                              overwrite = FALSE,
                              n_max = 2000,
                              details = FALSE,
                              min_size = 10,
                              rawcoords = FALSE){
  # Correct path to de-authenticate it, use https
  path <- match.arg(path)
  path <- gsub("^gs\\:\\/","https://storage.googleapis.com",path)

  # Check if the file exists, if not, create it
  file_path <- file.path(system.file(package = "bancr"), "data", "banc_data.sqlite")
  if (!file.exists(file_path)) {
    con <- DBI::dbConnect(RSQLite::SQLite(), file_path)
    DBI::dbDisconnect(con)
    table_exists <- FALSE
    message("Created: ", file_path)
  }else{
    # Check if the 'synapses' table exists
    con <- DBI::dbConnect(RSQLite::SQLite(), file_path)
    table_exists <- "synapses" %in% DBI::dbListTables(con)
  }

  # The column types for our read
  col.types <- readr::cols(
    .default = readr::col_character(),
    cleft_segid  = readr::col_character(),
    centroid_x = readr::col_number(),
    centroid_y = readr::col_number(),
    centroid_z = readr::col_number(),
    bbox_bx = readr::col_number(),
    bbox_by = readr::col_number(),
    bbox_bz = readr::col_number(),
    bbox_ex = readr::col_number(),
    bbox_ey = readr::col_number(),
    bbox_ez = readr::col_number(),
    presyn_segid = readr::col_character(),
    postsyn_segid  = readr::col_character(),
    presyn_x = readr::col_integer(),
    presyn_y = readr::col_integer(),
    presyn_z = readr::col_integer(),
    postsyn_x = readr::col_integer(),
    postsyn_y = readr::col_integer(),
    postsyn_z = readr::col_integer(),
    clefthash = readr::col_number(),
    partnerhash = readr::col_integer(),
    size = readr::col_number(),
    clefthash = readr::col_number(),
    partnerhash = readr::col_number()
  )

  # Are we just sampling or going for the full thing?
  if(!is.null(n_max)){
    if(details){
      syns <- readr::read_csv(file=path, col_types = col.types, lazy = TRUE, n_max = n_max)
    }else{
      syns <- readr::read_csv(file=path, col_types = col.types, lazy = TRUE, n_max = n_max,
                              col_select = c("presyn_segid", "presyn_x", "presyn_y", "presyn_z", "size"))
    }
    return(syns)
  }else if (!table_exists|overwrite){
    # Get all synapses
    if (details){
      syns <- readr::read_csv(file=path, col_types = col.types, lazy = TRUE)
    }else{
      syns <- readr::read_csv(file=path, col_types = col.types, lazy = TRUE,
                              col_select = c("presyn_segid", "presyn_x", "presyn_y", "presyn_z", "size"))
    }

    # # Process
    # if(!is.null(min_size)){
    #   syns <- syns %>%
    #     dplyr::filter(size>=min_size)
    # }
    # if(!rawcoords){
    #   syns[,c("presyn_x", "presyn_y", "presyn_z")] <- bancr::banc_raw2nm(syns[,c("presyn_x", "presyn_y", "presyn_z")])
    # }

    # Connect to the SQLite database
    con <- DBI::dbConnect(RSQLite::SQLite(), file_path)

    # Write the data frame to the 'synapses' table
    # If the table already exists, it will be overwritten
    DBI::dbWriteTable(con, "synapses", syns, overwrite = TRUE)
    DBI::dbDisconnect(con)
    message("Added tab synapses, no. rows: ", nrow(syns))
  }

  # Read
  dplyr::tbl(src = con, from = "synapses")

}
# Helpful scene: https://spelunker.cave-explorer.org/#!middleauth+https://global.daf-apis.com/nglstate/api/v1/4753860997414912

# # googleCloudStorageR
# banc_gcs_read <- function(path = "gs://zetta_lee_fly_cns_001_synapse/240623_run/assignment/final_edgelist.df"){

  # # Not sure hot to get this working in R, but looks useful
  # googleCloudStorageR::gcs_setup(token="none")
  # googleCloudStorageR::gcs_auth(token=NULL)
  # # OR?
  # scope <-c("https://www.googleapis.com/auth/cloud-platform")
  # token <- gargle::token_fetch(scopes = scope)
  # googleCloudStorageR::gcs_auth(token = token)
  # # Then get data?
  # path <- 'gs://zetta_lee_fly_cns_001_synapse/240529_run/240604_assignment/final_edgelist.df'
  # df <- googleCloudStorageR::gcs_get_object(path, parseFunction = function(x) read.csv(x, nrows = 1000))
# }


#' Add synapses to neuron objects
#'
#' This function family adds synaptic data to neuron objects or neuron lists.
#' It retrieves synaptic connections and attaches them to the neuron object(s).
#'
#' @param x A neuron object, neuronlist, or other object to add synapses to
#' @param id The root ID of the neuron. If `NULL`, it uses the ID from the neuron object
#' @param connectors A dataframe of synaptic connections. If `NULL`, it retrieves the data
#' @param size.threshold Minimum size threshold for synapses to include
#' @param remove.autapses Whether to remove autapses (self-connections)
#' @param update.id Logical, whether or not to use \code{banc_latestid} to update the neuron's `root_id` when fetching synapses.
#' @param ... Additional arguments passed to methods, \code{nat::nlapply}
#'
#' @return An object of the same type as `x`, with synapses added
#' @examples
#' \dontrun{
#' # Get BANC ID for DNA01
#' id <- "720575941572711675"
#' id <- banc_latestid(id)
#'
#' # Get the L2 skeletons
#' n <- banc_read_l2skel(id)
#'
#' # Re-root to soma
#' n.rerooted <- banc_reroot(n)
#'
#' # Add synapse information, stored at n.syn[[1]]$connectors
#' n.syn <- banc_add_synapses(n.rerooted)
#'
#' # Split neuron
#' n.split <- hemibrainr::flow_centrality(n.syn)
#'
#' # Visualise
#' banc_neuron_comparison_plot(n.split)
#' }
#' @export
banc_add_synapses <- function(x,
                              id = NULL,
                              connectors = NULL,
                              size.threshold = 5,
                              remove.autapses = TRUE,
                              update.id = TRUE,
                              ...) {
  UseMethod("banc_add_synapses")
}

#' @rdname banc_add_synapses
#' @export
banc_add_synapses.neuron <- function(x,
                              id = NULL,
                              connectors = NULL,
                              size.threshold = 5,
                              remove.autapses = TRUE,
                              update.id = TRUE,
                              ...){
  # Get valid root id
  if(is.null(id)){
    id <- x$id
  }
  if(update.id){
    id <- banc_latestid(id)
  }

  # Get synaptic data
  if(is.null(connectors)){
    connectors.in <- banc_partners(id, partners = "input")
    if(nrow(connectors.in)){
      connectors.in.xyz <- do.call(rbind,connectors.in$post_pt_position)
      connectors.in.xyz <- as.data.frame(connectors.in.xyz)
      colnames(connectors.in.xyz) <- c("X","Y","Z")
      connectors.in <- cbind(connectors.in,connectors.in.xyz)
      connectors.in <- connectors.in %>%
        dplyr::rename(connector_id = .data$id,
                      pre_id = .data$pre_pt_root_id,
                      pre_svid = .data$pre_pt_supervoxel_id,
                      post_id = .data$post_pt_root_id,
                      post_svid = .data$post_pt_supervoxel_id) %>%
        dplyr::filter(size>size.threshold) %>%
        dplyr::mutate(prepost = 1) %>%
        dplyr::select(.data$connector_id,
                      .data$pre_id, .data$post_id, .data$prepost,
                      .data$pre_svid, .data$post_svid, .data$size,
                      .data$X, .data$Y, .data$Z)
    }
    connectors.out <- banc_partners(id, partners = "output")
    if(nrow(connectors.out)){
      connectors.out.xyz <- do.call(rbind,connectors.out$pre_pt_position)
      connectors.out.xyz <- as.data.frame(connectors.out.xyz)
      colnames(connectors.out.xyz) <- c("X","Y","Z")
      connectors.out <- cbind(connectors.out,connectors.out.xyz)
      connectors.out <- connectors.out %>%
        dplyr::rename(connector_id = .data$id,
                      pre_id = .data$pre_pt_root_id,
                      pre_svid = .data$pre_pt_supervoxel_id,
                      post_id = .data$post_pt_root_id,
                      post_svid = .data$post_pt_supervoxel_id) %>%
        dplyr::filter(.data$size>size.threshold) %>%
        dplyr::mutate(prepost = 0) %>%
        dplyr::select(.data$connector_id,
                      .data$pre_id, .data$post_id, .data$prepost,
                      .data$pre_svid, .data$post_svid, .data$size,
                      .data$X, .data$Y, .data$Z)
    }
    connectors <- rbind(connectors.in,connectors.out)
  }else{
    connectors <- connectors %>%
      dplyr::filter(.data$post_id==id|.data$pre_id==id)
  }

  # Attach synapses
  if(nrow(connectors)){
    if(remove.autapses) {
      connectors=connectors[connectors$post_id!=connectors$pre_id,,drop=FALSE]
    }
    near <- nabor::knn(query = nat::xyzmatrix(connectors),
                       data = nat::xyzmatrix(x$d),k=1)
    connectors$treenode_id <- x$d[near$nn.idx,"PointNo"]
    x$connectors = as.data.frame(connectors, stringsAsFactors = FALSE)
  }else{
    connectors <- data.frame()
  }
  x$connectors <- connectors

  # Change class to work with connectivity functions in other packages
  class(x) <- union(c("synapticneuron"), class(x))

  # Return
  x
}

#' @rdname banc_add_synapses
#' @export
banc_add_synapses.neuronlist <- function(x,
                                         id = NULL,
                                         connectors = NULL,
                                         size.threshold = 5,
                                         remove.autapses = TRUE,
                                         update.id = TRUE,
                                         ...) {
  if(is.null(id)){
    x <- add_field_seq(x, entries= names(x), field = "id")
  }
  nat::nlapply(x,
               banc_add_synapses.neuron,
               id=NULL,
               connectors=connectors,
               size.threshold=size.threshold,
               remove.autapses=remove.autapses,
               ...)
}

#' @rdname banc_add_synapses
#' @export
banc_add_synapses.default <- function(x,
                                      id = NULL,
                                      connectors = NULL,
                                      size.threshold = 5,
                                      remove.autapses = TRUE,
                                      update.id = TRUE,
                                      ...) {
  stop("No method for class ", class(x))
}










