#' Download all of the BANC synapses as a .sqlite file that you can read lazily from later
#'
#' @details Downloads all automatic Zetta.ai synapse detections for the BANC and saves them as
#' a \code{banc_data.sqlite} file. Once this is done, in the future the function will read from this file
#' lazily so as not to throw the whole thing into system memory.
#'
#' @param path The google storage path to the desired synapses file. Read using \code{readr::read_csv}.
#' @param overwrite Logical, whether or not to overwrite an extant \code{banc_data.sqlite} file.
#' @param n_max Numeric, the maximum number of rows ot read from \code{path} if you just want to see
#' a taster of the file.
#' @param details Logical Whether or not to read all data columns in the target synapse \code{.csv}. Defaults to
#' \code{FALSE} in order to read only the essential presynapse position data.
#' @param min_size Numeric, filter parameter, the minimum size (in nm) of the detected synapse.
#' @param rawcoords Logical, whether or not yto convert from raw coordinates into nanometers. Default is `FALSE`.
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
banc_all_synapses <- function(path = "gs://zetta_lee_fly_cns_001_synapse/240623_run/assignment/final_edgelist.df",
                              overwrite = FALSE,
                              n_max = 2000,
                              details = FALSE,
                              min_size = 10,
                              rawcoords = FALSE){
  # Correct path to de-authenticate it, use https
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
    size = readr::col_number()
  )

  # Are we just sampling or going for the full thing?
  if(!is.null(n_max)){
    if(details){
      syns <- readr::read_csv(file=path, col_types = col.types, lazy = TRUE, n_max = n_max)
    }else{
      syns <- readr::read_csv(file=path, col_types = col.types, lazy = TRUE, n_max = n_max,
                              col_select = c(presyn_segid, presyn_x, presyn_y, presyn_z, size))
    }
    return(syns)
  }else if (!table_exists|overwrite){
    # Get all synapses
    if (details){
      syns <- readr::read_csv(file=path, col_types = col.types, lazy = TRUE)
    }else{
      syns <- readr::read_csv(file=path, col_types = col.types, lazy = TRUE,
                              col_select = c(presyn_segid, presyn_x, presyn_y, presyn_z, size))
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
