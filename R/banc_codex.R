
#' Generate FlyWireCodex network visualization URLs for BANC v626 dataset
#'
#' Creates URLs for the FlyWireCodex interactive connectivity browser, allowing 
#' users to visualize neural network diagrams and explore connectivity patterns 
#' in the Brain And Nerve Cord (BANC) connectome. The BANC dataset represents 
#' the first unified brain-and-nerve-cord connectome of a limbed animal, revealing 
#' distributed control architecture and behavior-centric neural modules across 
#' the entire central nervous system.
#'
#' @param cell.types Character vector of cell type names to include in the network.
#'   Cell types represent functionally and morphologically distinct neuron classes
#'   in the BANC connectome (e.g., "DNa02" for specific descending neurons).
#' @param ids Character vector of neuron root IDs to include in the network.
#'   These are unique identifiers for individual neurons in the dataset.
#' @param codex.url Character string specifying the base FlyWireCodex URL.
#'   Defaults to the BANC connectivity browser.
#' @param open Logical indicating whether to automatically open the URL in 
#'   the default web browser. Default is FALSE.
#' @param min_syn_cnt Integer specifying the minimum number of synapses required
#'   for connections to be displayed. Default is 3.
#' @param edge_syn_cap Integer specifying the maximum number of synapses to 
#'   display per connection for visualization clarity. Default is 50.
#'
#' @return Character string containing the FlyWireCodex URL, or invisible NULL 
#'   if \code{open = TRUE}.
#'
#' @details
#' FlyWireCodex (\url{https://codex.flywire.ai/?dataset=banc}) provides an 
#' interactive web interface for exploring connectome data. This function 
#' generates properly formatted URLs that pre-configure the visualization 
#' with specified neurons or cell types, enabling researchers to quickly 
#' examine connectivity patterns, synaptic strengths, and network topology.
#' 
#' The BANC v626 dataset contains approximately 160,000 neurons with complete 
#' synaptic connectivity across both brain and ventral nerve cord, making it 
#' ideal for studying sensorimotor integration and distributed neural computation.
#'
#' @seealso \code{\link{banc_codex_search}} for generating search URLs,
#'   \code{\link{banc_edgelist}} for programmatic connectivity data access
#'
#' @export
#' @examples
#' \dontrun{
#' # Generate URL for DNa02 descending neurons and open in browser
#' banc_codex_network(cell.types = "DNa02", open = TRUE)
#' 
#' # Create URL for multiple cell types with custom synapse threshold
#' url <- banc_codex_network(cell.types = c("DNa02", "PFL3"), min_syn_cnt = 5)
#' 
#' # Generate URL for specific neuron IDs
#' banc_codex_network(ids = c("720575941566983162", "720575941562355975"))
#' }
banc_codex_network <- function(cell.types = NULL,
                               ids = NULL,
                               codex.url = "https://codex.flywire.ai/app/connectivity?dataset=banc",
                               open = FALSE,
                               min_syn_cnt = 3,
                               edge_syn_cap = 50){
  if(!is.null(cell.types)){
    cell.types.search <- paste(cell.types,collapse="+%7C%7C+cell_type+%3D%3D+")
    cell.types.search <- paste0("cell_type+%3D%3D+",gsub("\\+\\%7C\\%7C\\+cell_type\\+\\%3D\\%3D\\+$","",cell.types.search))
    cell.types.search <- paste(cell.types,collapse="+%7C%7C+cell_type+%3D%3D+")
    cell.types.search <- paste0("cell_type+%3D%3D+",gsub("\\+\\%7C\\%7C\\+cell_type\\+\\%3D\\%3D\\+$","",cell.types.search))
  }
  if(!is.null(ids)){
    ids.search <- paste(ids,collapse="+%7C%7C+root_id+%3D%3D+")
    ids.search <- paste0("root_id+%3D%3D+",gsub("\\+\\%7C\\%7C\\+root_id\\+\\%3D\\%3D\\+$","",ids.search))
    ids.search <- paste(ids,collapse="+%7C%7C+root_id+%3D%3D+")
    ids.search <- paste0("root_id+%3D%3D+",gsub("\\+\\%7C\\%7C\\+root_id\\+\\%3D\\%3D\\+$","",ids.search))
  }
  if(is.null(ids)&!is.null(cell.types)){
    search <- cell.types.search
  }else if(!is.null(ids)&is.null(cell.types)){
    search <- ids.search
  }else if(!is.null(ids)&!is.null(cell.types)){
    search <- paste(ids.search,"+%7C%7C+",cell.types.search)
  }else{
    stop("please provide argument cell.types or ids")
  }
  url <- sprintf("%s&cell_names_or_ids=%s&download=&group_by=type&edge_filter=all&cap=%d&min_syn_cnt=%d",
                 codex.url,
                 search,
                 edge_syn_cap,
                 min_syn_cnt)
  if(open){
    utils::browseURL(url)
    invisible()
  }else{
    url
  }
}

#' Generate FlyWireCodex search URLs for BANC v626 dataset
#'
#' Creates URLs for the FlyWireCodex search interface, allowing users to search
#' and browse neuron metadata in the Brain And Nerve Cord (BANC) connectome.
#' This function enables researchers to quickly access detailed information about
#' specific cell types or individual neurons, including morphological data,
#' cell type classifications, and connectivity summaries. The BANC dataset 
#' represents the first unified brain-and-nerve-cord connectome of a limbed 
#' animal, providing unprecedented insight into distributed neural control.
#'
#' @param cell.types Character vector of cell type names to search for.
#'   Cell types represent functionally and morphologically distinct neuron classes
#'   in the BANC connectome (e.g., "DNa02" for specific descending neurons).
#' @param ids Character vector of neuron root IDs to search for.
#'   These are unique identifiers for individual neurons in the dataset.
#' @param codex.url Character string specifying the base FlyWireCodex search URL.
#'   Defaults to the BANC search interface.
#' @param open Logical indicating whether to automatically open the URL in 
#'   the default web browser. Default is FALSE.
#' @param page.size Integer specifying the number of search results to display
#'   per page. Default is 100.
#'
#' @return Character string containing the FlyWireCodex search URL, or invisible 
#'   NULL if \code{open = TRUE}.
#'
#' @details
#' FlyWireCodex (\url{https://codex.flywire.ai/?dataset=banc}) provides an 
#' interactive web interface for exploring connectome data. The search function
#' allows researchers to query the comprehensive metadata associated with BANC
#' neurons, including cell type annotations, morphological measurements, and
#' connectivity statistics.
#' 
#' The BANC v626 dataset contains detailed annotations for approximately 160,000
#' neurons across brain and ventral nerve cord regions, with rich metadata 
#' supporting studies of sensorimotor integration, distributed computation, 
#' and behavior-centric neural modules.
#'
#' @seealso \code{\link{banc_codex_network}} for generating network visualization URLs,
#'   \code{\link{banc_codex_annotations}} for programmatic metadata access
#'
#' @export
#' @examples
#' \dontrun{
#' # Search for DNa02 descending neurons and open in browser
#' banc_codex_search(cell.types = "DNa02", open = TRUE)
#' 
#' # Create search URL for multiple cell types with custom page size
#' url <- banc_codex_search(cell.types = c("DNa02", "PFL3"), page.size = 50)
#' 
#' # Generate search URL for specific neuron IDs
#' banc_codex_search(ids = c("720575941566983162", "720575941562355975"))
#' }
banc_codex_search <- function(cell.types = NULL,
                              ids = NULL,
                              codex.url = "https://codex.flywire.ai/app/search?dataset=banc",
                              open = FALSE,
                              page.size = 100){
  if(!is.null(cell.types)){
    cell.types.search <- paste(cell.types,collapse="+%7C%7C+cell_type+%3D%3D+")
    cell.types.search <- paste0("cell_type+%3D%3D+",gsub("\\+\\%7C\\%7C\\+cell_type\\+\\%3D\\%3D\\+$","",cell.types.search))
    cell.types.search <- paste(cell.types,collapse="+%7C%7C+cell_type+%3D%3D+")
    cell.types.search <- paste0("cell_type+%3D%3D+",gsub("\\+\\%7C\\%7C\\+cell_type\\+\\%3D\\%3D\\+$","",cell.types.search))
  }
  if(!is.null(ids)){
    ids.search <- paste(ids,collapse="+%7C%7C+root_id+%3D%3D+")
    ids.search <- paste0("root_id+%3D%3D+",gsub("\\+\\%7C\\%7C\\+root_id\\+\\%3D\\%3D\\+$","",ids.search))
    ids.search <- paste(ids,collapse="+%7C%7C+root_id+%3D%3D+")
    ids.search <- paste0("root_id+%3D%3D+",gsub("\\+\\%7C\\%7C\\+root_id\\+\\%3D\\%3D\\+$","",ids.search))
  }
  if(is.null(ids)&!is.null(cell.types)){
    search <- cell.types.search
  }else if(!is.null(ids)&is.null(cell.types)){
    search <- ids.search
  }else if(!is.null(ids)&!is.null(cell.types)){
    search <- paste(ids.search,"+%7C%7C+",cell.types.search)
  }else{
    stop("please provide argument cell.types or ids")
  }
  url <- sprintf("%s&filter_string=%s&sort_by=&page_size=%d",
                 codex.url,
                 search,
                 page.size)
  if(open){
    utils::browseURL(url)
    invisible()
  }else{
    url
  }
}
