#' Create or refresh cache of BANC meta information
#'
#' @description
#' `banc_meta_create_cache()` builds or refreshes an in-memory cache of BANC metadata
#' for efficient repeated lookups. You can choose the data source using `use_seatable`.
#' The main accessor function [banc_meta()] will always use the most recently created cache.
#'
#' @details
#' BANC meta queries can be slow; caching avoids repeated computation/database access.
#' Whenever labels are updated, simply rerun this function to update the cache.
#'
#' @param use_seatable Whether to build BANC meta data from the `codex_annotations` CAVE table
#' (production) or our internal seatable (development). Both require different types of authenticated
#' access, for details see `bancr` documentation.
#' @param return Logical; if `TRUE`, return the cache tibble/invisible.
#'
#' @return Invisibly returns the cache (data.frame) if `return=TRUE`; otherwise invisibly `NULL`.
#' @export
#'
#' @examples
#' \dontrun{
#' #' # Requires authenticated access to BANC CAVE
#' banc_meta_cache(use_seatable=FALSE)
#'
#' banc_meta_create_cache(use_seatable=TRUE) # create cache
#' ## BANCTABLE_TOKEN must be set, see bancr package
#' result <- banc_meta() # use cache
#'
#' # use cache to quickly make plot
#' register_banc_coconat()
#' cf_cosine_plot(cf_ids('/type:LAL0(08|09|10|42)', datasets = c("banc", "hemibrain")))
#' }
banc_meta_create_cache <- NULL # Placeholder, assigned below

#' Query cached BANC meta data
#'
#' @description
#' Returns results from the in-memory cache, filtered by `ids` if given.
#' Cache must be created first using [banc_meta_create_cache()].
#'
#' @details
#' `banc_meta()` never queries databases directly.
#' If `ids` are given, filters the meta table by root_id.
#'
#' @param ids Vector of neuron/root IDs to select, or `NULL` for all.
#' @return tibble/data.frame, possibly filtered by ids.
#' @export
#' @seealso [banc_meta_create_cache()]
#'
#' @examples
#' \dontrun{
#' banc_meta_create_cache() # build the cache
#' all_meta <- banc_meta()  # retrieve all
#' }
banc_meta <- NULL # Placeholder, assigned below

# hidden
banc_meta <- local({
  .banc_meta_cache <- NULL

  .refresh_cache <- function(use_seatable=FALSE) {
    if (use_seatable) {
      # Read from seatable
      banc.meta <- banctable_query(
        "SELECT root_id, side, cell_type, cell_class, cell_sub_class from banc_meta"
      )
      banc.meta %>%
        dplyr::rename(
          id = root_id,
          class = cell_class,
          type = cell_type,
          side = side,
          subclass = cell_sub_class
        ) %>%
        dplyr::mutate(id = as.character(id))
    } else {
      banc.community.meta <- banc_cell_info() %>%
        dplyr::filter(valid == 't') %>%
        dplyr::arrange(pt_root_id, tag) %>%
        dplyr::distinct(pt_root_id, tag2, tag, .keep_all = TRUE) %>%
        dplyr::group_by(pt_root_id, tag2) %>%
        dplyr::summarise(
          tag = {
            if (length(tag) > 1 && any(grepl("?", tag, fixed = TRUE))) {
              usx = unique(sub("?", "", tag, fixed = TRUE))
              if (length(usx) < length(tag)) tag = usx
            }
            paste0(tag, collapse = ";")
          },
          .groups = 'drop'
        ) %>%
        tidyr::pivot_wider(
          id_cols = pt_root_id,
          names_from = tag2,
          values_from = tag,
          values_fill = ""
        ) %>%
        dplyr::select(
          id = pt_root_id,
          class = `primary class`,
          type = `neuron identity`,
          side = `soma side`,
          subclass = `anterior-posterior projection pattern`
        ) %>%
        dplyr::mutate(class = gsub(" ","_", class))

      banc.codex.meta <- banc_codex_annotations() %>%
        dplyr::distinct(pt_root_id, .keep_all = TRUE) %>%
        dplyr::select(
          id = pt_root_id,
          class = cell_class,
          type = cell_type,
          side = side,
          subclass = cell_sub_class
        )

      rbind(
        banc.codex.meta,
        banc.community.meta
      ) %>%
        dplyr::distinct(id, .keep_all = TRUE) %>%
        dplyr::mutate(id = as.character(id))
    }
  }

  list(
    create_cache = function(use_seatable=FALSE, return = FALSE) {
      meta <- .refresh_cache(use_seatable=use_seatable)
      .banc_meta_cache <<- meta
      if (return) meta else invisible()
    },
    get_meta = function(ids = NULL) {
      if (is.null(.banc_meta_cache)){
        message("No BANC meta cache loaded. Creating with banc_meta_create_cache(use_seatable=FALSE)")
        banc_meta_create_cache(use_seatable=FALSE)
      }
      meta <- .banc_meta_cache
      if (!is.null(ids)) {
        ids <- extract_ids(unname(unlist(ids)))
        ids <- tryCatch(banc_ids(ids), error = function(e) NULL)
        meta %>% dplyr::filter(id %in% ids)
      } else {
        meta
      }
    }
  )
})

# Exported user-friendly functions
banc_meta_create_cache <- banc_meta$create_cache
banc_meta <- banc_meta$get_meta

# banc_coconat.R
coconat_banc_meta <- function(ids) {
  banc_meta(ids)
}

# hidden
coconat_banc_ids <- function(ids=NULL) {
  if(is.null(ids)) return(NULL)
  # extract numeric ids if possible
  ids <- extract_ids(ids)
  if(is.character(ids) && length(ids)==1 && !fafbseg:::valid_id(ids)) {
    # query
    metadf=banc_meta()
    if(isTRUE(ids=='all')) return(banc_ids(metadf$id, integer64 = F))
    if(isTRUE(ids=='neurons')) {
      ids <- metadf %>%
        dplyr::filter(is.na(class) | class!='glia') %>%
        dplyr::pull(id)
      return(banc_ids(ids, integer64 = F))
    }
    if(isTRUE(substr(ids, 1, 1)=="/"))
      ids=substr(ids, 2, nchar(ids))
    else warning("All BANC queries are regex queries. ",
                 "Use an initial / to suppress this warning!")
    if(!grepl(":", ids)) ids=paste0("type:", ids)
    qsplit=stringr::str_match(ids, pattern = '[/]{0,1}(.+):(.+)')
    field=qsplit[,2]
    value=qsplit[,3]
    if(!field %in% colnames(metadf)) {
      stop(glue("BANC queries only work with these fields: ",
                paste(colnames(metadf)[-1], collapse = ',')))
    }
    ids <- metadf %>%
      dplyr::filter(grepl(value, .data[[field]])) %>%
      dplyr::pull(id)
  } else if(length(ids)>0) {
    # check they are valid for current materialisation
    banc_latestid(ids, version = banc_version())
  }
  return(banc_ids(ids, integer64 = F))
}

# minimal version of this function
coconat_banc_partners <- function(ids,
                                        partners,
                                        threshold,
                                        version=banc_version(),
                                        ...) {
  tres=banc_partner_summary(banc_ids(ids),
                                   partners = partners,
                                   threshold = threshold-1L,
                                   version=version,
                                   ...)
  tres$side=substr(toupper(tres$side),1,1)
  partner_col=grep("_id", colnames(tres), value = T)
  for(pc in partner_col){
    tres[[pc]] <- as.character(tres[[pc]])
  }
  # nb coconatfly can looks after adding metadata
  tres
}

#' Use BANC data with coconat for connectivity similarity
#'
#' @description
#' Register the BANC dataset with \href{https://github.com/natverse/coconat}{coconat},
#' a natverse R package for between and within dataset connectivity comparisons using cosine similarity.
#'
#' @details
#' `register_banc_coconat()` registers `bancr`-backed functionality for use with
#'
#' @param showerror Logically, error-out silently or not.
#' @export
#' @seealso [banc_meta_create_cache()]
#'
#' @examples
#' \dontrun{
#' library(coconat)
#' banc_meta_create_cache(use_seatable=TRUE)
#' register_banc_coconat()
#' cf_cosine_plot(cf_ids('/type:LAL0(08|09|10|42)', datasets = c("banc", "hemibrain","flywire")))
#' }
register_banc_coconat <- function(showerror=TRUE){
  if (!requireNamespace("coconat", quietly = TRUE)) {
    stop("Package 'coconat' is required for this function. Please install it with: devtools::install_github(natverse/coconat)")
  }
  if(requireNamespace('coconat', quietly = !showerror))
    coconat::register_dataset(
      name = 'banc',
      shortname = 'bc',
      namespace = 'coconatfly',
      metafun = coconat_banc_meta,
      idfun = coconat_banc_ids,
      partnerfun = coconat_banc_partners
    )
}
