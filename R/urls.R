#' Return a sample Neuroglancer scene URL for banc dataset
#'
#' @details See
#'   \href{https://banc-reconstruction.slack.com/archives/C01RZP5JH9C/p1616522511001900}{banc
#'    slack} for details.
#'
#' @param url a spelunker neuroglancer URL.
#' @param ids A set of root ids to include in the scene. Can also be a
#'   data.frame.
#' @param layer the segmentation layer for which `ids` intended. Defaults to 'segmentation proofreading',
#' but could point to another dataset layer.
#' @param open Whether to open the URL in your browser (see
#'   \code{\link{browseURL}})
#' @return A character vector containing a single Neuroglancer URL (invisibly
#'   when \code{open=TRUE}).
#' @seealso \code{\link{bancsee}}
#' @export
#' @importFrom utils browseURL
#' @examples
#' \dontrun{
#' browseURL(banc_scene())
#' banc_scene(open=T)
#' banc_scene("720575941545083784", open=T)
#' }
banc_scene <- function(ids=NULL,
                       open=FALSE,
                       layer = NULL,
                       url="https://spelunker.cave-explorer.org/#!middleauth+https://global.daf-apis.com/nglstate/api/v1/6431332029693952") {
  #url="https://spelunker.cave-explorer.org/#!middleauth+https://global.daf-apis.com/nglstate/api/v1/6283844278812672"
  #url="https://spelunker.cave-explorer.org/#!middleauth+https://global.daf-apis.com/nglstate/api/v1/6431332029693952"
  url=sub("#!middleauth+", "?", url, fixed = T)
  parts=unlist(strsplit(url, "?", fixed = T))
  json=try(fafbseg::flywire_fetch(parts[2], token=banc_token(), return = 'text', cache = TRUE))
  if(inherits(json, 'try-error')) {
    badtoken=paste0("You have a token but it doesn't seem to be authorised for banc.\n",
                    "Have you definitely used `banc_set_token()` to make a token for the banc dataset?")
    if(grepl(500, json))
      stop("There seems to be a (temporary?) problem with the zetta server!")
    else if(grepl(401, json))
      stop(badtoken)

    token=try(banc_token(), silent = T)
    if(inherits(token, 'try-error'))
      stop("It looks like you do not have a stored token. Please use `banc_set_token()` to make one.")
    else
      stop(badtoken)
  }
  u=ngl_encode_url(json, baseurl = parts[1])
  if(!is.null(ids)){
    banc_ngl_segments(u, layer=layer) <- banc_ids(ids)
  }
  if(open) {
    browseURL(u)
    invisible(u)
  } else (u)
}

#' Visualise neurons across multiple Drosophila connectomic datasets in BANC spelunker
#'
#' @description
#' This function constructs a Neuroglancer scene that visualizes neurons from multiple
#' co-registered Drosophila connectomic datasets, including BANC, FAFB, hemibrain, and MANC.
#' It allows for simultaneous visualization of corresponding neurons across these datasets.
#'
#' @param url a spelunker neuroglancer URL.
#' @param banc_ids A vector of neuron IDs from the BANC dataset. Default is NULL.
#' @param fafb_ids A vector of neuron IDs from the FAFB dataset. Default is NULL.
#' @param hemibrain_ids A vector of neuron IDs from the hemibrain dataset. Default is NULL.
#' @param manc_ids A vector of neuron IDs from the MANC dataset. Default is NULL.
#' @param nuclei_ids A vector of nuclei IDs for the BANC dataset. Default is NULL.
#' @param open Logical; if TRUE, the function will open the Neuroglancer scene in a web browser. Default is FALSE.
#' @param banc.cols Vector of hex codes describing a colour spectrum of colors to be interpolated for BANC neurons. Defaults are cyan-purple.
#' @param fafb.cols Vector of hex codes describing a colour spectrum of colors to be interpolated for BANC neurons. Defaults are red hues.
#' @param hemibrain.cols Vector of hex codes describing a colour spectrum of colors to be interpolated for BANC neurons. Defaults  green hues.
#' @param hemibrain.mirrored.cols Vector of hex codes describing a colour spectrum of colors to be interpolated for BANC neurons. Defaults are yellow hues.
#' @param manc.cols Vector of hex codes describing a colour spectrum of colors to be interpolated for MANC neurons. Defaults are orange hues.
#' @param nulcei.col Hex code for the colour in which nuclei will be plotted. Default is pink.
#'
#' @return
#' If `open = FALSE`, returns a character string containing the URL for the Neuroglancer scene.
#' If `open = TRUE`, opens the Spelunker Neuroglancer scene in a web browser and invisibly returns the URL.
#'
#' @details
#' The function creates a Neuroglancer scene with multiple layers, each corresponding to a different dataset:
#' - BANC: "segmentation proofreading" layer
#' - FAFB: "fafb v783 imported" layer
#' - Hemibrain: "hemibrain v1.2.1 imported" and "hemibrain v1.2.1 imported, mirrored" layers
#' - MANC: "manc v1.2.1 imported" layer
#' - BANC nuclei: "nuclei (v1)" layer
#'
#' Each dataset is assigned a unique color palette to distinguish neurons from different sources:
#' - BANC: Blue to purple spectrum
#' - FAFB: Red spectrum
#' - Hemibrain: Green spectrum (original) and Yellow spectrum (mirrored)
#' - MANC: Orange spectrum
#' - BANC nuclei: Pink
#'
#' @note
#' This function suppresses all warnings during execution. While this ensures smooth operation,
#' it may hide important messages. Use with caution and refer to individual function documentation
#' if unexpected behavior occurs.
#'
#' @examples
#' \dontrun{
#' # Visualize cell type DNa01 across datasets
#' bancsee(banc_ids = c("720575941493078142","720575941455137261"),
#'         fafb_ids = c("720575940644438551","720575940627787609"),
#'         hemibrain_ids = c("1170939344"),
#'         manc_ids = c("10751","10760"),
#'         open = TRUE)
#'
#' # Get URL without opening browser
#' url <- bancsee(banc_ids = c("720575941493078142"),
#'                fafb_ids = c("720575940644438551"),
#'                open = FALSE)
#' }
#'
#' @importFrom grDevices colorRampPalette
#' @importFrom utils browseURL
#' @seealso \code{\link{banc_scene}}
#'
#' @export
bancsee <- function(banc_ids = NULL,
                    fafb_ids = NULL,
                    hemibrain_ids = NULL,
                    manc_ids = NULL,
                    nuclei_ids = NULL,
                    open = FALSE,
                    banc.cols = c("#54BCD1", "#0000FF", "#8A2BE2"),
                    fafb.cols = c("#C41E3A", "#FF3131", "#F88379"),
                    hemibrain.cols = c("#00FF00", "#32CD32", "#006400"),
                    hemibrain.mirrored.cols = c("#FFFF00", "#FFD700", "#FFA500"),
                    manc.cols = c("#FFA07A", "#FF4500", "#FF8C00"),
                    nulcei.col = "#FC6882",
                    url="https://spelunker.cave-explorer.org/#!middleauth+https://global.daf-apis.com/nglstate/api/v1/6431332029693952"){

  # Do not get the neuroglancer warnings
  old_warn <- options(warn = -1)  # Suppress all warnings
  on.exit(options(old_warn))  # Restore warning settings when function exits
  banc_ids <- banc_ids[!is.na(banc_ids)]
  banc_ids <- banc_ids[!banc_ids%in%c("0","","NA")]
  fafb_ids <- fafb_ids[!is.na(fafb_ids)]
  fafb_ids <- fafb_ids[!fafb_ids%in%c("0","","NA")]
  hemibrain_ids <- hemibrain_ids[!is.na(hemibrain_ids)]
  hemibrain_ids <- hemibrain_ids[!hemibrain_ids%in%c("0","","NA")]
  manc_ids <- manc_ids[!is.na(manc_ids)]
  manc_ids <- manc_ids[!manc_ids%in%c("0","","NA")]
  nuclei_ids <- nuclei_ids[!is.na(nuclei_ids)]
  nuclei_ids <- nuclei_ids[!nuclei_ids%in%c("0","","NA")]

  # Get BANC IDs
  if(length(banc_ids)){
    u1=banc_scene(url = url, ids = banc_ids, open=F, layer = "segmentation proofreading")
    colourdf1 = data.frame(ids = banc_ids,
                           col=grDevices::colorRampPalette(banc.cols)(length(banc_ids)))
    sc1<-fafbseg::ngl_add_colours(u1, colourdf1, layer = "segmentation proofreading")
  }else{
    sc1 = fafbseg::ngl_decode_scene(banc_scene(url = url))
    banc_ngl_segments(sc1) <- NULL
  }

  if(length(fafb_ids)){
    u2=banc_scene(url = url, ids = fafb_ids, open=F, layer = "fafb v783 imported")
    colourdf2 = data.frame(ids = fafb_ids,
                           col=grDevices::colorRampPalette(fafb.cols)(length(fafb_ids)))
    sc2<-fafbseg::ngl_add_colours(u2, colourdf2, layer = "fafb v783 imported")
    fafbseg::ngl_layers(sc1)$`fafb v783 imported` <- fafbseg::ngl_layers(sc2)$`fafb v783 imported`
  }

  if(length(hemibrain_ids)){
    u3=banc_scene(url = url, ids = hemibrain_ids, open=F, layer = "hemibrain v1.2.1 imported")
    colourdf3 = data.frame(ids = hemibrain_ids,
                           col=grDevices::colorRampPalette(hemibrain.cols)(length(hemibrain_ids)))
    sc3<-fafbseg::ngl_add_colours(u3, colourdf3, layer = "hemibrain v1.2.1 imported")
    fafbseg::ngl_layers(sc1)$`hemibrain v1.2.1 imported` <- fafbseg::ngl_layers(sc3)$`hemibrain v1.2.1 imported`
    u4=banc_scene(url = url, ids = hemibrain_ids, open=F, layer = "hemibrain v1.2.1 imported, mirrored")
    colourdf4 = data.frame(ids = hemibrain_ids,
                           col=grDevices::colorRampPalette(hemibrain.mirrored.cols)(length(hemibrain_ids)))
    sc4<-fafbseg::ngl_add_colours(u4, colourdf4, layer = "hemibrain v1.2.1 imported, mirrored")
    fafbseg::ngl_layers(sc1)$`hemibrain v1.2.1 imported, mirrored` <- fafbseg::ngl_layers(sc4)$`hemibrain v1.2.1 imported, mirrored`
  }

  if(length(manc_ids)){
    u5=banc_scene(url = url, ids = manc_ids, open=F, layer = "manc v1.2.1 imported")
    colourdf5 = data.frame(ids = manc_ids,
                           col=grDevices::colorRampPalette(manc.cols)(length(manc_ids)))
    sc5<-fafbseg::ngl_add_colours(u5, colourdf5, layer = "manc v1.2.1 imported")
    fafbseg::ngl_layers(sc1)$`manc v1.2.1 imported` <- fafbseg::ngl_layers(sc5)$`manc v1.2.1 imported`
  }

  if(length(nuclei_ids)){
    u6=banc_scene(url = url, ids = nuclei_ids, open=F, layer = "nuclei (v1)")
    colourdf6 = data.frame(ids = nuclei_ids,
                           col=nulcei.col)
    sc6<-fafbseg::ngl_add_colours(u6, colourdf6, layer = "nuclei (v1)")
    fafbseg::ngl_layers(sc1)$`nuclei (v1)` <- fafbseg::ngl_layers(sc6)$`nuclei (v1)`
  }

  # see
  u<-as.character(sc1)
  if(open) {
    utils::browseURL(u)
    invisible(u)
  } else {
    banc_shorturl(u)
  }
}

# hidden
banc_shorturl <- function (x,
                           baseurl = NULL,
                           cache = TRUE,
                           ...)
{
  if (fafbseg:::is.ngscene(x)) {
    sc <- x
    x <- fafbseg::ngl_encode_url(sc, ...)
  }
  else {
    stopifnot(is.character(x))
    if (length(x) > 1) {
      res = pbapply::pbsapply(x,
                              banc_shorturl,
                              baseurl = baseurl,
                              cache = cache,
                              ...)
      return(res)
    }
    sc = fafbseg::ngl_decode_scene(x)
  }
  state_server = "https://global.daf-apis.com/nglstate/post"
  json = fafbseg::ngl_decode_scene(x, return.json = TRUE)
  res = fafbseg::flywire_fetch(state_server, body = json, cache = cache)
  sprintf("https://spelunker.cave-explorer.org/#!middleauth+https://global.daf-apis.com/nglstate/api/v1/%s",basename(res))
}

#' Choose or (temporarily) use the banc autosegmentation
#'
#' @details \code{bancr} inherits a significant amount of infrastructure from
#'   the \code{\link{fafbseg}} package. This has the concept of the
#'   \emph{active} autosegmentation, which in turn defines one or more R options
#'   containing URIs pointing to voxel-wise segmentation, mesh etc data. These
#'   are normally contained within a single neuroglancer URL which points to
#'   multiple data layers. For banc this is the neuroglancer scene returned by
#'   \code{\link{banc_scene}}.
#' @param set Whether or not to permanently set the banc autosegmentation as the
#'   default for \code{\link{fafbseg}} functions.

#'
#' @return If \code{set=TRUE} a list containing the previous values of the
#'   relevant global options (in the style of \code{\link{options}}. If
#'   \code{set=FALSE} a named list containing the option values.
#' @export
#'
#' @examples
#' \dontrun{
#' choose_banc()
#' options()[grep("^fafbseg.*url", names(options()))]
#' }
choose_banc <- function(set=TRUE) {
  fafbseg::choose_segmentation(
    banc_scene(),
    set=set,
    moreoptions=list(
      fafbseg.cave.datastack_name=banc_datastack_name()
      ))
}

#' @param expr An expression to evaluate while banc is the default
#'   autosegmentation
#' @rdname choose_banc
#' @export
#' @examples
#' \donttest{
#' with_banc(fafbseg::flywire_islatest('648518346498254576'))
#' }
#' \dontrun{
#' with_banc(fafbseg::flywire_latestid('648518346498254576'))
#' with_banc(fafbseg::flywire_latestid('648518346494405175'))
#' }
with_banc <- function(expr) {
  op <- choose_banc(set = TRUE)
  on.exit(options(op))
  force(expr)
}

# hidden
banc_fetch <- function(url, token=banc_token(), ...) {
  flywire_fetch(url, token=token, ...)
}

# hidden
`banc_ngl_segments<-` <- function (x, layer = NULL, value) {
  was_char <- is.character(x)
  baseurl <- if (was_char)
    x
  else NULL
  x = fafbseg::ngl_decode_scene(x)
  layers = fafbseg::ngl_layers(x)
  nls = fafbseg:::ngl_layer_summary(layers)
  sel = which(nls$type == "segmentation_with_graph")
  if (length(sel) == 0)
    sel = which(nls$visible & grepl("^segmentation", nls$type))
  if (length(sel) == 0)
    stop("Could not find a visible segmentation layer!")
  if (length(sel) > 1) {
    if(is.null(layer)){
      sel = 1
    }else{
      sel = match(layer,nls$name)
    }
  }
  if (is.null(value))
    value <- character()
  newsegs = fafbseg::ngl_segments(value, as_character = TRUE, must_work = FALSE)
  if (!all(fafbseg:::valid_id(newsegs)))
    warning("There are ", sum(!fafbseg:::valid_id(newsegs)), " invalid segments")
  x[["layers"]][[sel]][["segments"]] = newsegs
  nls$nhidden[sel] <- 1
  if (nls$nhidden[sel] > 0)
    x[["layers"]][[sel]][["hiddenSegments"]] = NULL
  if (was_char)
    as.character(x, baseurl = baseurl)
  else x
}

# Make a neuroglancer layer with 3D points in it, for synapse review
banc_annotation_layer <- function(data,
                                  layer = "sample",
                                  open = FALSE,
                                  rawcoords = NA,
                                  colpal = NULL){
  #data <- read_csv('/Users/abates/projects/flyconnectome/bancpipeline/tracing/2024-08-12_banc_synapse_sample_v1.csv')
  #data$pt_position<-data$`Coordinate 1`
  data$layer <- "synapse_sample"
  al <- ngl_annotation_layers(data[,c("pt_position", "layer")], rawcoords=rawcoords, colpal=colpal)
  sc<-fafbseg::ngl_decode_scene(banc_scene())
  sc2<-sc+al
  u<-as.character(sc2)
  su<-banc_shorturl(u)
  # # Extract annotations!
  # ngl_annotations(su)
  if(open){
    browseURL(u)
  }else{
    su
  }
}
