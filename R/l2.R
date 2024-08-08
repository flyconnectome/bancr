#' Read L2 skeleton or dotprops for BANC neurons using fafbseg-py
#'
#' @description \code{banc_read_l2skel} reads one or more neurons as simplified
#'   L2 skeletons.
#'
#' @description \code{banc_read_l2dp} reads one or more neurons as simplified
#'   dotprops format. See \code{\link[fafbseg]{read_l2skel}}.
#'
#' @param dataset An optional CAVE dataset name (expert use only, by default
#'   will choose the standard banc dataset). See details.
#' @inheritParams fafbseg::read_l2skel
#'
#' @details \code{banc_read_l2dp} uses a special data structure for rapid
#'   download of the dotprops version of neurons required for NBLASTing. It
#'   leverages the python navis / fafbseg-py packages and you will need to
#'   install these, typically using the \code{\link[fafbseg]{simple_python}}
#'   function.
#'
#'   \code{banc_read_l2skel} treats the dataset argument a little differently
#'   than \code{banc_read_l2dp} because it actually needs to identify two data sources
#'   a CAVE data
#'
#'   See \code{\link[fafbseg]{read_l2skel}} for additional details of
#'
#' @return a \code{\link{neuronlist}} containing one or more
#'   \code{\link{neuron}} or \code{\link{dotprops}} objects. Note that neurons
#'   will be calibrated in \code{nm} while dotprops will be calibrated in microns.
#' @export
#'
#' @examples
#'
#' \dontrun{
#' # one time install of necessary python packages
#' fafbseg::simple_python(pkgs="fafbseg")
#'
#' dna02=c("720575941478275714", "720575941512946243")
#' dna02.latest=banc_latestid(dna02)
#' dna02.dps <- banc_read_l2dp(dna02.latest)
#'
#' # plot those
#' nclear3d()
#' plot3d(dna02.dps, lwd=3)
#' # nb dotprops are always in microns
#' wire3d(banc.surf/1e3, col='grey')
#'
#' nclear3d()
#' dna02.skel <- banc_read_l2skel(dna02.latest)
#' plot3d(dna02.skel, lwd=2)
#' # nb neuron skeletons are in nm
#' wire3d(banc.surf, col='grey')
#' }
banc_read_l2dp <- function(id, OmitFailures=TRUE, dataset=NULL, ...) {
  id=banc_ids(id)
  if(is.null(dataset))
    dataset=with_banc(getOption("fafbseg.cave.datastack_name"))
  fafbseg::read_l2dp(id, dataset=dataset, OmitFailures=OmitFailures, ...)
}

#' @export
#' @rdname banc_read_l2dp
banc_read_l2skel <- function(id, OmitFailures=TRUE, dataset=NULL, ...) {
  id=banc_ids(id)
  # partial duplication of fafbseg::read_l2skel for older versions of fafbseg-py
  fp=fafbseg:::check_fafbsegpy()

  if("set_default_dataset" %in% names(fp$flywire)) {
    # new fafbseg-py, everything is simpler
    return(with_banc(fafbseg::read_l2skel(id)))
  }
  # this used to work with banc and older fafbseg-py, but not sure if it still works ...
  # manually set the cloudvolume url / cave datastack name
  # see https://flyconnectome.slack.com/archives/C342Q3H4Y/p1694682637820519
  if(is.null(dataset)) {
    ops=choose_banc(set=F)
    fp$flywire$utils$FLYWIRE_URLS$banc = ops$fafbseg.cloudvolume.url
    fp$flywire$utils$CAVE_DATASETS$banc = ops$fafbseg.cave.datastack_name
    dataset="banc"
  }
  sk=fp$flywire$l2_skeleton(id, omit_failures = OmitFailures, dataset=dataset, ...)
  fafbseg:::navis2nat_neuronlist(sk)
}

#' @title Re-root BANC neuron skeleton at soma
#'
#' @description
#' This function re-roots a neuron skeleton represented as a `neuron`
#' object at the location of the corresponding soma in the `banc_nuclei` data
#' frame. It uses the `root_id` in the skeleton object to identify the soma
#' location.
#'
#' @param x A `banc.neurite` object representing the neuron skeleton.
#' @param id (Optional) The `root_id` of the neuron in the `banc_nuclei` data
#' frame. If NULL, it will be taken from the `x$root_id` slot.
#' @param banc_nuclei A data frame containing information about nuclei
#' obtained using `bancr::banc_nuclei()`. This data frame is assumed to have
#' columns named `root_id` and `nucleus_position_nm`, where `nucleus_position_nm`
#' specifies the 3D coordinates of the soma for each `root_id`.
#' @param estimate if \code{TRUE} and nucleus position is not in `banc_nuclei`,
#' then root is estimated as a leaf node furthest outside of the brain neuropil.
#' @param ... Methods passed to \code{nat::nlapply}.
#'
#' @return The function returns the re-rooted `neuron` object.
#'
#' @examples
#' \dontrun{
#' x <- banc_read_l2skel(..., simplify = FALSE)
#' banc_nuclei <- banc_nuclei()
#' re-rooted_neuron <- banc_reroot(x, banc_nuclei = banc_nuclei)
#' }
#' @export
#' @rdname banc_reroot
banc_reroot <- function(x, id = NULL, banc_nuclei = bancr::banc_nuclei(rawcoords = FALSE), estimate = TRUE, ...) UseMethod("banc_reroot")

#' @rdname banc_reroot
#' @method banc_reroot neuron
#' @export
banc_reroot.neuron <- function(x, id = NULL, banc_nuclei = bancr::banc_nuclei(rawcoords = FALSE), estimate = TRUE, ...){
  if(is.null(id)){
    id <- x$id
  }
  if(is.null(id)){
    stop("a root_id in banc_nuclei must be given")
  }
  df <- subset(banc_nuclei,banc_nuclei$root_id==id & nucleus_id!="0" & !is.na(banc_nuclei$nucleus_position_nm))
  if(nrow(df)){
    soma <- nat::xyzmatrix(df$nucleus_position_nm)[1,]
    x <- nat::reroot(x = x, point = c(soma))
    x$tags$soma <- nat::rootpoints(x)
  }else if (estimate){ # As best we can
    warning(sprintf("no valid nucleus ID detecting for %s, estimating root point"),id)
    leaves <- nat::endpoints(x)
    npoints1 <- nat::xyzmatrix(x)[leaves,]
    if(nrow(npoints1)){npoints=npoints1}
    pin <- nat::pointsinside(x = npoints, surf = bancr::banc_neuropil.surf)
    npoints2 <- data.frame(npoints[!pin,])
    if(nrow(npoints2)){npoints=npoints2}
    pin <- nat::pointsinside(x = npoints, surf = bancr::banc_neck_connective.surf)
    npoints3 <- data.frame(npoints[!pin,])
    if(nrow(npoints3)){npoints=npoints3}
    npoints$nucleus_id <- 0
    npoints$root_id <- id
    nearest <- nabor::knn(query = nat::xyzmatrix(npoints),
                          data = rbind(nat::xyzmatrix(bancr::banc_neuropil.surf),
                                       nat::xyzmatrix(bancr::banc_neck_connective.surf)), k = 1)
    soma <-nat::xyzmatrix(npoints)[which.max(nearest$nn.dists),]
    x <- nat::reroot(x = x, point = c(soma))
  }else{
    warning(sprintf("no valid nucleus ID detecting for %s, no action taken"),id)
  }
  x
}

#' @rdname banc_reroot
#' @method banc_reroot neuronlist
#' @export
banc_reroot.neuronlist <- function(x, id = NULL, banc_nuclei = bancr::banc_nuclei(), estimate = TRUE, ...){
  if(is.null(id)){
    id <- names(x)
  }
  nat::nlapply(x, FUN = banc_reroot.neuron, banc_nuclei = banc_nuclei, id = id, ...)
}



