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
#' dnp42=c("648518346507131167", "648518346485772414")
#' dnp42.latest=banc_latestid(dnp42)
#' dnp42.dps <- banc_read_l2dp(dnp42.latest)
#'
#' # plot those
#' nclear3d()
#' plot3d(dnp42.dps, lwd=3)
#' # nb dotprops are always in microns
#' wire3d(banc.surf/1e3, col='grey')
#'
#' nclear3d()
#' dnp42.skel <- banc_read_l2skel(dnp42.latest)
#' plot3d(dnp42.skel, lwd=2)
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
