% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/l2.R
\name{banc_read_l2dp}
\alias{banc_read_l2dp}
\alias{banc_read_l2skel}
\title{Read L2 skeleton or dotprops for BANC neurons using fafbseg-py}
\usage{
banc_read_l2dp(id, OmitFailures = TRUE, dataset = NULL, ...)

banc_read_l2skel(id, OmitFailures = TRUE, dataset = NULL, ...)
}
\arguments{
\item{id}{One or more flywire ids}

\item{OmitFailures}{Whether or not to drop neurons that cannot be read from
the results (rather than erroring out). Default \code{TRUE}.}

\item{dataset}{An optional CAVE dataset name (expert use only, by default
will choose the standard banc dataset). See details.}

\item{...}{Additional arguments passed to the
\code{fafbseg.flywire.l2_skeleton} or
\code{fafbseg.flywire.l2_dotprops}functions.}
}
\value{
a \code{\link{neuronlist}} containing one or more
\code{\link{neuron}} or \code{\link{dotprops}} objects. Note that neurons
will be calibrated in \code{nm} while dotprops will be calibrated in microns.
}
\description{
\code{banc_read_l2skel} reads one or more neurons as simplified
L2 skeletons.

\code{banc_read_l2dp} reads one or more neurons as simplified
dotprops format. See \code{\link[fafbseg]{read_l2skel}}.
}
\details{
\code{banc_read_l2dp} uses a special data structure for rapid
download of the dotprops version of neurons required for NBLASTing. It
leverages the python navis / fafbseg-py packages and you will need to
install these, typically using the \code{\link[fafbseg]{simple_python}}
function.

\code{banc_read_l2skel} treats the dataset argument a little differently
than \code{banc_read_l2dp} because it actually needs to identify two data sources
a CAVE data

See \code{\link[fafbseg]{read_l2skel}} for additional details of
}
\examples{

\dontrun{
# one time install of necessary python packages
fafbseg::simple_python(pkgs="fafbseg")

dna02=c("720575941478275714", "720575941512946243")
dna02.latest=banc_latestid(dna02)
dna02.dps <- banc_read_l2dp(dna02.latest)

# plot those
nclear3d()
plot3d(dna02.dps, lwd=3)
# nb dotprops are always in microns
wire3d(banc.surf/1e3, col='grey')

nclear3d()
dna02.skel <- banc_read_l2skel(dna02.latest)
plot3d(dna02.skel, lwd=2)
# nb neuron skeletons are in nm
wire3d(banc.surf, col='grey')
}
}
