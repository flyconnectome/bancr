% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ids.R
\name{banc_leaves}
\alias{banc_leaves}
\title{Find the supervoxel identifiers of a banc neuron}
\usage{
banc_leaves(x, integer64 = TRUE, ...)
}
\arguments{
\item{x}{One or more FlyWire segment ids}

\item{integer64}{Whether to return ids as integer64 type (the default, more
compact but a little fragile) rather than character (when \code{FALSE}).}

\item{...}{additional arguments passed to \code{\link{flywire_leaves}}}
}
\value{
A vector of supervoxel ids
}
\description{
Find the supervoxel identifiers of a banc neuron
}
\examples{
\dontrun{
svids=banc_leaves("720575941478275714")
head(svids)
}
}
\seealso{
\code{\link{flywire_leaves}}

Other banc-ids: 
\code{\link{banc_cellid_from_segid}()},
\code{\link{banc_ids}()},
\code{\link{banc_islatest}()},
\code{\link{banc_latestid}()},
\code{\link{banc_rootid}()},
\code{\link{banc_xyz2id}()}
}
\concept{banc-ids}
