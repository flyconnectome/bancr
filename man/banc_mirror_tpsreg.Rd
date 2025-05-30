% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{banc_mirror_tpsreg}
\alias{banc_mirror_tpsreg}
\title{Thin-Plate Spline Registration for Mirroring in BANC Space}
\format{
An object of class \code{tpsCoeff} created using \code{Morpho::computeTransform}.
It contains the following components:
\describe{
\item{Lw}{Matrix of TPS coefficients}
\item{refmat}{Reference matrix (source landmarks)}
\item{tarmat}{Target matrix (mirrored landmarks)}
\item{lattice}{3D array representing the deformation grid}
\item{lambda}{Smoothing parameter used in the TPS computation}
\item{scale}{Logical indicating whether scaling was used}
\item{reflection}{Logical indicating whether reflection was allowed}
}
}
\source{
Derived from landmark correspondences between original and mirrored points
in the BANC space, possibly utilizing transformations to and from the
JRC2018F space for accurate mirroring.
}
\usage{
data(banc_mirror_tpsreg)
}
\description{
A thin-plate spline (TPS) registration object that mirrors 3D points
within the BANC (Buhmann et al. Adult Neural Connectome) coordinate system.
}
\details{
This TPS registration was computed to allow mirroring of points directly within
the BANC coordinate system. It provides a smooth, interpolated transformation
for any point in the BANC space to its mirrored counterpart, accounting for
any asymmetries in the BANC reference brain.
}
\examples{
\dontrun{
data(banc_mirror_tpsreg)

# Mirror BANC points using the TPS registration
banc_points <- matrix(rnorm(300), ncol=3)  # Example BANC points
mirrored_points <- Morpho::tps3d(banc_points, banc_mirror_tpsreg)
}

}
\seealso{
\code{\link{banc_mirror}} for the function that performs mirroring
of points in BANC space.
\code{\link{banc_to_jrc2018f_tpsreg}} for the TPS registration between
BANC and JRC2018F spaces.
\code{\link[Morpho]{computeTransform}} for details on creating tpsCoeff objects.
}
\keyword{datasets}
