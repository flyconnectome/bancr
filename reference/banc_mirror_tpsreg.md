# Thin-Plate Spline Registration for Mirroring in BANC Space

A thin-plate spline (TPS) registration object that mirrors 3D points
within the BANC (Buhmann et al. Adult Neural Connectome) coordinate
system.

## Usage

``` r
data(banc_mirror_tpsreg)
```

## Format

An object of class `tpsCoeff` created using
[`Morpho::computeTransform`](https://rdrr.io/pkg/Morpho/man/computeTransform.html).
It contains the following components:

- Lw:

  Matrix of TPS coefficients

- refmat:

  Reference matrix (source landmarks)

- tarmat:

  Target matrix (mirrored landmarks)

- lattice:

  3D array representing the deformation grid

- lambda:

  Smoothing parameter used in the TPS computation

- scale:

  Logical indicating whether scaling was used

- reflection:

  Logical indicating whether reflection was allowed

## Source

Derived from landmark correspondences between original and mirrored
points in the BANC space, possibly utilizing transformations to and from
the JRC2018F space for accurate mirroring.

## Details

This TPS registration was computed to allow mirroring of points directly
within the BANC coordinate system. It provides a smooth, interpolated
transformation for any point in the BANC space to its mirrored
counterpart, accounting for any asymmetries in the BANC reference brain.

## See also

[`banc_mirror`](https://flyconnectome.github.io/bancr/reference/banc_mirror.md)
for the function that performs mirroring of points in BANC space.
[`banc_to_jrc2018f_tpsreg`](https://flyconnectome.github.io/bancr/reference/banc_to_jrc2018f_tpsreg.md)
for the TPS registration between BANC and JRC2018F spaces.
[`computeTransform`](https://rdrr.io/pkg/Morpho/man/computeTransform.html)
for details on creating tpsCoeff objects.

## Examples

``` r
if (FALSE) { # \dontrun{
data(banc_mirror_tpsreg)

# Mirror BANC points using the TPS registration
banc_points <- matrix(rnorm(300), ncol=3)  # Example BANC points
mirrored_points <- Morpho::tps3d(banc_points, banc_mirror_tpsreg)
} # }
```
