# Thin-Plate Spline Registration from BANC to JRC2018F template brain

A thin-plate spline (TPS) registration object that transforms 3D points
from the BANC nanometer coordinate system to the D. melanogaster
template brain JRC2018F coordinate system.

## Usage

``` r
data(banc_to_jrc2018f_tpsreg)

jrc2018f_to_banc_tpsreg

jrcvnc2018f_to_banc_tpsreg

banc_to_jrcvnc2018f_tpsreg
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

  Target matrix (target landmarks)

- lattice:

  3D array representing the deformation grid

- lambda:

  Smoothing parameter used in the TPS computation

- scale:

  Logical indicating whether scaling was used

- reflection:

  Logical indicating whether reflection was allowed

An object of class `tpsCoeff` of length 4.

An object of class `tpsCoeff` of length 4.

An object of class `tpsCoeff` of length 4.

## Source

Derived from Elastix registration results using the `banc_to_JRC2018F`
function and landmark correspondences extracted from that registration.

## Details

This TPS registration was computed based on landmark correspondences
derived from an Elastix registration between the BANC and JRC2018F
spaces. It provides a smooth, interpolated transformation for any point
in the BANC space to its corresponding location in the JRC2018F space.

## See also

[`banc_to_JRC2018F`](https://flyconnectome.github.io/bancr/reference/banc_to_JRC2018F.md)
for the function that performs transformations between BANC and JRC2018F
spaces.
[`computeTransform`](https://rdrr.io/pkg/Morpho/man/computeTransform.html)
for details on creating tpsCoeff objects.

## Examples

``` r
if (FALSE) { # \dontrun{
data(banc_to_jrc2018f_tpsreg)

# Transform BANC points to JRC2018F using the TPS registration
banc_points <- matrix(rnorm(300), ncol=3)  # Example BANC points
jrc2018f_points <- Morpho::tps3d(banc_points, banc_to_jrc2018f_tpsreg)
} # }
```
