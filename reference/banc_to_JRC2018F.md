# Transform Points between BANC Connectome and JRC2018F Template Brain

This function transforms 3D points between the BANC (Buhmann et al.
Adult Neural Connectome) coordinate system and the D. melanogaster
template brain JRC2018F coordinate system. Transforming to JRC2018F
helps move data from BANC into a more standard reference space for
comparison with light level data. JRC2018F is used by Janelia Research
Campus for their light-level registered data, including with
Neuronbridge (<https://neuronbridge.janelia.org/>). This transformation
is a first step in enabling users to match BANC connectome
reconstructions with genetic resources for wetlab experimentation.

## Usage

``` r
banc_to_JRC2018F(
  x,
  region = c("brain", "vnc"),
  banc.units = c("nm", "um", "raw"),
  subset = NULL,
  inverse = FALSE,
  transform_file = NULL,
  method = c("tpsreg", "elastix", "navis_elastix_xform")
)
```

## Arguments

- x:

  An object containing 3D points (must be compatible with
  nat::xyzmatrix).

- region:

  Whether this transform is for the JRC2018F brainspace (default) ot the
  JRCVNC2018F VNC template (only alternative).

- banc.units:

  Character string specifying the units of the BANC space data (input or
  output, depending on the inverse argument). Must be one of "nm"
  (nanometers), "um", or "raw" (BANC raw banc.units). Default is "nm".

- subset:

  Optional. A logical vector or expression to subset the input object.

- inverse:

  Logical. If TRUE, performs the inverse transformation (JRC2018F to
  BANC). Default is FALSE.

- transform_file:

  Optional. Path to a custom transform file. If NULL, uses default
  files.

- method:

  Character string specifying the transformation method. Must be either
  "elastix" or "tpsreg". Default is "elastix".

## Value

The input object with transformed 3D points.

## Details

This function applies either an Elastix transform or a thin-plate spline
registration to convert points between the BANC and JRC2018F coordinate
systems. It handles unit conversions as necessary.

The default transformation files are included with the package and are
located in the 'inst/extdata/brain_240721' directory.

## See also

[`elastix_xform`](https://flyconnectome.github.io/bancr/reference/elastix_xform.md)
for the underlying Elastix transformation function.
[`banc_raw2nm`](https://flyconnectome.github.io/bancr/reference/banc_voxdims.md)
and
[`banc_nm2raw`](https://flyconnectome.github.io/bancr/reference/banc_voxdims.md)
for unit conversion functions.

## Examples

``` r
if (FALSE) { # \dontrun{
### BRAIN EXAMPLE ####
# Transform points from BANC to JRC2018F
transformed_points <- banc_to_JRC2018F(points, banc.units = "nm")

# Use a custom transform file
custom_transformed <- banc_to_JRC2018F(points, transform_file = "path/to/custom/transform.txt")

# Where the default transform files are located:
banc_to_JRC2018F_file <- system.file(file.path("extdata","brain_240721"),
"BANC_to_template.txt", package="bancr")
JRC2018F_to_banc_file <- system.file(file.path("extdata","brain_240721"),
"template_to_BANC.txt", package="bancr")

### VNC EXAMPLE ####
library(malevnc)
library(nat.jrcbrains)
nat.jrcbrains::register_saalfeldlab_registrations()

# Get DNa02 axons from the MANC project
DNa02s=read_manc_meshes('DNa02')
plot3d(JRCVNC2018U)

# Transform into JRCVNC2918F
## nb convert from nm to microns
DNa02s.jrcvnc2018f=xform_brain(DNa02s/1e3, reference = "JRCVNC2018F", sample="MANC")
plot3d(DNa02s.jrcvnc2018f, co = "red")
plot3d(JRCVNC2018F)

# Transform into the BANC
DNa02s.banc <- banc_to_JRC2018F(DNa02s.jrcvnc2018f, region="VNC", method="tpsreg")
open3d()
plot3d(DNa02s.banc, co = "blue")
plot3d(banc_vnc_neuropil.surf)

} # }
```
