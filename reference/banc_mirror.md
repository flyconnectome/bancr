# Mirror BANC Connectome Points

This function mirrors 3D points in the BANC (Buhmann et al. Adult Neural
Connectome) coordinate system by transforming to JRC2018F, mirroring,
and transforming back.

## Usage

``` r
banc_mirror(
  x,
  banc.units = c("nm", "um", "raw"),
  subset = NULL,
  inverse = FALSE,
  transform_files = NULL,
  method = c("tpsreg", "elastix", "navis_elastix_xform"),
  ...
)
```

## Arguments

- x:

  An object containing 3D points (must be compatible with
  nat::xyzmatrix).

- banc.units:

  Character string specifying the banc.units of the input points. Must
  be one of "nm" (nanometers), "um", or "raw" (BANC raw banc.units).
  Default is "nm".

- subset:

  Optional. A logical vector or expression to subset the input object.

- inverse:

  Logical. Not used in this function, kept for compatibility with
  banc_to_JRC2018F.

- transform_files:

  Optional. A vector of two file paths for custom transform files. If
  NULL, uses default files.

- method:

  Character string specifying the transformation method. Must be either
  "elastix" or "tpsreg". Default is "elastix".

- ...:

  Additional arguments passed to
  [`nat.templatebrains::mirror_brain`](https://natverse.org/nat.templatebrains/reference/mirror_brain.html).

## Value

The input object with mirrored 3D points.

## Details

This function performs mirroring of BANC points by first transforming
them to the JRC2018F coordinate system, applying the mirroring
operation, and then transforming them back to BANC. It can use either
Elastix transforms or thin-plate spline registration for the coordinate
system transformations.

## See also

[`banc_to_JRC2018F`](https://flyconnectome.github.io/bancr/reference/banc_to_JRC2018F.md)
for the underlying transformation function.
[`nat.templatebrains::mirror_brain`](https://natverse.org/nat.templatebrains/reference/mirror_brain.html)
for the mirroring operation in JRC2018F space.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example using saved tpsreg
banc_neuropil.surf.m <- banc_mirror(banc_neuropil.surf, method = "tpsreg")
clear3d()
banc_view()
plot3d(banc_neuropil.surf, alpha = 0.5, col = "lightgrey")
plot3d(banc_neuropil.surf.m, alpha = 0.5, col = "green")

# Example using custom Elastix transforms
choose_banc()
rootid <- "720575941626035769"
neuron.mesh <- banc_read_neuron_meshes(rootid)

# Show neuron in BANC neuropil
banc_view()
plot3d(neuron.mesh, col = "red")
plot3d(banc_neuropil.surf, alpha = 0.1, col = "lightgrey")

# Show only the portion in the brain
neuron.mesh.brain <- banc_decapitate(neuron.mesh, invert = TRUE)

# Mirror in BANC space
neuron.mesh.mirror <- banc_mirror(neuron.mesh.brain,
transform_files = c("brain_240721/BANC_to_template.txt",
 "brain_240721/template_to_BANC.txt"))
plot3d(neuron.mesh.mirror, col = "cyan")
} # }
```
