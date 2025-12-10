# Apply Elastix Transform using Navis

Applies an Elastix transform to 3D points using the Navis Python
library.

## Usage

``` r
navis_elastix_xform(x, transform_file)
```

## Arguments

- x:

  Matrix or data frame of 3D points.

- transform_file:

  Path to the Elastix transform file, usually a `.txt` file.

## Value

A matrix of transformed 3D points.

## Details

This function requires the reticulate R package and the Navis Python
library.

## Examples

``` r
if (FALSE) { # \dontrun{
neuron.mesh <- banc_read_neuron_meshes("720575941478275714")
points <- nat::xyzmatrix(neuron.mesh)
transformed_points <- navis_elastix_xform(points,
transform_file = "brain_240721/BANC_to_template.txt")
points3d(points)
plot3d(nat.flybrains::JRC2018F)
} # }
```
