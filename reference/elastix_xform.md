# Perform Elastix Transform on 3D Points

This function applies an Elastix spatial transform to a set of 3D
points.

## Usage

``` r
elastix_xform(points, transform_file, copy_files = c(), return_logs = FALSE)
```

## Arguments

- points:

  A matrix with 3 columns or a data frame with x, y, z columns
  representing 3D points.

- transform_file:

  Path to the Elastix transform file, usually a `.txt` file, usually a
  `.txt` file.

- copy_files:

  Vector of additional file paths to copy to the temporary directory.

- return_logs:

  Logical, if TRUE, returns the Elastix log instead of transformed
  points.

## Value

A matrix of transformed 3D points, or Elastix logs if return_logs is
TRUE.

## Details

This function requires Elastix to be installed and added to the system
PATH. It creates a temporary directory for processing, applies the
Elastix transform, and cleans up afterwards.

## Examples

``` r
if (FALSE) { # \dontrun{
points <- matrix(rnorm(30), ncol = 3)
transformed_points <- elastix_xform(points, "path/to/transform.txt")
} # }
```
