# Handle raw and nm calibrated banc coordinates

`banc_voxdims` returns the image voxel dimensions which are normally
used to scale between **raw** and **nm** coordinates.

## Usage

``` r
banc_voxdims(url = choose_banc(set = FALSE)[["fafbseg.sampleurl"]])

banc_nm2raw(x, vd = banc_voxdims())

banc_raw2nm(x, vd = banc_voxdims())
```

## Arguments

- url:

  Optional neuroglancer URL containing voxel size. Defaults to
  `getOption("fafbseg.sampleurl")` as set by
  [`choose_banc`](https://flyconnectome.github.io/bancr/reference/choose_banc.md).

- x:

  3D coordinates in any form compatible with
  [`xyzmatrix`](https://rdrr.io/pkg/nat/man/xyzmatrix.html)

- vd:

  The voxel dimensions in nm. Expert use only. Normally found
  automatically.

## Value

For `banc_voxdims` A 3-vector

for `banc_raw2nm` and `banc_nm2raw` an Nx3 matrix of coordinates

## Details

relies on nat \>= 1.10.4

## Examples

``` r
banc_voxdims()
#> Warning: running command ''/home/runner/.cache/R/reticulate/uv/cache/archive-v0/b0-oak1mtUw0zpTvaBS0m/bin/python' -m pip freeze' had status 1
#> Warning: using default setting
#> Warning: Multiple segmentation layers. Choosing first!
#>  x  y  z 
#>  4  4 45 
#> attr(,"units")
#> [1] "nm"
banc_raw2nm(c(159144, 22192, 3560))
#> Warning: running command ''/home/runner/.cache/R/reticulate/uv/cache/archive-v0/b0-oak1mtUw0zpTvaBS0m/bin/python' -m pip freeze' had status 1
#> Warning: using default setting
#>           X     Y      Z
#> [1,] 636576 88768 160200
banc_raw2nm('159144 22192 3560')
#> Warning: running command ''/home/runner/.cache/R/reticulate/uv/cache/archive-v0/b0-oak1mtUw0zpTvaBS0m/bin/python' -m pip freeze' had status 1
#> Warning: using default setting
#>           X     Y      Z
#> [1,] 636576 88768 160200
if (FALSE) { # \dontrun{
banc_nm2raw(clipr::read_clip())
} # }
```
