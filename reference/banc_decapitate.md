# Subset points to be in the brain or in the VNC

Subset points to be in the brain or in the VNC

## Usage

``` r
banc_decapitate(x, y.cut = 325000, invert = FALSE, ...)

# S3 method for class '`NULL`'
banc_decapitate(x, y.cut = 325000, invert = FALSE, ...)

# S3 method for class 'neuron'
banc_decapitate(x, y.cut = 325000, invert = FALSE, ...)

# S3 method for class 'neuronlist'
banc_decapitate(x, y.cut = 325000, invert = FALSE, ...)

# S3 method for class 'matrix'
banc_decapitate(x, y.cut = 325000, invert = FALSE, ...)

# S3 method for class 'data.frame'
banc_decapitate(x, y.cut = 325000, invert = FALSE, ...)

# S3 method for class 'mesh3d'
banc_decapitate(x, y.cut = 325000, invert = FALSE, ...)

# S3 method for class 'hxsurf'
banc_decapitate(x, y.cut = 325000, invert = FALSE, ...)
```

## Arguments

- x:

  an object with 3d points to be subsetted, e.g. an xyz matrix, a
  `neuron`, `neuronlist` or a `mesh3d` object. Points must be in native
  BANC space, i.e. plottable inside `banc.surf`.

- y.cut:

  Numeric, the Y-axis cut point, in nanometers, in BANC space, that
  separates the head from the neck and ventral nerve cord. For fitting
  to the MANC data set, a cut height of `y.cut=5e05` seems good.

- invert:

  if `TRUE` returns brain points, if `FALSE` returns VNC points.

- ...:

  Additional arguments passed to
  [`nlapply`](https://rdrr.io/pkg/nat/man/nlapply.html) and then
  [`prune_vertices`](https://rdrr.io/pkg/nat/man/prune_vertices.html)

## Value

Remove points above or below the midsection of the neck connective of
BANC.

## See also

[`banc.surf`](https://flyconnectome.github.io/bancr/reference/banc.surf.md)

## Examples

``` r
# \donttest{
# DNa02
m = banc_read_neuron_meshes("720575941478275714")
#> Warning: running command ''/home/runner/.cache/R/reticulate/uv/cache/archive-v0/Zu4cyt5m5nm-eRob5h0Pk/bin/python' -m pip freeze' had status 1
#> Warning: using default setting
#>   downloading meshes
#> No module named 'cloudvolume'
#> Error: Please install the python cloudvolume package:
#> This should normally work:
#> fafbseg::simple_python('basic')
#> For more details see ?simple_python or the cloud-volume docshttps://github.com/seung-lab/cloud-volume#setup
#> If you have already installed cloudvolume but it is not found
#> then R probably can't find the relevant version of Python
#> Do:
#> usethis::edit_r_environ()
#>  to point to the right python
#> e.g. RETICULATE_PYTHON="/opt/miniconda3/envs/r-reticulate/bin/python"
m.brain = banc_decapitate(m)
#> Error: object 'm' not found
m.vnc = banc_decapitate(m, invert = TRUE)
#> Error: object 'm' not found
# }
if (FALSE) { # \dontrun{
plot3d(m.brain, col = "red")
plot3d(m.vnc, col = "cyan")
plot3d(banc.surf, col = "grey", alpha = 0.1)
} # }
```
