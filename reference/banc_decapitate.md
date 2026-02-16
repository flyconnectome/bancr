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
if (FALSE) { # \dontrun{
# DNa02
m = banc_read_neuron_meshes("720575941478275714")
m.brain = banc_decapitate(m)
m.vnc = banc_decapitate(m, invert = TRUE)
} # }
if (FALSE) { # \dontrun{
plot3d(m.brain, col = "red")
plot3d(m.vnc, col = "cyan")
plot3d(banc.surf, col = "grey", alpha = 0.1)
} # }
```
