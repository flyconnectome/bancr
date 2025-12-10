# Re-root BANC neuron skeleton at soma

This function re-roots a neuron skeleton represented as a `neuron`
object at the location of the corresponding soma in the `roots` data
frame. It uses the `root_id` in the skeleton object to identify the soma
location.

## Usage

``` r
banc_reroot(x, id = NULL, roots = NULL, estimate = TRUE, ...)

# S3 method for class 'neuron'
banc_reroot(x, id = NULL, roots = NULL, estimate = TRUE, ...)

# S3 method for class 'neuronlist'
banc_reroot(x, id = NULL, roots = NULL, estimate = TRUE, ...)
```

## Arguments

- x:

  A `banc.neurite` object representing the neuron skeleton.

- id:

  (Optional) The `root_id` of the neuron in the `roots` data frame. If
  NULL, it will be taken from the `x$root_id` slot.

- roots:

  A data frame containing information about root points, i.e. nuclei
  obtained using `bancr:::banc_roots()`. This data frame is assumed to
  have columns named `root_id` and `root_position_nm`, where
  `root_position_nm` specifies the 3D coordinates of the soma for each
  `root_id`.

- estimate:

  if `TRUE` and nucleus position is not in `roots`, then root is
  estimated as a leaf node furthest outside of the brain neuropil.

- ...:

  Methods passed to
  [`nat::nlapply`](https://rdrr.io/pkg/nat/man/nlapply.html).

## Value

The function returns the re-rooted `neuron` object.

## Examples

``` r
if (FALSE) { # \dontrun{
x <- banc_read_l2skel(..., simplify = FALSE)
roots <- bancr:::banc_roots()
re-rooted_neuron <- banc_reroot(x, roots = roots)
} # }
```
