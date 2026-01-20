# Choose or (temporarily) use the banc autosegmentation

Choose or (temporarily) use the banc autosegmentation

## Usage

``` r
choose_banc(set = TRUE)

with_banc(expr)
```

## Arguments

- set:

  Whether or not to permanently set the banc autosegmentation as the
  default for
  [`fafbseg`](https://rdrr.io/pkg/fafbseg/man/fafbseg-package.html)
  functions.

- expr:

  An expression to evaluate while banc is the default autosegmentation

## Value

If `set=TRUE` a list containing the previous values of the relevant
global options (in the style of
[`options`](https://rdrr.io/r/base/options.html). If `set=FALSE` a named
list containing the option values.

## Details

`bancr` inherits a significant amount of infrastructure from the
[`fafbseg`](https://rdrr.io/pkg/fafbseg/man/fafbseg-package.html)
package. This has the concept of the *active* autosegmentation, which in
turn defines one or more R options containing URIs pointing to
voxel-wise segmentation, mesh etc data. These are normally contained
within a single neuroglancer URL which points to multiple data layers.
For banc this is the neuroglancer scene returned by
[`banc_scene`](https://flyconnectome.github.io/bancr/reference/banc_scene.md).

## Examples

``` r
if (FALSE) { # \dontrun{
choose_banc()
options()[grep("^fafbseg.*url", names(options()))]
} # }
# \donttest{
with_banc(fafbseg::flywire_islatest('648518346498254576'))
#> Warning: running command ''/home/runner/.cache/R/reticulate/uv/cache/archive-v0/r9AbzfZXekaq7Yq_R12Y8/bin/python' -m pip freeze' had status 1
#> Warning: using default setting
#> [1] FALSE
# }
if (FALSE) { # \dontrun{
with_banc(fafbseg::flywire_latestid('648518346498254576'))
with_banc(fafbseg::flywire_latestid('648518346494405175'))
} # }
```
