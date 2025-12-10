# Read BANC euroglancer meshes, e.g., ROI meshes

Read BANC euroglancer meshes, e.g., ROI meshes

## Usage

``` r
banc_read_neuroglancer_mesh(
  x = 1,
  url = paste0("https://www.googleapis.com/storage/v1/b/",
    "zetta_lee_fly_cns_001_kisuk/o/final%2Fv2%2F",
    "volume_meshes%2Fmeshes%2F{x}%3A0.drc?alt=media",
    "&neuroglancer=610000b05b6497edcf20b78f29516970"),
  ...
)
```

## Arguments

- x:

  the numeric identifier that specifies the mesh to read, defaults to
  `1` the BANC outline mesh.

- url:

  the URL that directs `bancr` to where BANC meshes are stored.

- ...:

  additional arguments to `GET`

## Value

a mesh3d object for the specified mesh.

## See also

[`banc_read_neuron_meshes`](https://flyconnectome.github.io/bancr/reference/banc_read_neuron_meshes.md)

## Examples

``` r
if (FALSE) { # \dontrun{
banc.mesh  <- banc_read_neuroglancer_mesh()
} # }
```
