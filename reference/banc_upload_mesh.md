# Upload Mesh to Google Cloud Storage for BANC Neuroglancer

This function uploads a mesh.obj file to a Google Cloud Storage bucket
that can be accessed by BANC neuroglancer. It uses reticulate to call
Python code and manages a conda environment.

The underlying Python code is based on the bikinibottom package:
<https://github.com/jasper-tms/bikini-bottom/tree/main>

## Usage

``` r
banc_upload_mesh(
  mesh,
  mesh_id,
  vol = "precomputed://gs://lee-lab_brain-and-nerve-cord-fly-connectome/volume_meshes/",
  compress = TRUE,
  overwrite = FALSE
)
```

## Arguments

- mesh:

  Character or mesh3d object. Either a path to the .obj file or an R
  mesh object.

- mesh_id:

  Integer. ID for the mesh.

- vol:

  Character. URL for the precomputed volume.

- compress:

  Logical. Whether to compress the mesh. Default is `TRUE`.

- overwrite:

  Logical. Whether or not to overwrite an extant file of the same name
  and path at the location specified by `vol`+`mesh_id`. Default is
  `FALSE`.

## Value

Invisible `NULL`. The function is called for it upload effect.

## Details

For this function to work the global environment variable
`GOOGLE_APPLICATION_CREDENTIALS` to something like
`"$HOME/.config/gcloud/legacy_credentials/firstname.secondname@gmail.com/adc.json"`.
This file can be obtains by installing `gsutil` and running
`gsutil init` in the terminal.

## References

bikinibottom package:
<https://github.com/jasper-tms/bikini-bottom/tree/main>

## Examples

``` r
if (FALSE) { # \dontrun{
banc_upload_mesh(
mesh = "/Users/abates/projects/flyconnectome/bancpipeline/deformetrica/obj/banc_brain_neuropil.obj",
mesh_id = 3,
vol = "precomputed://gs://lee-lab_brain-and-nerve-cord-fly-connectome/volume_meshes/meshes"
)
} # }
```
