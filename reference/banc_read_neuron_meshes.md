# Read one or more BANC neuron and nuclei meshes

Read one or more BANC neuron and nuclei meshes

## Usage

``` r
banc_read_neuron_meshes(ids, savedir = NULL, format = c("ply", "obj"), ...)

banc_read_mitochondria_mesh(
  ids,
  lod = 0L,
  savedir = NULL,
  method = c("vf", "ply"),
  ...
)

banc_read_nuclei_mesh(
  ids,
  lod = 0L,
  savedir = NULL,
  method = c("vf", "ply"),
  ...
)
```

## Arguments

- ids:

  One or more root ids

- savedir:

  An optional location to save downloaded meshes. This acts as a simple
  but effective cache since flywire neurons change id whenever they are
  edited.

- format:

  whether to save meshes in Wavefront obj or Stanford poly format. obj
  is the default but ply is a simpler and more compact format.

- ...:

  Additional arguments passed to
  `fafbseg::`[`read_cloudvolume_meshes`](https://rdrr.io/pkg/fafbseg/man/read_cloudvolume_meshes.html)

- lod:

  The level of detail (highest resolution is 0, default of 2 gives a
  good overall morphology while 3 is also useful and smaller still).

- method:

  How to treat the mesh object returned from neuroglancer, i.e. as a
  `mesh3d` object or a `ply` mesh.

## Value

A [`neuronlist`](https://rdrr.io/pkg/nat/man/neuronlist.html) containing
one or more `mesh3d` objects. See
`nat::`[`read.neurons`](https://rdrr.io/pkg/nat/man/read.neurons.html)
for details.

## See also

`fafbseg::`[`read_cloudvolume_meshes`](https://rdrr.io/pkg/fafbseg/man/read_cloudvolume_meshes.html)

## Examples

``` r
if (FALSE) { # \dontrun{
neuron.mesh <- banc_read_neuron_meshes("720575941478275714")
plot3d(neuron.mesh, alpha = 0.1)
nucleus.mesh <- banc_read_nuclei_mesh("72903876004544795")
plot3d(nucleus.mesh, col = "black")
} # }
```
