# Simplified tissue and neuropil surfaces for BANC

`banc.surf` is unsymmetrical and not a normalised version of the mesh.
It is the outline of the dataset in nanometers.`banc_neuropil.surf`
represents the synaptic neuropil. Built from the BANC synapse cloud, but
not optimised to include 100% of bona fide presynapses.
`banc_neuropils.surf` contains the standard Ito et al., 2014 brain
neuropil volumes transformed into BANC space. `banc_al.surf` contains
the standard Bates and Schlegel et al, 2021 right antennal lobe
glomeruli brain neuropil volumes transformed into BANC space. These
neuropils may also be seen in neuroglancer,
[here](https://spelunker.cave-explorer.org/#!middleauth+https://global.daf-apis.com/nglstate/api/v1/6237404072509440).

A BANC neuroglaner scene can be directed to a google cloud storage
location, where BANC-transformed standard neuropil meshes reside. The
source is
`precomputed://gs://lee-lab_brain-and-nerve-cord-fly-connectome/volume_meshes`
They can be plotted in neuroglancer by adding
[this](https://spelunker.cave-explorer.org/#!middleauth+https://global.daf-apis.com/nglstate/api/v1/6237404072509440)
location, entering the `Seg.` pane and entering the number or name that
corresponds to the correct mesh. This data frame gives the mesh name to
number correspondences.

## Usage

``` r
banc.surf

banc_neuropil.surf

banc_brain_neuropil.surf

banc_vnc_neuropil.surf

banc_neck_connective.surf

banc_brain_neuropils.surf

banc_al.surf

banc_vnc_neuropils.surf

banc_vnc_tracts.surf

banc_vnc_nerves.surf

banc_hemibrain.surf

banc_volumes.df
```

## Format

An object of class `hxsurf` (inherits from `list`) of length 4.

An object of class `hxsurf` (inherits from `list`) of length 4.

An object of class `hxsurf` (inherits from `list`) of length 4.

An object of class `hxsurf` (inherits from `list`) of length 4.

An object of class `hxsurf` (inherits from `list`) of length 4.

An object of class `hxsurf` (inherits from `list`) of length 4.

An object of class `hxsurf` of length 4.

An object of class `hxsurf` (inherits from `list`) of length 4.

An object of class `hxsurf` (inherits from `list`) of length 4.

An object of class `hxsurf` (inherits from `list`) of length 4.

An object of class `hxsurf` (inherits from `list`) of length 4.

An object of class `data.frame` with 311 rows and 3 columns.

## See also

`banc.surf` for the available neuropil objects for BANC. These are
`hxsruf` objects, names for subregions can be found as so:
`banc_brain_neuropil.surf$RegionList`

## Examples

``` r
if (FALSE) { # \dontrun{
# Depends on nat
library(nat)
rgl::wire3d(banc.surf)
} # }
```
