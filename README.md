# bancr

<!-- badges: start -->
[![natverse](https://img.shields.io/badge/natverse-Part%20of%20the%20natverse-a241b6)](https://natverse.github.io)
[![Docs](https://img.shields.io/badge/docs-100%25-brightgreen.svg)](https://flyconnectome.github.io/bancr/reference/)
[![R-CMD-check](https://github.com/flyconnectome/banc/workflows/R-CMD-check/badge.svg)](https://github.com/flyconnectome/banc/actions)
[![R-CMD-check](https://github.com/flyconnectome/bancr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/flyconnectome/bancr/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of **bancr** is to support analysis of the Brain And
Nerve Cord dataset aka (BANC), especially autosegmentation data. 
Those data are made available by the BANC project led by Wei-Chung Allen Lee (Harvard) 
and collaborators including Zetta.ai and the FlyWire team at Princeton. 

To access banc resources, you must have permissions to access the [banc
autosegmentation
dataset](https://fanc-reconstruction.slack.com/archives/C01RZP5JH9C/p1616522511001900)
and have [confirmed your
acceptance](https://fanc-reconstruction.slack.com/archives/C01RZP5JH9C/p1617404290005300)
of the banc proofreading and data ownership guidelines. At this point you should
have a linked Google account that will be authorised (see below) for access to
banc online resources.

Broadly speaking the **bancr** package is a thin wrapper over the 
[fafbseg](https://github.com/natverse/fafbseg) package setting up necessary 
default paths etc. It is based on another wrapper for a separate project, 
[fancr](https://github.com/flyconnectome/fancr).

## Installation

You can install the development version of bancr from github:

```r
remotes::install_github('flyconnectome/bancr')
```

To do anything useful with the bancr package, you need authorisation to access
banc resources. To prove your authorisation for programmatic access you must
generate and store a token in your web browser after logging in to an approved
Google account. This should be streamlined by running the following command in R
(which will also set you up for Pythonic access via cloudvolume).

```r
# set up token - will open your browser to generate a new token
banc_set_token()


# if you already have one do 
# banc_set_token("<my token>")
```

To check that everything is set up properly, try:

```r
# diagnose issues
dr_banc()

# confirm functionality
banc_xyz2id(cbind(34495, 82783, 1954), rawcoords=TRUE)
svids=banc_leaves("720575941650432785")
head(svids)
```

Some functions rely on underlying Python code by Philipp Schlegel, 
called using the `reticulate` package.
You can install full set of recommended libraries including `fafbseg-py`:

```
simple_python("full")
```

### Updating

You can just repeat the install instructions, but this ensures
that all dependencies are updated:

```r
remotes::install_github('flyconnectome/bancr')
```

If you need to update a specific Python library dependent, you can do:

```r

```

### Vignette

## Plotting related ascending neurons

### Load the code we need

First we need to load the package, and direct ourselves to the BANC data set.

```r
library(bancr)
choose_banc()
```

### Identify the neurons we care about

Next, let's query a BANC CAVE table in order to get the neurons users have 
annotated as 'ascending' neurons, i.e. neurons that have their cell bodies
and dendrites in the ventral nerve cord, and their axons in the brain.

```r
banc.neck.connective.neurons <- banc_neck_connective_neurons()
head(banc.neck.connective.neurons)
```

After considering these neurons, I have decided I would like to 
plot two of them. They are both members of the same cell type. They are 
identified by a 16-bit `root_id`.

```
an1.right <- "720575941566983162"
an1.left <- "720575941562355975"
```

This ID changes each time a neuron is edited, so while the BANC is an 
active project they are unstable. Likely by the time you read this, they
have changed a little, although they describe the same cells.

Therefore, let us make sure we have the most up to date IDs.

```
an1.right <- banc_latestid(an1.right)
an1.left <- banc_latestid(an1.left)
an1.ids <- c(an1.right, an1.left)
an1.ids
```

Sometimes a more stable way to track a neuron (as long as it has a cell body
within the BANC volume) is to consider its `nucleus_id`. 

We can get a table of nucleus ids from CAVE and find ours. The `root_id`
column in these CAVE tables automatically update.

```
banc.nulcei <- banc_nuclei()
banc.nuclei.an1 <- subset(banc.nulcei, banc.nulcei$pt_root_id %in% an1.ids)
banc.nuclei.an1.ids <- as.character(banc.nuclei.an1$id)
banc.nuclei.an1
```

### Obtain neuron segmentation data

Great. Next, we want to read the mesh objects of our neurons.

```r
an1.right.mesh <- banc_read_neuron_meshes(an1.right)
an1.left.mesh <- banc_read_neuron_meshes(an1.left)
```

These neurons will be in 'BANC coordinates', in nanometers. They are read as 
`mesh3d` objects which describe triangular meshes.

But we can also get proxy 'L2' skeletons from the segmentation graph for each
neuron. 

These functions depend on Philipp Schlegel's `fafbseg-py` library. 
You can install this using `fafbseg::simple_python`. See above.

```r
an1.right.skel <- banc_read_l2skel(an1.right)
an1.left.skel <- banc_read_l2skel(an1.left)
```

We can also get `mesh3d` objects for our nuclei.

```r
an1.right.nucleus <- banc_read_nuclei_mesh(an1.right)
an1.left.nucleus <- banc_read_nuclei_mesh(an1.left)
```

### Plot our BANC neurons

We can plot our neurons in 3D using the `rgl` package.

First, we can plot the BANC volume mesh which shows all the brain tissue.

```r
nopen3d()
banc_view()
plot3d(banc.surf, col = "lightgrey", alpha = 0.1)
```

We can also see the synaptic neuropil inside of it.

```r
plot3d(banc_neuropil.surf, col = "lightgrey", alpha = 0.25)
```

And now our neurons, their skeletons and their nuclei.

```r
# Plot neuron meshes
plot3d(an1.right.mesh, col = "coral", alpha = 0.75)
plot3d(an1.left.mesh, col = "chartreuse", alpha = 0.75)

# Plot neuron skeletons
plot3d(an1.right.skel, col = "darkred", alpha = 1)
plot3d(an1.left.skel, col = "darkgreen", alpha = 1)

# Plot nuclei meshes
plot3d(an1.right.nucleus, col = "darkred", alpha = 0.75)
plot3d(an1.left.nucleus, col = "darkgreen", alpha = 0.75)
```

We can also make a 2D image of multiple views using `ggplot2`.

```r
# Simplify meshes to  make plotting faster
banc_neuropil <- Rvcg::vcgQEdecim(as.mesh3d(banc_neuropil.surf), percent = 0.1)
banc_brain_neuropil <- Rvcg::vcgQEdecim(as.mesh3d(banc_brain_neuropil.surf), percent = 0.1)
banc_vnc_neuropil <- Rvcg::vcgQEdecim(as.mesh3d(banc_vnc_neuropil.surf), percent = 0.1)
an1.right.mesh.simp <- Rvcg::vcgQEdecim(an1.right.mesh[[1]], percent = 0.1)
an1.left.mesh.simp <- Rvcg::vcgQEdecim(an1.left.mesh[[1]], percent = 0.1)

# Pot!
g <- banc_neuron_comparison_plot(neuron1 = an1.right.mesh.simp,
                                 neuron2 = an1.left.mesh.simp,
                                 neuron1.info = "AN1_right",
                                 neuron2.info = "AN1_left",
                                 banc_neuropil = banc_neuropil,
                                 banc_brain_neuropil = banc_brain_neuropil,
                                 banc_vnc_neuropil = banc_vnc_neuropil)

# Tip: You may need to hit 'zoom' on the RStudio plot pane, to see finer meshes.
plot(g)
```

### Left-right mirror BANC neurons

Using a bridge to the symmetric templatebrain (see below), 
we can 'mirror' neurons in BANC even though it is an unsymmetric space.

```r
an1.right.skel.m <- banc_mirror(an1.right.skel, method = "tpsreg")
an1.left.skel.m <- banc_mirror(an1.left.skel, , method = "tpsreg")
```

And now plot the mirrored skeletons, and the non-mirrored meshes for
comparison:

```r
# Set up 3D plot
nopen3d()
banc_view()
plot3d(banc_neuropil.surf, col = "lightgrey", alpha = 0.1)

# Plot native neuron meshes
plot3d(an1.right.mesh, col = "coral", alpha = 0.3)
plot3d(an1.left.mesh, col = "chartreuse", alpha = 0.3)

# Plot mirrored neuron skeletons
plot3d(an1.right.skel.m, col = "darkred", alpha = 1)
plot3d(an1.left.skel.m, col = "darkgreen", alpha = 1)
```

### Co-plot FAFB-FlyWire and Hemibrain neurons

Jasper Phelps has made a BANC-to-JRC2018F and JRC2018F-to-BANC transform using 
the software Elastix. Therefore, we can use this transform to move data 
transformed first into `JRC2018F`, into the BANC. Or data out of the BANC, 
into `JRC2018` and then any other template brain to which `JRC2018F` can be
bridged.

We can either use the Elastix transform directly if you have Elastix installed
on your machine (implemented as `method="elastix"`). 
This can be a bit of a journey, so I have also implemented
a thinplate-spine registration that is based on the Elastix transform, 
made and applied using the R package `Morpho` (implemented as `method="tpsreg"`)
. The end result of the two methods can be very slightly different.

```r
```

```r
```

We can also see the difference between the Elastix registration and the `Morpho`
based on.

```r
```


