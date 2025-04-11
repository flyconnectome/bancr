<!-- badges: start -->
[![natverse](https://img.shields.io/badge/natverse-Part%20of%20the%20natverse-a241b6)](https://natverse.github.io)
[![Docs](https://img.shields.io/badge/docs-100%25-brightgreen.svg)](https://flyconnectome.github.io/bancr/reference/)
[![R-CMD-check](https://github.com/flyconnectome/bancr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/flyconnectome/bancr/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

bancr
===========

The goal of **bancr** is to support analysis of the Brain And
Nerve Cord data set (*BANC*), especially auto-segmentation data. 
Those data are made available by the *BANC* project led by Wei-Chung Allen Lee (Harvard) 
and collaborators including Zetta.ai and the FlyWire team at Princeton. 
Anyone can request access to the data [here](https://flywire.ai/banc_access).
Learn more on the [BANC wiki](https://github.com/jasper-tms/the-BANC-fly-connectome/wiki/).

To access banc resources, you must have permissions to access the [banc
autosegmentation
dataset](https://fanc-reconstruction.slack.com/archives/C01RZP5JH9C/p1616522511001900)
and have [confirmed your
acceptance](https://fanc-reconstruction.slack.com/archives/C01RZP5JH9C/p1617404290005300)
of the banc proofreading and data ownership guidelines. At this point you should
have a linked Google account that will be authorised (see below) for access to
banc online resources.

Broadly speaking the **bancr** package is largely a wrapper over the 
[fafbseg](https://github.com/natverse/fafbseg) package setting up necessary 
default paths etc. It is based on another wrapper for a separate project, 
[fancr](https://github.com/flyconnectome/fancr).

If you have access, you can view *BANC* data in this helpful 
[neuroglancer scene](https://spelunker.cave-explorer.org/#!middleauth+https://global.daf-apis.com/nglstate/api/v1/4753860997414912).

The BANC project uses [CAVE tables](https://global.daf-apis.com/info/) to store many sorts of annotation information. 
You can see the available CAVE tables for BANC [here](https://cave.fanc-fly.com/annotation/views/aligned_volume/brain_and_nerve_cord). CAVE tables can be joined into useful [CAVE views](https://cave.fanc-fly.com/materialize/views/datastack/brain_and_nerve_cord) pinned to a materialisation, which can provide very useful objects such as the whole BANC edgelist. BANC view names and SQL formulas are given [here](https://github.com/jasper-tms/the-BANC-fly-connectome/wiki/CAVE-Views). 

## Installation

You can install the development version of `bancr` from github:

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

# confirm functionality, should return FALSE
banc_islatest("720575941562355975")
```

Some functions rely on underlying Python code by 
[Philipp Schlegel](https://www.zoo.cam.ac.uk/directory/dr-philip-schlegel), 
called using the `reticulate` package.
You can install full set of recommended libraries including `fafbseg-py`:

```
fafbseg::simple_python("full")
```

Note that this package is designed to play nicely with `fafbseg`, which has
been used mainly for the *FAFB-FlyWire* project, but could be used to work with
data from many neuroglancer/CAVE based projects.

If you get an error related to not finding cloud-volume or the 
cloud-volume version, the solution may be to update cloud-volume, as so:

```r
fafbseg::simple_python('none', pkgs='cloud-volume~=8.32.1')
```

Use `with_banc()` to wrap many additional `fafbseg::flywire_*` functions 
for use with the *BANC*. Alternatively `choose_banc()` to set all 
`flywire_*` functions from `fafbseg` to target the *BANC*. Not all functions
will work.

### Updating

You can just repeat the install instructions, but this ensures
that all dependencies are updated:

```r
remotes::install_github('flyconnectome/bancr')
```

If you need to update a specific Python library dependent, you can do:

```r
fafbseg::simple_python(pkgs='fafbseg')
```

Ascending Neuron Vignette
--------------------------------

### Load the code we need

First we need to load the package, and direct ourselves to the *BANC* data set.

```r
library(bancr)
# choose_banc()
```

### Identify the neurons we care about

Next, let us query a *BANC* CAVE table in order to get the neurons users have 
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
an1.left <- "720575941566983162"
an1.right <- "720575941562355975"
```

This ID changes each time a neuron is edited, so while the *BANC* is an 
active project they are unstable. Likely by the time you read this, they
have changed a little, although they describe the same cells.

Therefore, let us make sure we have the most up to date IDs.

```
an1.left <- banc_latestid(an1.left)
an1.right <- banc_latestid(an1.right)
an1.ids <- c(an1.left, an1.right)
an1.ids
```

Sometimes a more stable way to track a neuron (as long as it has a cell body
within the *BANC* volume) is to consider its `nucleus_id`. 

We can get a table of nucleus ids from CAVE and find ours. The `root_id`
column in these CAVE tables automatically update.

```
banc.nuclei <- banc_nuclei()
banc.nuclei.an1 <- banc.nuclei[as.character(banc.nuclei$pt_root_id) %in% an1.ids,]
banc.nuclei.an1.ids <- as.character(banc.nuclei.an1$id)
banc.nuclei.an1
```

### Obtain neuron segmentation data

Great. Next, we want to read the mesh objects of our neurons.

```r
an1.left.mesh <- banc_read_neuron_meshes(an1.left)
an1.right.mesh <- banc_read_neuron_meshes(an1.right)
```

These neurons will be in '*BANC* coordinates', in nanometers. They are read as 
`mesh3d` objects which describe triangular meshes.

But we can also get proxy 'L2' skeletons from the segmentation graph for each
neuron. 

These functions depend on Philipp Schlegel's `fafbseg-py` library. 
You can install this using `fafbseg::simple_python`. See above.

```r
an1.left.skel <- banc_read_l2skel(an1.left)
an1.right.skel <- banc_read_l2skel(an1.right)
```

We can also get `mesh3d` objects for our nuclei.

```r
an1.left.nucleus <- banc_read_nuclei_mesh(banc.nuclei.an1.ids[1])
an1.right.nucleus <- banc_read_nuclei_mesh(banc.nuclei.an1.ids[2])
```

### Plot our *BANC* neurons

We can plot our neurons in 3D using the `rgl` package.

First, we can plot the *BANC* volume mesh which shows all the brain tissue.

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
plot3d(an1.left.mesh, col = "coral", alpha = 0.75)
plot3d(an1.right.mesh, col = "chartreuse", alpha = 0.75)

# Plot neuron skeletons
plot3d(an1.left.skel, col = "darkred", alpha = 1)
plot3d(an1.right.skel, col = "darkgreen", alpha = 1)

# Plot nuclei meshes
plot3d(an1.left.nucleus, col = "black", alpha = 1, add = TRUE)
plot3d(an1.right.nucleus, col = "black", alpha = 1, add = TRUE)
```
![banc_an1](https://github.com/flyconnectome/bancr/blob/main/inst/images/banc_an1.png?raw=true)

We can also make a 2D image of multiple views using `ggplot2`.

```r
# Simplify meshes to  make plotting faster
banc_neuropil <- Rvcg::vcgQEdecim(as.mesh3d(banc_neuropil.surf), percent = 0.1)
banc_brain_neuropil <- Rvcg::vcgQEdecim(as.mesh3d(banc_brain_neuropil.surf), percent = 0.1)
banc_vnc_neuropil <- Rvcg::vcgQEdecim(as.mesh3d(banc_vnc_neuropil.surf), percent = 0.1)
an1.left.mesh.simp <- Rvcg::vcgQEdecim(an1.left.mesh[[1]], percent = 0.1)
an1.right.mesh.simp <- Rvcg::vcgQEdecim(an1.right.mesh[[1]], percent = 0.1)

# Plot! Saves as a PNG file
banc_neuron_comparison_plot(neuron1 = an1.left.mesh.simp,
                           neuron2 = an1.right.mesh.simp,
                           neuron1.info = "AN1_right",
                           neuron2.info = "AN1_left",
                           banc_neuropil = banc_neuropil,
                           banc_brain_neuropil = banc_brain_neuropil,
                           banc_vnc_neuropil = banc_vnc_neuropil,
                           filename = "banc_an_comparison_ggplot2.png?raw=true")

# Tip: You may need to hit 'zoom' on the RStudio plot pane, to see finer meshes,
# when filename = NULL.
```
![banc_an_comparison_ggplot2](https://github.com/flyconnectome/bancr/blob/main/inst/images/banc_an_comparison_ggplot2.png?raw=true)

### Left-right mirror *BANC* neurons

Using a bridge to the symmetric 'template brain' (see below), 
we can 'mirror' neurons in *BANC* even though it is an asymmetric space.

Here we can see the normal (grey) and mirrored mesh (green). At the moment, 
this works less well in the VNC than the brain.

![banc_neuropil_mirrored](https://github.com/flyconnectome/bancr/blob/main/inst/images/banc_neuropil_mirrored.png?raw=true)

```r
an1.left.skel.m <- banc_mirror(an1.left.skel, method = "tpsreg")
an1.right.skel.m <- banc_mirror(an1.right.skel, , method = "tpsreg")
```

And now plot the mirrored skeletons, and the non-mirrored meshes for
comparison:

```r
# Set up 3D plot
nopen3d()
banc_view()
plot3d(banc_neuropil.surf, col = "lightgrey", alpha = 0.1)

# Plot native neuron meshes
plot3d(an1.left.mesh, col = "coral", alpha = 0.5)
plot3d(an1.right.mesh, col = "chartreuse", alpha = 0.5)

# Plot mirrored neuron skeletons
plot3d(an1.left.skel.m, col = "darkred", alpha = 1)
plot3d(an1.right.skel.m, col = "darkgreen", alpha = 1)
```
![banc_ans_mirrored](https://github.com/flyconnectome/bancr/blob/main/inst/images/banc_ans_mirrored.png?raw=true)

We can also change the view to see, for example, the brain more clearly.

```r
banc_front_view()
```
![banc_ans_mirrored_brain](https://github.com/flyconnectome/bancr/blob/main/inst/images/banc_ans_mirrored_brain.png?raw=true)

Or the ventral nerve cord.

```r
banc_vnc_view()
```
![banc_ans_mirrored_vnc](https://github.com/flyconnectome/bancr/blob/main/inst/images/banc_ans_mirrored_vnc.png?raw=true)

### Co-plot *FAFB-FlyWire* and Hemibrain neurons

[Jasper Phelps](https://people.epfl.ch/jasper.phelps/?lang=en) 
has made a *BANC*-to-JRC2018F and JRC2018F-to-*BANC* transform using 
the software [Elastix](https://elastix.dev/download.php). 
Therefore, we can use this transform to move data 
transformed first into `JRC2018F`, into the *BANC*. Or data out of the *BANC*, 
into `JRC2018` and then any other template brain to which `JRC2018F` can be
bridged.

We can either use the Elastix transform directly if you have Elastix installed
on your machine (implemented as `method="elastix"`). 
This can be a bit of a journey, so I have also implemented
a [thin plate spine registration](https://rdrr.io/cran/Morpho/man/computeTransform.html)
that is based on the Elastix transform, 
made and applied using the R package `Morpho` (implemented as `method="tpsreg"`)
. The end result of the two methods can be very slightly different.

Firstly, let's just take the brain part of the ANs we have, as `JRC2018F` only
include the brain.

```r
# Show only the portion in the brain
an1.mesh.simp <- neuronlist(an1.left.mesh.simp, an1.right.mesh.simp)
an1.mesh.simp.brain  <- banc_decapitate(an1.mesh.simp, invert = TRUE)

# Convert to JRC2018F space
an1.mesh.simp.brain.jrc2018f <- banc_to_JRC2018F(an1.mesh.simp.brain,
                                                  method = "tpsreg",
                                                  banc.units = "nm")

# Plot in JRC2018 space
nopen3d()
plot3d(JRC2018F.surf, col = "lightgrey", alpha = 0.1)
plot3d(an1.mesh.simp.brain.jrc2018f, col = c("turquoise", "navy"), alpha = 0.75, add = TRUE)
```
![an_banc_jrc2018f](https://github.com/flyconnectome/bancr/blob/main/inst/images/an_banc_jrc2018f.png?raw=true)

We can now read a neuron from *FAFB-FlyWire*. I already know the ID of the 
comparable *FAFB-FlyWire* neurons to fetch.

We need to load a new R package first.

```r
if (!requireNamespace("remotes")) install.packages("remotes")
remotes::install_github('natverse/nat.jrcbrains')

library(nat.jrcbrains)

## You may need to download the relevant registration, if you have not already:
# download_saalfeldlab_registrations()
```

And then get known *FAFB-FlyWire* neurons.

```r
## if you previously ran choose_banc()
## now run: 
# choose_segmentation("flywire31")
# Which directs you towards the active FAFB-FlyWire segmentation

# Define the IDs we wish to fetch
# these are from the 783 materialisation (i.e.) published version
fw.an1.ids <- c("720575940626768442", "720575940636821616")

# Get neuron meshes
fw.an1.meshes <- read_cloudvolume_meshes(fw.an1.ids)

# Convert to JRC2018F
fw.an1.meshes.jrc2018f <- xform_brain(fw.an1.meshes, sample = "FAFB14",
reference = "JRC2018F")

# Add to plot
plot3d(fw.an1.meshes.jrc2018f, col = c("red","orange"), alpha = 1, add = TRUE)
```
![an_banc_fafb](https://github.com/flyconnectome/bancr/blob/main/inst/images/an_banc_fafb.png?raw=true)

We can do the same with the *Hemibrain*.

We need to load a new R package first.

```r
if (!requireNamespace("remotes")) install.packages("remotes")
remotes::install_github('flyconnectome/hemibrainr')

library(hemibrainr)
``` 

And now we can get the equivalent *Hemibrain* neuron.

```r
# Read hemibrain neuron
hb.an1 <- "706176085"

# Read mesh, divide by 1000 to reach microns
hb.an1.mesh <- hemibrain_neuron_meshes(hb.an1)

# Transforms to JRC2018F, divide by 1000 to reach microns for JRCFIB2018F
hb.an1.mesh.jrc2018f <- xform_brain(hb.an1.mesh/1000, sample = "JRCFIB2018F", reference = "JRC2018F")

# Add to plot
plot3d(hb.an1.mesh.jrc2018f , col = c("chartreuse"), alpha = 1, add = TRUE)
```

![an_banc_fafb_hemibrain](https://github.com/flyconnectome/bancr/blob/main/inst/images/an_banc_fafb_hemibrain.png?raw=true)

Now we see all related neurons from three data sets in one space. Awesome!

We can also see the difference between the Elastix registration and the `Morpho`
based on.

You will first need to download and install Elastix. To do so, you can follow
the instructions [here](https://elastix.dev/download/elastix-5.1.0-manual.pdf). Remember that the Elastix binaries must be on your
`PATH` and your system must be able to see its libraries. On MacOSX it is tricky
to just install and use the [Elastix binaries directly](https://github.com/SuperElastix/elastix/releases), 
you need instead to compile 
ITK then Elastix yourself.

```r
# Transform with Elastix
an1.mesh.simp.brain.jrc2018f.elastix <- banc_to_JRC2018F(an1.mesh.simp.brain,
                                                  method = "elastix",
                                                  banc.units = "nm")

# Plot in JRC2018 space
nopen3d()
plot3d(JRC2018F.surf, col = "lightgrey", alpha = 0.1)
plot3d(an1.mesh.simp.brain.jrc2018f, col = "blue", alpha = 0.75, add = TRUE)
plot3d(an1.mesh.simp.brain.jrc2018f.elastix, col = "green", alpha = 0.75, add = TRUE)
```

### Get neurons connectivity

```r
an.upstream <- banc_partner_summary(an1.ids, partners = "input")
an.downstream <- banc_partner_summary(an1.ids, partners = "output")

# Combine the two data frames and add a source column
combined_data <- bind_rows(
  mutate(an.upstream, source = "Upstream"),
  mutate(an.downstream, source = "Downstream")
) %>%
dplyr::filter(weight>2)

# Create the histogram
ggplot(combined_data, aes(x = weight, fill = source)) +
  geom_histogram(binwidth = 1, color = "black", alpha = 0.7, position = "identity") +
  scale_fill_manual(values = c("Upstream" = "skyblue", "Downstream" = "orange")) +
  labs(title = "histogram of weights: upstream vs downstream",
       x = "Weight",
       y = "Frequency",
       fill = "Source") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "top"
  )
```


Acknowledging the data and tools
--------------------------------

BANC data needs to be acknowledged in accordance to the [BANC community guidelines](https://github.com/jasper-tms/the-BANC-fly-connectome/wiki/)
and in agreement with the BANC consortium. 
If you use this package, you should cite the eventual BANC paper as well (TBD) as our *natverse* paper [(Bates et al. 2020)](https://elifesciences.org/articles/53350)
and the R package as follows:

``` r
citation(package = "bancr")
```

**Bates A, Jefferis G** (2024). _bancr: R Client Access to the Brain And Nerve Cord (BANC) Dataset_. R package version 0.1.0,
  <https://github.com/flyconnectome/bancr>.


Acknowledgements
----------------

The BANC data set was collected at Harvard Medical School in the laboratory of Wei-Chung Allen Lee, by Minsu Kim and Jasper Phelps.
The segmentation and synapse prediction was built by [Zetta.ai](https://zetta.ai/). 
The neuron reconstruction effort has been hosted and supported by [FlyWire](https://flywire.ai/).
This R package was initialised using the [fancr](https://github.com/flyconnectome/fancr) package developed by Greg Jefferis at the MRC Laboratory of Molecular Biology, Cambridge.
Alex Bates worked on this R package while in the laboratory of Rachel Wilson
at Harvard Medical School.

References
----------

**Bates, Alexander Shakeel, James D. Manton, Sridhar R. Jagannathan, Marta Costa, Philipp Schlegel, Torsten Rohlfing, and Gregory SXE Jefferis**. 2020. *The Natverse, a Versatile Toolbox for Combining and Analysing Neuroanatomical Data.* eLife 9 (April). https://doi.org/10.7554/eLife.53350.
