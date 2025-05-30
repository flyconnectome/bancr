% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggplot2.R
\name{banc_neuron_comparison_plot}
\alias{banc_neuron_comparison_plot}
\title{Compare two neurons from the BANC connectome dataset}
\usage{
banc_neuron_comparison_plot(
  neuron1 = NULL,
  neuron2 = NULL,
  neuron3 = NULL,
  neuron1.info = NULL,
  neuron2.info = NULL,
  neuron3.info = NULL,
  volume = NULL,
  region = c("both", "brain", "vnc"),
  banc_brain_neuropil = NULL,
  banc_vnc_neuropil = NULL,
  banc_neuropil = NULL,
  alpha = 0.5,
  filename = NULL,
  width = 16,
  height = 16
)
}
\arguments{
\item{neuron1}{A neuronlist object representing the first neuron to be compared.}

\item{neuron2}{A neuronlist object representing the second neuron to be compared.}

\item{neuron3}{A neuronlist object representing the third neuron to be compared.}

\item{neuron1.info}{Character, a vector to be printed on the outputted ggplot in reference to neuron1.}

\item{neuron2.info}{Character, a vector to be printed on the outputted ggplot in reference to neuron2.}

\item{neuron3.info}{Character, a vector to be printed on the outputted ggplot in reference to neuron3.}

\item{volume}{A mesh3d  or hxsurf object in BANC space that you wish you co-plot.}

\item{region}{Character, whether to plot the brain area, VNC area or both (default).}

\item{banc_brain_neuropil}{A mesh object representing the brain neuropil. Default is banc_brain_neuropil.surf.}

\item{banc_vnc_neuropil}{A mesh object representing the ventral nerve cord (VNC) neuropil. Default is banc_vnc_neuropil.surf.}

\item{banc_neuropil}{A mesh object representing the entire neuropil. Default is banc_neuropil.surf.}

\item{alpha}{Vector of alpha values for neurons 1, 2 and 3. If a single value, applied to all.}

\item{filename}{Optional. If provided, the plot will be saved to this file. If NULL (default), the plot will be displayed but not saved.}

\item{width}{Numeric. The width of the output plot in inches. Default is 16.}

\item{height}{Numeric. The height of the output plot in inches. Default is 16.}
}
\value{
If filename is NULL, the function returns a ggplot object. If filename is provided, the function saves the plot and returns invisibly.
}
\description{
This function creates a visual comparison of two neurons from the BANC (Bilateral Antennal lobe Neuron Connectome)
dataset. It generates four different anatomical views to facilitate comparison.
}
\details{
The function generates four views of the neurons:
\enumerate{
\item Main view
\item Side view
\item Front view (brain only)
\item VNC (Ventral Nerve Cord) view
}

Each view applies appropriate rotations and uses different neuropil meshes as backgrounds.
The two neurons are plotted in different colors (navy/turquoise for neuron1, red/darkred for neuron2) for easy comparison.
}
\examples{
\dontrun{
# Assuming neuron1 and neuron2 are valid neuron objects
banc.meta <- banctable_query()

# Get some neurons to plot
banc.meta.dnao1 <- subset(banc.meta, cell_type=="DNa01")
dna01 <- banc_read_neuron_meshes(unique(banc.meta.dnao1$root_id))
banc.meta.dnao2 <- subset(banc.meta, cell_type=="DNa02")
dna02 <- banc_read_neuron_meshes(unique(banc.meta.dnao2$root_id))

# Simplify neurons to make them easier to plot
dna01 <- nat::nlapply(dna01,Rvcg::vcgQEdecim,percent = 0.1)
dna02 <- nat::nlapply(dna02,Rvcg::vcgQEdecim,percent = 0.1)

# Make plot!
banc_neuron_comparison_plot(dna01, dna02, neuron1.info = "DNa01",
neuron2.info = "DNa02", filename = "neuron_comparison.png")
}

}
