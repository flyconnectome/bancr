# Visualise neurons across multiple Drosophila connectomic datasets in BANC spelunker

This function constructs a Neuroglancer scene that visualizes neurons
from multiple co-registered Drosophila connectomic datasets, including
BANC, FAFB, hemibrain, and MANC. It allows for simultaneous
visualization of corresponding neurons across these datasets.

## Usage

``` r
bancsee(
  banc_ids = NULL,
  banc_static_ids = NULL,
  fafb_ids = NULL,
  hemibrain_ids = NULL,
  manc_ids = NULL,
  nuclei_ids = NULL,
  open = FALSE,
  banc.cols = c("#54BCD1", "#0000FF", "#8A2BE2"),
  fafb.cols = c("#C41E3A", "#FF3131", "#F88379"),
  hemibrain.cols = c("#00FF00", "#32CD32", "#006400"),
  hemibrain.mirrored.cols = c("#FFFF00", "#FFD700", "#FFA500"),
  manc.cols = c("#FFA07A", "#FF4500", "#FF8C00"),
  nulcei.col = "#FC6882",
  url = NULL,
  shorturl = TRUE
)
```

## Arguments

- banc_ids:

  A vector of neuron IDs from the BANC dataset. Default is NULL.

- banc_static_ids:

  A vector of neuron IDs from the static v626 BANC release. This dataset
  version represents the snapshot for our preprint. Default is NULL.

- fafb_ids:

  A vector of neuron IDs from the FAFB dataset. Default is NULL.

- hemibrain_ids:

  A vector of neuron IDs from the hemibrain dataset. Default is NULL.

- manc_ids:

  A vector of neuron IDs from the MANC dataset. Default is NULL.

- nuclei_ids:

  A vector of nuclei IDs for the BANC dataset. Default is NULL.

- open:

  Logical; if TRUE, the function will open the Neuroglancer scene in a
  web browser. Default is FALSE.

- banc.cols:

  Vector of hex codes describing a colour spectrum of colours to be
  interpolated for BANC neurons. Defaults are cyan-purple.

- fafb.cols:

  Vector of hex codes describing a colour spectrum of colours to be
  interpolated for BANC neurons. Defaults are red hues.

- hemibrain.cols:

  Vector of hex codes describing a colour spectrum of colors to be
  interpolated for BANC neurons. Defaults green hues.

- hemibrain.mirrored.cols:

  Vector of hex codes describing a colour spectrum of colors to be
  interpolated for BANC neurons. Defaults are yellow hues.

- manc.cols:

  Vector of hex codes describing a colour spectrum of colors to be
  interpolated for MANC neurons. Defaults are orange hues.

- nulcei.col:

  Hex code for the colour in which nuclei will be plotted. Default is
  pink.

- url:

  a spelunker neuroglancer URL.

- shorturl:

  Logical, whether or not to return a shortened URL

## Value

If `open = FALSE`, returns a character string containing the URL for the
Neuroglancer scene. If `open = TRUE`, opens the Spelunker Neuroglancer
scene in a web browser and invisibly returns the URL.

## Details

The function creates a Neuroglancer scene with multiple layers, each
corresponding to a different dataset:

- BANC: "segmentation proofreading" layer

- FAFB: "fafb v783 imported" layer

- Hemibrain: "hemibrain v1.2.1 imported" and "hemibrain v1.2.1 imported,
  mirrored" layers

- MANC: "manc v1.2.1 imported" layer

- BANC nuclei: "nuclei (v1)" layer

Each dataset is assigned a unique color palette to distinguish neurons
from different sources:

- BANC: Blue to purple spectrum

- FAFB: Red spectrum

- Hemibrain: Green spectrum (original) and Yellow spectrum (mirrored)

- MANC: Orange spectrum

- BANC nuclei: Pink

## Note

This function suppresses all warnings during execution. While this
ensures smooth operation, it may hide important messages. Use with
caution and refer to individual function documentation if unexpected
behavior occurs.

## See also

[`banc_scene`](https://flyconnectome.github.io/bancr/reference/banc_scene.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Visualize cell type DNa01 across datasets
bancsee(banc_ids = c("720575941493078142","720575941455137261"),
        fafb_ids = c("720575940644438551","720575940627787609"),
        hemibrain_ids = c("1170939344"),
        manc_ids = c("10751","10760"),
        open = TRUE)

# Get URL without opening browser
url <- bancsee(banc_ids = c("720575941493078142"),
               fafb_ids = c("720575940644438551"),
               open = FALSE)
} # }
```
