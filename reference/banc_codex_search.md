# Generate FlyWireCodex search URLs for BANC v626 dataset

Creates URLs for the FlyWireCodex search interface, allowing users to
search and browse neuron metadata in the Brain And Nerve Cord (BANC)
connectome. This function enables researchers to quickly access detailed
information about specific cell types or individual neurons, including
morphological data, cell type classifications, and connectivity
summaries. The BANC dataset represents the first unified
brain-and-nerve-cord connectome of a limbed animal, providing
unprecedented insight into distributed neural control.

## Usage

``` r
banc_codex_search(
  cell.types = NULL,
  ids = NULL,
  codex.url = "https://codex.flywire.ai/app/search?dataset=banc",
  open = FALSE,
  page.size = 100
)
```

## Arguments

- cell.types:

  Character vector of cell type names to search for. Cell types
  represent functionally and morphologically distinct neuron classes in
  the BANC connectome (e.g., "DNa02" for specific descending neurones).

- ids:

  Character vector of neuron root IDs to search for. These are unique
  identifiers for individual neurons in the dataset.

- codex.url:

  Character string specifying the base FlyWireCodex search URL. Defaults
  to the BANC search interface.

- open:

  Logical indicating whether to automatically open the URL in the
  default web browser. Default is FALSE.

- page.size:

  Integer specifying the number of search results to display per page.
  Default is 100.

## Value

Character string containing the FlyWireCodex search URL, or invisible
NULL if `open = TRUE`.

## Details

FlyWireCodex (<https://codex.flywire.ai/?dataset=banc>) provides an
interactive web interface for exploring connectome data. The search
function allows researchers to query the comprehensive metadata
associated with BANC neurons, including cell type annotations,
morphological measurements, and connectivity statistics.

The BANC v626 dataset contains detailed annotations for approximately
160,000 neurons across brain and ventral nerve cord regions, with rich
metadata supporting studies of sensorimotor integration, distributed
computation, and behavior-centric neural modules.

## See also

[`banc_codex_network`](https://flyconnectome.github.io/bancr/reference/banc_codex_network.md)
for generating network visualization URLs,
[`banc_codex_annotations`](https://flyconnectome.github.io/bancr/reference/banc_codex_annotations.md)
for programmatic metadata access

## Examples

``` r
if (FALSE) { # \dontrun{
# Search for DNa02 descending neurons and open in browser
banc_codex_search(cell.types = "DNa02", open = TRUE)

# Create search URL for multiple cell types with custom page size
url <- banc_codex_search(cell.types = c("DNa02", "PFL3"), page.size = 50)

# Generate search URL for specific neuron IDs
banc_codex_search(ids = c("720575941566983162", "720575941562355975"))
} # }
```
