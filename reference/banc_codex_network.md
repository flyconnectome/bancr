# Generate FlyWireCodex network visualization URLs for BANC v626 dataset

Creates URLs for the FlyWireCodex interactive connectivity browser,
allowing users to visualize neural network diagrams and explore
connectivity patterns in the Brain And Nerve Cord (BANC) connectome. The
BANC dataset represents the first unified brain-and-nerve-cord
connectome of a limbed animal, revealing distributed control
architecture and behaviour-centric neural modules across the entire
central nervous system.

## Usage

``` r
banc_codex_network(
  cell.types = NULL,
  ids = NULL,
  codex.url = "https://codex.flywire.ai/app/connectivity?dataset=banc",
  open = FALSE,
  min_syn_cnt = 3,
  edge_syn_cap = 50
)
```

## Arguments

- cell.types:

  Character vector of cell type names to include in the network. Cell
  types represent functionally and morphologically distinct neuron
  classes in the BANC connectome (e.g., "DNa02" for specific descending
  neurons).

- ids:

  Character vector of neuron root IDs to include in the network. These
  are unique identifiers for individual neurons in the dataset.

- codex.url:

  Character string specifying the base FlyWireCodex URL. Defaults to the
  BANC connectivity browser.

- open:

  Logical indicating whether to automatically open the URL in the
  default web browser. Default is FALSE.

- min_syn_cnt:

  Integer specifying the minimum number of synapses required for
  connections to be displayed. Default is 3.

- edge_syn_cap:

  Integer specifying the maximum number of synapses to display per
  connection for visualization clarity. Default is 50.

## Value

Character string containing the FlyWireCodex URL, or invisible NULL if
`open = TRUE`.

## Details

FlyWireCodex (<https://codex.flywire.ai/?dataset=banc>) provides an
interactive web interface for exploring connectome data. This function
generates properly formatted URLs that pre-configure the visualization
with specified neurons or cell types, enabling researchers to quickly
examine connectivity patterns, synaptic strengths, and network topology.

The BANC v626 dataset contains approximately 160,000 neurons with
complete synaptic connectivity across both brain and ventral nerve cord,
making it ideal for studying sensorimotor integration and distributed
neural computation.

## See also

[`banc_codex_search`](https://flyconnectome.github.io/bancr/reference/banc_codex_search.md)
for generating search URLs,
[`banc_edgelist`](https://flyconnectome.github.io/bancr/reference/banc_cave_tables.md)
for programmatic connectivity data access

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate URL for DNa02 descending neurons and open in browser
banc_codex_network(cell.types = "DNa02", open = TRUE)

# Create URL for multiple cell types with custom synapse threshold
url <- banc_codex_network(cell.types = c("DNa02", "PFL3"), min_syn_cnt = 5)

# Generate URL for specific neuron IDs
banc_codex_network(ids = c("720575941566983162", "720575941562355975"))
} # }
```
