# Download all of the BANC synapses as a .sqlite file that you can read lazily from later

Download all of the BANC synapses as a .sqlite file that you can read
lazily from later

## Usage

``` r
banc_all_synapses(
  path =
    c("gs://lee-lab_brain-and-nerve-cord-fly-connectome/synapses/v2.0/final_edgelist.csv",
    "gs://zetta_lee_fly_cns_001_synapse/240623_run/assignment/final_edgelist.df"),
  overwrite = FALSE,
  n_max = 2000,
  details = FALSE,
  min_size = 10,
  rawcoords = FALSE
)
```

## Arguments

- path:

  The google storage path to the desired synapses file. Read using
  [`readr::read_csv`](https://readr.tidyverse.org/reference/read_delim.html).

- overwrite:

  Logical, whether or not to overwrite an extant `banc_data.sqlite`
  file.

- n_max:

  Numeric, the maximum number of rows to read from `path` if you just
  want to see a taster of the file.

- details:

  Logical Whether or not to read all data columns in the target synapse
  `.csv`. Defaults to `FALSE` in order to read only the essential
  presynapse position data.

- min_size:

  Numeric, filter parameter, the minimum size (in nm) of the detected
  synapse.

- rawcoords:

  Logical, whether or not to convert from raw coordinates into
  nanometers. Default is `FALSE`.

## Value

a data.frame

## Details

Downloads all automatic Zetta.ai synapse detections for the BANC and
saves them as a `banc_data.sqlite` file. Once this is done, in the
future the function will read from this file lazily so as not to throw
the whole thing into system memory.

## See also

[`banc_partner_summary`](https://flyconnectome.github.io/bancr/reference/banc_partner_summary.md),
[`banc_partners`](https://flyconnectome.github.io/bancr/reference/banc_partner_summary.md)

## Examples

``` r
if (FALSE) { # \dontrun{
syns <- banc_all_synapses()
} # }
```
