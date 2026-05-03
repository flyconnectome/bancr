# Query BANC all-to-all influence scores

Retrieve pre-computed influence scores between BANC neurons from
partitioned parquet files stored on Google Cloud Storage. Influence
scores quantify how much a "seed" (upstream) neuron's activity affects a
"target" (downstream) neuron's steady-state response, based on the
connectome's synaptic weight matrix (see Bates et al. 2020).

## Usage

``` r
banc_influence(
  upstream_ids = NULL,
  downstream_ids = NULL,
  const = 24,
  min_score = 0,
  method = c("arrow", "duckdb"),
  local_path = NULL,
  force_download = FALSE
)
```

## Arguments

- upstream_ids:

  Character vector of upstream (seed) neuron root IDs. If `NULL`, all
  upstream neurons are included (use with caution — dataset is very
  large).

- downstream_ids:

  Character vector of downstream (target) neuron root IDs. If `NULL`,
  all downstream neurons are included.

- const:

  Numeric constant for adjusted influence calculation. Adjusted
  influence = `max(0, log(raw_influence) + const)`. Default 24
  corresponds to a minimum meaningful influence of `exp(-24)` (approx
  3.78e-11).

- min_score:

  Minimum adjusted influence score to return. Pairs with adjusted
  influence below this threshold are filtered out. Default 0 returns all
  pairs with non-zero adjusted influence.

- method:

  Character, either `"arrow"` or `"duckdb"`. Controls which backend is
  used to read the parquet files.

- local_path:

  Path to a local directory containing the parquet files. If `NULL`
  (default), uses a cache directory under
  `tools::R_user_dir("bancr", "cache")`.

- force_download:

  Logical. If `TRUE`, re-download parquet files from GCS even if a local
  cache exists. Default `FALSE`.

## Value

A `data.frame` with columns:

- upstream_id:

  Character. Root ID of the upstream (seed) neuron.

- downstream_id:

  Character. Root ID of the downstream (target) neuron.

- raw_influence:

  Numeric. Raw steady-state influence score.

- adjusted_influence:

  Numeric. `max(0, log(raw_influence) + const)`.

## Details

The parquet files contain columns `upstream_id`, `downstream_id`, and
`raw_influence`. Adjusted influence is computed on-the-fly as
`log(raw_influence) + const`, floored at 0.

Data is read from a local cache directory. On first use (or when
`force_download=TRUE`), parquet files are downloaded from GCS using
`gsutil`. Subsequent calls read directly from the cache.

The `"arrow"` method uses
[`open_dataset`](https://arrow.apache.org/docs/r/reference/open_dataset.html)
with predicate pushdown to scan only relevant chunks. The `"duckdb"`
method registers a DuckDB view over the parquet files for fast SQL-based
filtering. For small queries (few upstream/downstream IDs), both methods
are fast; for large scans, `"duckdb"` may be faster.

## References

Bates, A.S., Schlegel, P., Roberts, R.J.V. et al. Complete Connectomic
Reconstruction of Olfactory Projection Neurons in the Fly Brain. *Curr
Biol* 30, 3183-3199.e6 (2020).

## Examples

``` r
if (FALSE) { # \dontrun{
# Get influence of one neuron on all targets
inf <- banc_influence(upstream_ids = "720575941521131930")
head(inf)

# Get influence between specific pairs
inf <- banc_influence(
  upstream_ids = c("720575941521131930", "720575941478275714"),
  downstream_ids = c("720575941555924992")
)

# Use duckdb backend
inf <- banc_influence(
  upstream_ids = "720575941521131930",
  method = "duckdb"
)

# Only return strong connections
inf <- banc_influence(
  upstream_ids = "720575941521131930",
  min_score = 5
)
} # }
```
