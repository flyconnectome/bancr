# Resolve local path to influence parquet files, downloading from GCS if needed

Resolve local path to influence parquet files, downloading from GCS if
needed

## Usage

``` r
banc_influence_path(local_path = NULL, force_download = FALSE, gs_url = NULL)
```

## Arguments

- local_path:

  Local directory for cached parquet files. Default uses
  `tools::R_user_dir("bancr", "cache")`.

- force_download:

  If TRUE, re-download even if local files exist.

- gs_url:

  GCS path to the influence parquet directory. Defaults to the v850
  all-to-all influence data.
