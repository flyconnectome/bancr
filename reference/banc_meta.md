# Query cached BANC meta data

Returns results from the in-memory cache, filtered by `ids` if given.
Cache must be created first using
[`banc_meta_create_cache()`](https://flyconnectome.github.io/bancr/reference/banc_meta_create_cache.md).

## Usage

``` r
banc_meta(ids = NULL)
```

## Arguments

- ids:

  Vector of neuron/root IDs to select, or `NULL` for all.

## Value

tibble/data.frame, possibly filtered by ids.

## Details

`banc_meta()` never queries databases directly. If `ids` are given,
filters the meta table by root_id.

## See also

[`banc_meta_create_cache()`](https://flyconnectome.github.io/bancr/reference/banc_meta_create_cache.md)

## Examples

``` r
if (FALSE) { # \dontrun{
banc_meta_create_cache() # build the cache
all_meta <- banc_meta()  # retrieve all
} # }
```
