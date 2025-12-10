# Low level access to banc's CAVE annotation infrastructure

Low level access to banc's CAVE annotation infrastructure

## Usage

``` r
banc_cave_client()
```

## Value

A reticulate R object wrapping the python CAVEclient.

## Examples

``` r
if (FALSE) { # \dontrun{
fcc=banc_cave_client()
tables=fcc$annotation$get_tables()
fcc$materialize$get_table_metadata(tables[1])
} # }
```
