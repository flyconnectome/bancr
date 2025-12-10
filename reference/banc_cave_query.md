# Query banc tables in the CAVE annotation system

Query banc tables in the CAVE annotation system

## Usage

``` r
banc_cave_query(table, live = TRUE, ...)
```

## Arguments

- table:

  The name of the table (or view, see views section) to query

- live:

  Whether to use live query mode, which updates any root ids to their
  current value (or to another `timestamp` when provided). Values of
  `TRUE` or `1` select CAVE's *Live* mode, while `2` selects `Live live`
  mode which gives access even to annotations that are not part of a
  materialisation version. See section **Live and Live Live queries**
  for details.

- ...:

  Additional arguments passed to
  [`flywire_cave_query`](https://rdrr.io/pkg/fafbseg/man/flywire_cave_query.html)
  [`flywire_cave_query`](https://rdrr.io/pkg/fafbseg/man/flywire_cave_query.html)

## Value

A data.frame

## See also

[`flywire_cave_query`](https://rdrr.io/pkg/fafbseg/man/flywire_cave_query.html)

## Examples

``` r
if (FALSE) { # \dontrun{
library(dplyr)
cell_info=banc_cave_query('cell_info')
cell_info %>%
  filter(tag2=='anterior-posterior projection pattern') %>%
  count(tag)
} # }
```
