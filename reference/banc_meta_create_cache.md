# Create or refresh cache of BANC meta information

`banc_meta_create_cache()` builds or refreshes an in-memory cache of
BANC metadata for efficient repeated lookups. You can choose the data
source using `use_seatable`. The main accessor function
[`banc_meta()`](https://flyconnectome.github.io/bancr/reference/banc_meta.md)
will always use the most recently created cache.

## Usage

``` r
banc_meta_create_cache(use_seatable = FALSE, return = FALSE)
```

## Arguments

- use_seatable:

  Whether to build BANC meta data from the `codex_annotations` CAVE
  table (production) or our internal seatable (development). Both
  require different types of authenticated access, for details see
  `bancr` documentation.

- return:

  Logical; if `TRUE`, return the cache tibble/invisible.

## Value

Invisibly returns the cache (data.frame) if `return=TRUE`; otherwise
invisibly `NULL`.

## Details

BANC meta queries can be slow; caching avoids repeated
computation/database access. Whenever labels are updated, simply rerun
this function to update the cache.

## See also

Other coconatfly:
[`register_banc_coconat()`](https://flyconnectome.github.io/bancr/reference/register_banc_coconat.md)

## Examples

``` r
if (FALSE) { # \dontrun{
#' # Requires authenticated access to BANC CAVE
banc_meta_cache(use_seatable=FALSE)

banc_meta_create_cache(use_seatable=TRUE) # create cache
## BANCTABLE_TOKEN must be set, see bancr package
result <- banc_meta() # use cache

# use cache to quickly make plot
library(coconatfly)
# only needed once per session
register_banc_coconat()
cf_cosine_plot(cf_ids('/type:LAL0(08|09|10|42)', datasets = c("banc", "hemibrain")))
} # }
```
