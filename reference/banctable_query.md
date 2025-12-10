# Read and write to the seatable for draft BANC annotations

These functions use the logic and wrap some code from the `flytable_.*`
functions in the `fafbseg` R package. `banctable_set_token` will obtain
and store a permanent seatable user-level API token. `banctable_query`
performs a SQL query against a banctable database. You can omit the
`base` argument unless you have tables of the same name in different
bases. `banctable_base` returns a `base` object (equivalent to a mysql
database) which allows you to access one or more tables, logging in to
the service if necessary. The returned base object give you full access
to the Python
[`Base`](https://seatable.github.io/seatable-scripts/python/base/) API
allowing a range of row/column manipulations. `banctable_update_rows`
updates existing rows in a table, returning TRUE on success.
`banctable_move_to_bigdata` moves rows between normal backend and big
data backend. When `invert=FALSE` (archive), it moves all rows from a
specified view to big data storage. When `invert=TRUE` (unarchive), it
moves specific rows by row_id from big data storage back to normal
backend. Note: The big data backend must be enabled in your base for
this function to work.

## Usage

``` r
banctable_query(
  sql = "SELECT * FROM banc_meta",
  limit = 200000L,
  base = NULL,
  python = FALSE,
  convert = TRUE,
  ac = NULL,
  token_name = "BANCTABLE_TOKEN",
  workspace_id = "57832",
  retries = 10,
  table.max = 10000L
)

banctable_set_token(
  user,
  pwd,
  url = "https://cloud.seatable.io/",
  token_name = "BANCTABLE_TOKEN"
)

banctable_login(
  url = "https://cloud.seatable.io/",
  token_name = "BANCTABLE_TOKEN"
)

banctable_update_rows(
  df,
  table,
  base = NULL,
  append_allowed = FALSE,
  chunksize = 1000L,
  workspace_id = "57832",
  token_name = "BANCTABLE_TOKEN",
  ...
)

banctable_move_to_bigdata(
  table = "banc_meta",
  base = "banc_meta",
  url = "https://cloud.seatable.io/",
  workspace_id = "57832",
  token_name = "BANCTABLE_TOKEN",
  view_name = "archive",
  view_id = NULL,
  where = NULL,
  invert = FALSE,
  row_ids = NULL
)

franken_meta(sql = "SELECT * FROM franken_meta", base = "cns_meta", ...)

banctable_append_rows(
  df,
  table,
  bigdata = FALSE,
  base = NULL,
  chunksize = 1000L,
  workspace_id = "57832",
  token_name = "BANCTABLE_TOKEN",
  ...
)
```

## Arguments

- sql:

  A SQL query string. See examples and [seatable
  docs](https://seatable.github.io/seatable-scripts/python/query/).

- limit:

  An optional limit, which only applies if you do not specify a limit
  directly in the `sql` query. By default seatable limits SQL queries to
  100 rows. We increase the limit to 100000 rows by default.

- base:

  Character vector specifying the `base`

- python:

  Logical. Whether to return a Python pandas DataFrame. The default of
  FALSE returns an R data.frame

- convert:

  Expert use only: Whether or not to allow the Python seatable module to
  process raw output from the database. This is is principally for
  debugging purposes. NB this imposes a requirement of seatable_api
  \>=2.4.0.

- ac:

  A seatable connection object as returned by `banctable_login`.

- token_name:

  The name of the token in your .Renviron file, should be
  `BANCTABLE_TOKEN`.

- workspace_id:

  A numeric id specifying the workspace. Advanced use only

- retries:

  if a request to the seatable API fails, the number of times to re-try
  with a 0.1 second pause.

- table.max:

  the maximum number of rows to read from the seatable at one time,
  which is capped at 10000L by seatable.

- user, pwd:

  banctable user and password used by `banctable_set_token` to obtain a
  token

- url:

  Optional URL to the server

- df:

  A data.frame containing the data to upload including an `_id` column
  that can identify each row in the remote table.

- table:

  Character vector specifying a table foe which you want a `base`
  object.

- append_allowed:

  Logical. Whether rows without row identifiers can be appended.

- chunksize:

  To split large requests into smaller ones with max this many rows.

- ...:

  Additional arguments passed to the underlying parallel processing
  functions which might include cl=2 to specify a number of parallel
  jobs to run.

- view_name:

  Character, the name of the view containing rows to archive (required
  for archive operation). Mutually exclusive with view_id.

- view_id:

  Character, the ID of the view containing rows to archive (alternative
  to view_name). Mutually exclusive with view_name.

- where:

  DEPRECATED. The API no longer supports WHERE clauses. Use view_name or
  view_id instead.

- invert:

  Logical. If `FALSE` (default), archives rows from normal backend to
  big data backend (requires view_name or view_id). If `TRUE`,
  unarchives rows from big data backend back to normal backend (requires
  row_ids).

- row_ids:

  Character vector of seatable row IDs. Required for unarchive operation
  (when invert=TRUE). These are the specific rows to move from big data
  backend back to normal backend. Use the table_id (not table_name) for
  unarchive operations.

- bigdata:

  logical, if `TRUE` new rows are added to the bigdata archive rather
  than the 'normal' seatable.

## Value

a `data.frame` of results. There should be 0 rows if no rows matched
query.

## See also

`fafbseg::`[`flytable_query`](https://rdrr.io/pkg/fafbseg/man/flytable-queries.html)

## Examples

``` r
if (FALSE) { # \dontrun{
# Do this once
banctable_set_token(user="MY_EMAIL_FOR_SEATABLE.com",
                    pwd="MY_SEATABLE_PASSWORD",
                    url="https://cloud.seatable.io/")

# Query a table:
banc.meta <- banctable_query()

# Archive rows to big data backend (requires a view):
banctable_move_to_bigdata(
  table = "banc_meta",
  base = "banc_meta",
  view_name = "optic_region_view"
)

# Alternative: use view_id instead of view_name:
banctable_move_to_bigdata(
  table = "banc_meta",
  view_id = "0000"
)

# Unarchive specific rows from big data backend:
banctable_move_to_bigdata(
  table = "banc_meta",
  invert = TRUE,
  row_ids = c("FoDxhChYQSycLm88JZ11RA", "AnotherRowId123")
)
} # }
```
