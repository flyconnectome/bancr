library(testthat)
library(bancr)
if (identical(Sys.getenv("BANCR_ON_CI"), "true")) {
  # Only run this on CI
  simple_python("basic")
  bancr::banc_set_token()
  bancr::dr_banc()
}
test_check("bancr")
