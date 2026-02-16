library(testthat)
library(bancr)
if (identical(Sys.getenv("BANCR_ON_CI"), "true")) {
  # Only run this on CI
  simple_python("basic")
  # Token file is written by the GitHub Actions workflow;
  # only call banc_set_token() if it is missing
  secret_path <- file.path(Sys.getenv("HOME"), ".cloudvolume", "secrets", "chunkedgraph-secret.json")
  if (!file.exists(secret_path) || file.size(secret_path) == 0) {
    message("No chunkedgraph secret found at: ", secret_path)
  }
  tryCatch(bancr::dr_banc(), error = function(e) {
    message("dr_banc() failed: ", e$message)
  })
}
test_check("bancr")
