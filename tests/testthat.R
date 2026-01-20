library(testthat)
library(bancr)
if (identical(Sys.getenv("BANCR_ON_CI"), "true")) {
  # Only run this on CI
  simple_python("basic")

  # Use the path where your workflow writes the secret
  token_path <- path.expand("~/.cloudvolume/secrets/chunkedgraph-secret.json")
  if (file.exists(token_path)) {
    try(bancr::banc_set_token(token_path), silent = TRUE)
  }
}
test_check("bancr")
