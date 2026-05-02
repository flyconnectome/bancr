library(testthat)
library(bancr)
if (identical(Sys.getenv("BANCR_ON_CI"), "true")) {
  # On CI we use the Python provided by the workflow (RETICULATE_PYTHON)
  # rather than letting fafbseg install miniconda. With miniconda=FALSE,
  # simple_python uses the existing interpreter and pip-installs the
  # required Python packages into it.
  simple_python("basic", miniconda = FALSE)
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
