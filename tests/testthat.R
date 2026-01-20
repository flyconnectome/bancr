library(testthat)
library(bancr)

simple_python("basic")
banc_set_token()
test_check("bancr")
