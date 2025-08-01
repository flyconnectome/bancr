test_that("banc_scene works", {
  # This test should work as it just generates URLs without network access
  expect_type(sc <- banc_scene("720575941566983162"), 'character')
})
