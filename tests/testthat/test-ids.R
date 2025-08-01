test_that("banc_xyz2id works", {

  skip('Skipping banc_xyz2id test due to missing supervoxel infrastructure!')
  expect_equal(banc_xyz2id(cbind(34495, 82783, 1954), rawcoords=TRUE),
               "648518346499897667")

  expect_equal(
    banc_xyz2id(cbind(34495, 82783, 1954), rawcoords=TRUE, root=F),
    "73186243730767724")
})

test_that("banc_islatest works", {
  skip_if_not_installed("reticulate")
  skip_if(!reticulate::py_module_available("caveclient"), "caveclient not available")
  skip_on_ci() # Skip on continuous integration 
  
  expect_false(banc_islatest("720575941480769421"))
  expect_false(isTRUE(all.equal(
    banc_latestid("720575941480769421"), "720575941480769421")))
})

test_that("banc_ids works", {
  expect_equal(banc_ids("720575941480769421"), "720575941480769421")
  expect_equal(banc_ids("720575941480769421", integer64 = T), bit64::as.integer64("720575941480769421"))

  df1=data.frame(pt_root_id=bit64::as.integer64("720575941480769421"))
  df2=data.frame(id=bit64::as.integer64("720575941480769421"))

  expect_equal(banc_ids(df1, integer64 = F), "720575941480769421")
  expect_equal(banc_ids(df1), df1$pt_root_id)
  expect_equal(banc_ids(df2, integer64 = F), "720575941480769421")
})

test_that("banc_cellid_from_segid", {
  skip("Skipping banc_cellid_from_segid as BANC doesn't yet have a proper cell_id table")
  rid=banc_latestid("720575941480769421")
  expect_equal(banc_cellid_from_segid(rid),12967L)
})
