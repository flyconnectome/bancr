skip_if_not_installed('coconatfly')

test_that("coconatfly integration works", {
  register_banc_coconat()
  library(coconatfly)
  # suppress warnings as there is one about parsing ngl scenes
  bm=try(suppressWarnings(banc_meta()))
  skip_if(inherits(bm, 'try-error'),
          message = 'Skipping coconatfly tests as failed to read banc metadata!')

  expect_s3_class(dna02meta <- cf_meta(cf_ids(banc='/DNa02')), 'data.frame')
  expect_true(nrow(dna02meta)==2)

  expect_s3_class(dna02ds <- cf_partner_summary(
      dna02meta, partners = 'out',threshold = 10),
    'data.frame')

  expect_s3_class(
    cf_cosine_plot(cf_ids(banc='/type:DNa.+'), partners = 'o', heatmap = F),
    'hclust')
})
