test_that('hic_loess errors work', {
  library(HiCdiff)
  data('HMEC.chr22')
  data('NHEK.chr22')
  tab = create.hic.table(HMEC.chr22, NHEK.chr22, chr = 'chr22')
  expect_error(hic_loess(tab, degree = 5), "'degree' must be 0, 1 or 2")
  expect_error(hic_loess(tab, span = 0), "Enter a larger value for span")
  expect_error(hic_loess(tab, span = 2), "Enter a value <= 1 for span")
  expect_error(hic_loess(tab, loess.criterion = 'gkv'))
  expect_error(hic_loess(tab, diff.thresh = 'a'))
})


test_that('hic_diff errors work', {
  library(HiCdiff)
  data('HMEC.chr22')
  data('NHEK.chr22')
  tab = create.hic.table(HMEC.chr22, NHEK.chr22, chr = 'chr22')
  expect_error(hic_diff(tab))
  tab = hic_loess(tab)
  expect_error(hic_diff(tab, iterations = 50))
  expect_error(hic_diff(tab, diff.thresh = 'gwo'))
  expect_error(hic_diff(tab, diff.thresh = 0))
})
