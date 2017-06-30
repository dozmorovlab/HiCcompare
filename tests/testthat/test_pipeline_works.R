

test_that('Same output when seed set and input the same', {
  library(HiCdiff)
  data('HMEC.chr22')
  data('NHEK.chr22')
  tab = create.hic.table(HMEC.chr22, NHEK.chr22, chr = 'chr22')
  set.seed(1)
  norm1 = hic_loess(tab, check.differences = T)
  set.seed(1)
  norm2 = hic_loess(tab, check.differences = T)
  expect_equal(norm1, norm2)
})


test_that('Pipeline works', {
  library(HiCdiff)
  data('HMEC.chr22')
  data('NHEK.chr22')
  tab = create.hic.table(HMEC.chr22, NHEK.chr22, chr = 'chr22')
  set.seed(1)
  norm1 = hic_loess(tab, check.differences = F)
  norm2 = hic_diff(norm1)
  norm3 = hic_diff(norm1, diff.thresh = 1, Plot = F)
  norm4 = hic_diff(norm1, diff.thresh = 0.6, Plot =F)
  norm5 = hic_diff(norm1, diff.thresh = NA, Plot = F)
  expect_equal(class(norm1)[1], "data.table")
  expect_equal(class(norm2)[1], "data.table")
  expect_equal(class(norm3)[1], "data.table")
  expect_equal(class(norm4)[1], "data.table")
  expect_equal(class(norm5)[1], "data.table")
})


# test_that('Settings for hic_loess work', {
#   library(HiCdiff)
#   data('HMEC.chr22')
#   data('NHEK.chr22')
#   tab = create.hic.table(HMEC.chr22, NHEK.chr22, chr = 'chr22')
#   norm = hic_loess(tab, span = 0.05)
#   #norm = hic_loess(tab, span = 1)
#   #norm = hic_loess(tab, span = 2, Plot = T)
#   norm = hic_loess(tab, loess.criterion = 'aicc', Plot= T)
#
# })




