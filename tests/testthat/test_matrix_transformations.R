
test_that('matrix transformations work', {
  library(HiCcompare)
  data('HMEC.chr22')
  original = as.data.table(HMEC.chr22)
  full = HiCcompare::sparse2full(HMEC.chr22)
  sparse = HiCcompare::full2sparse(full)
  expect_equal(original, sparse)
})
