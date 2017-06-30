
test_that('matrix transformations work', {
  library(HiCcompare)
  data('HMEC.chr22')
  original = as.data.table(HMEC.chr22)
  full = sparse2full(HMEC.chr22)
  sparse = full2sparse(full)
  expect_equal(original, sparse)
})
