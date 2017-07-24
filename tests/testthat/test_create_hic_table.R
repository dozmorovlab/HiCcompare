test_that('create.hic.table produces same result for sparse matrix and BEDPE inputs', {
  library(testthat)
  library(HiCcompare)
  # test default settings
  data("HMEC.chr22")
  data("NHEK.chr22")
  chr22.table = HiCcompare::create.hic.table(HMEC.chr22, NHEK.chr22, chr = 'chr22')
  h22 = as.data.table(HMEC.chr22)
  n22 = as.data.table(NHEK.chr22)
  bed1 = h22[, ':=' (chr1 = 'chr22', start1 = region1, end1 = region1 + 500000, chr2 = 'chr22', start2 = region2, end2 = region2 + 500000, IF1 = IF)]
  bed1[, ':=' (region1 = NULL, region2 = NULL, IF = NULL)]
  bed2 = n22[, ':=' (chr1 = 'chr22', start1 = region1, end1 = region1 + 500000, chr2 = 'chr22', start2 = region2, end2 = region2 + 500000, IF1 = IF)]
  bed2[, ':=' (region1 = NULL, region2 = NULL, IF = NULL)]
  tab = HiCcompare::create.hic.table(bed1, bed2)
  expect_equal(chr22.table, tab)

  # test include zero
  tab.sparse.include.zero = HiCcompare::create.hic.table(HMEC.chr22, NHEK.chr22, chr = 'chr22', include.zeros = TRUE)
  tab.bed.include.zero = HiCcompare::create.hic.table(bed1, bed2, include.zeros = TRUE)
  expect_equal(tab.sparse.include.zero, tab.bed.include.zero)

  # test scale
  s.scale = HiCcompare::create.hic.table(HMEC.chr22, NHEK.chr22, chr = 'chr22', scale = FALSE)
  b.scale = HiCcompare::create.hic.table(bed1, bed2, scale = FALSE)
  expect_equal(s.scale, b.scale)

  # test subset.dist
  s.subset.dist = HiCcompare::create.hic.table(HMEC.chr22, NHEK.chr22, chr = 'chr22', subset.dist = 100)
  b.subset.dist = HiCcompare::create.hic.table(bed1, bed2, subset.dist = 100)
  expect_equal(s.subset.dist, b.subset.dist)

  # test subset.index
  s.subset.index = HiCcompare::create.hic.table(HMEC.chr22, NHEK.chr22, chr = 'chr22', subset.index = c(1, 50, 1, 50))
  b.subset.index = HiCcompare::create.hic.table(bed1, bed2, subset.index = c(1, 50, 1, 50))
  expect_equal(s.subset.dist, b.subset.dist)
})


test_that('subsetting works', {
  library(HiCcompare)
  library(testthat)
  q = matrix(1:100, 10, 10)
  colnames(q) = 1:10
  w = HiCcompare::full2sparse(q)
  tab = HiCcompare::create.hic.table(w, w, chr = 'test')
  s.dist = HiCcompare::create.hic.table(w, w, chr = 'test', subset.dist = 5)
  # e = sparse2full(s.dist[, c(2,5,7), with=F])

  expect_equal(max(s.dist$D), 5)
  s.index = HiCcompare::create.hic.table(w, w, chr = 'test', subset.index = c(1,5, 1,5))
  e = HiCcompare::sparse2full(s.index[, c(2,5,7), with=F])
  r = q %>% HiCcompare::full2sparse(.) %>% HiCcompare::sparse2full(.)
  expect_identical(e, r[1:5, 1:5])

  s.index2 = HiCcompare::create.hic.table(w, w, chr='test', subset.index = c(5, 9, 5, 9))
  y = HiCcompare::sparse2full(s.index2[, c(2,5,7), with=F])
  expect_identical(y, r[5:9, 5:9])
})


test_that('Input errors are correct', {
  library(HiCcompare)
  library(testthat)
  data("HMEC.chr22")
  data("nhek.IS")
  expect_error(HiCcompare::create.hic.table(HMEC.chr22, nhek.IS, chr = 'chr22'), "Make sure the classes of the sparse matrices match")
  expect_error(HiCcompare::create.hic.table(nhek.IS, HMEC.chr22, chr = 'chr22'), "Make sure the classes of the sparse matrices match")
  data("NHEK.chr22")
  chr22.table = HiCcompare::create.hic.table(HMEC.chr22, NHEK.chr22, chr = 'chr22')
  NHEK.chr22_BEDPE <- chr22.table[, c(1:6, 8), with=FALSE]
  expect_error(HiCcompare::create.hic.table(NHEK.chr22_BEDPE, HMEC.chr22, chr='chr22'), "Enter both sparse matrices in the same format; either 7 column BEDPE or 3 column sparse upper triangular matrix")
})
