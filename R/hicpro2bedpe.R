#' Convert HiC-Pro results to BEDPE format
#' 
#' @param mat The 3 column sparse upper triangular matrix
#'     from HiC-Pro. 
#' @param bed The BED file containing the mappings for
#'     for the matrix.
#'     
#' @details HiC-Pro will produce a .matrix file
#'     and a .bed file as the final aligned product
#'     of the alignment process. These files should
#'     be read into R using read.table() or similar
#'     and then entered as the mat and bed inputs 
#'     to this function. The function will convert
#'     the data into a format useable for HiCcompare.
#'     The cis matrices in the results can be 
#'     directly input into create.hic.table() as
#'     sparse matrices.
#' @return A list with two items. Item 1, "cis"
#'     contains the intra-chromosomal contact matrices,
#'     one per chromosome.
#'     Item 2, "trans" contains the inter-chromsomal
#'     contact matrix.
#'    
#' @examples 
#' \dontrun{
#'  # read in data
#'  mat <- read.table("hic_1000000.matrix")
#'  bed <- read.table("hic_1000000_abs.bed")
#'  # convert to BEDPE
#'  dat <- hicpro2bedpe(mat, bed)
#' }
#' @export
hicpro2bedpe <- function(mat, bed) {
  # name columns
  colnames(mat) <- c('i', 'j', 'IF')
  colnames(bed) <- c('chr1', 'start1', 'end1', 'id')
  # merge to BEDPE format
  new_mat <- left_join(mat, bed, by = c('i' = 'id'))
  colnames(bed) <- c('chr2', 'start2', 'end2', 'id')
  new_mat <- left_join(new_mat, bed, by = c('j' = 'id'))
  # reorganize columns
  new_mat <- new_mat[, c(4:9, 3)]
  # split between intra & inter matrices
  trans_mat <- subset(new_mat, chr1 != chr2)
  cis_mat <- subset(new_mat, chr1 == chr2)
  # split cis_mat up by chr
  cis_list <- S4Vectors::split(cis_mat, cis_mat$chr1)
  # order list to be natural order
  cis_list <- cis_list[gtools::mixedorder(names(cis_list))]
  results <- list(cis = cis_list, trans = trans_mat)
  return(results)
}