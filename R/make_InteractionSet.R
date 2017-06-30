#' Convert HiCdiff results to InteractionSet object
#'
#' @export
#' @param hic.table A hic.table object.
#'
#' @details This function will convert data from HiCdiff results in the
#'     hic.table object format to the InteractionSet format which makes
#'     use of GRanges objects.
#'
#' @return An object of class InteractionSet
#' @examples
#' # create hic.table
#' data(HMEC.chr22)
#' data(NHEK.chr22)
#' hic.table <- create.hic.table(HMEC.chr22, NHEK.chr22, chr='chr22')
#' # convert to InteractionSet
#' gi <- make_InteractionSet(hic.table)


make_InteractionSet <- function(hic.table) {
  # calculate bin size
  bins <- unique(c(hic.table$start1, hic.table$start2))
  bins <- bins[order(bins)]
  bin.size <- min(diff(bins))
  # make GRanges objects
  gr1 <- GRanges(seqnames = hic.table$chr1,
                 ranges = IRanges(start = hic.table$start1,
                                  end = hic.table$start1 + bin.size - 1))
  gr2 <- GRanges(seqnames = hic.table$chr2,
                 ranges = IRanges(start = hic.table$start2,
                                  end = hic.table$start2 + bin.size - 1))
  gi <- GInteractions(gr1, gr2)
  S4Vectors::values(gi) <- cbind(S4Vectors::values(gi), hic.table[, 7:ncol(hic.table), with = FALSE])
  return(gi)
}
