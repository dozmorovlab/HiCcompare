# - We need heatmap visualization of Hi-C matrices (like Fig 1A), overlapping the p-values of dots identified as differentially expressed. The goal is to see whether HiCcompare-detected differential interactions follow some patterns on the HiC heatmaps, e.g., overlap with visible TAD boundaries.  
# - Differentially interacting "hotspots" - individual regions that are frequently differentially interactiong with others
# - Are they cluster? There may be stretches of differential hotspots


#' Function to visualize p-values from HiCcompare results
#' 
#' @param hic.table A hic.table object that has been
#'     normalized and has had differences detected.
#' @details The goal of this function is to visualize
#'     where in the Hi-C matrix the differences are
#'     occuring between two experimental conditions.
#'     The function will produce a heatmap of the
#'     -log10(p-values) * sign(adj.M) 
#'     to visualize where the
#'     significant differences between the datasets
#'     are occuring on the genome. 
#' @return A heatmap
#' @examples 
#' # Create hic.table object using included Hi-C data in sparse upper triangular
#' # matrix format
#' data('HMEC.chr22')
#' data('NHEK.chr22')
#' hic.table <- create.hic.table(HMEC.chr22, NHEK.chr22, chr = 'chr22')
#' # Plug hic.table into hic_loess()
#' result <- hic_loess(hic.table, Plot = TRUE)
#' # perform difference detection
#' diff.result <- hic_diff(result, diff.thresh = 'auto', Plot = TRUE)
#' # visualize p-values
#' visualize_pvals(diff.result)
#' @export

visualize_pvals <- function(hic.table) {
  # check that hic.table is entered
  if (!is(hic.table, "data.table")) {
    stop("Enter a hic.table object")
  }
  # check that hic.table has necessary columns
  if (ncol(hic.table) < 16) {
    stop('Enter a hic.table that you have already performed normalization and difference detection on')
  }
  # convert to full matrix
  m <- sparse2full(hic.table, hic.table = TRUE, column.name = 'p.value') 
  # get fold change
  fc <- sparse2full(hic.table, hic.table = TRUE, column.name = 'adj.M')
  # plot heatmap
  pheatmap::pheatmap(-log10(m) * sign(fc), cluster_rows = FALSE, cluster_cols = FALSE)
}