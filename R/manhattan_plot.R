#' Create a Manhattan plot for the results of HiCcompare
#'
#' @param hic.table a hic.table object that has been normalized
#'     and had differences detected.
#' @param adj.p Logical, should the adjusted p-value be used (TRUE)
#'     of the raw p-value (FALSE)?
#' @param alpha The alpha level for calling a p-value significant.
#' @param return_df Logical, should the data.frame built to be used
#'    for plotting be returned? If TRUE then the data.frame will be
#'    returned and the plot will only be printed.
#' @details This function will produce a manhattan plot of the results
#'     of hic_compare(). Can be used to display which regions around found
#'     to be significantly different on the linear genome.
#' @return A manhattan plot.
#' @examples 
#' # Create hic.table object using included Hi-C data in 
#' # sparse upper triangular matrix format
#' data('HMEC.chr22')
#' data('NHEK.chr22')
#' hic.table <- create.hic.table(HMEC.chr22, NHEK.chr22, chr = 'chr22')
#' # Plug hic.table into hic_loess()
#' result <- hic_loess(hic.table, Plot = TRUE)
#' # perform difference detection
#' diff.result <- hic_compare(result, Plot = TRUE)
#' # make manhattan plot
#' manhattan_plot(diff.result)



manhattan_plot <- function(hic.table, adj.p = TRUE, alpha = 0.05, return_df = FALSE) {
  # get regions on chr
  regions <- unique(c(hic.table$start1, hic.table$start2))
  # count number of times each region found significant
  if (adj.p) {
    num_sig <- sapply(regions, function(x) sum(hic.table$start1 == x & hic.table$p.adj < alpha) + sum(hic.table$start2 == x & hic.table$p.adj < alpha))
  } else {
    num_sig <- sapply(regions, function(x) sum(hic.table$start1 == x & hic.table$p.value < alpha) + sum(hic.table$start2 == x & hic.table$p.value < alpha))
  }
  # make data.frame for plotting
  df <- data.frame(chr = hic.table$chr1[1], bp = regions, count = num_sig)
  # make plot
  p1 <- ggplot(df, aes(x = bp, y = count)) + geom_point() + labs(y = 'Number of times significant', x = paste0(df$chr[1], ' position'))
  if (return_df) {
    print(p1)
    return(df)
  } else {
    return(p1)
  }
}   