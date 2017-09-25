#' function to plot distribution of A to get input values for weight_M
#' 
#' @param hic.table A hic.table or list of hic.tables. You must run 
#'     hic_loess on the hic.table before running plot_A.
#' @details This function will plot 4 veiws of the distribution of
#'     A for use in determining the A.min and quant parameters for
#'     the weight_M function. The first plot is the densities of A
#'     for each chromosome cut off at the 75th quantile. The second
#'     plot is the densities of each chromosome cut at the 90th
#'     quantile. The third plot is the log2 histogram of the 
#'     distribution of A. The fourth plot is the histogram of
#'     A cut at the 75th quantile. You should try to choose a
#'     value of A.min where the distribution starts to only
#'     descend after that point. You should choose the value
#'     for quant at the point where the density flattens
#'     out at the lowest values.
#' @return Plots of A.
#' @examples 
#' # Create hic.table object using included Hi-C data in sparse upper
#' # triangular matrix format
#' data("HMEC.chr22")
#' data("NHEK.chr22")
#' hic.table <- create.hic.table(HMEC.chr22, NHEK.chr22, chr= 'chr22')
#' # Plug hic.table into hic_loess()
#' hic.table <- hic_loess(hic.table, Plot = TRUE)
#' # plot_A
#' plot_A(hic.table)
#' 
plot_A <- function(hic.table) {
  # check if list or single hic.table
  if (is.data.table(hic.table)) hic.table = list(hic.table)
  # calculate A
  hic.table <- lapply(hic.table, function(x) x[, A := (adj.IF1 + adj.IF2) / 2])
  # pull out As
  A <- lapply(hic.table, function(x) x$A)
  # make density plots for 75th and 90th quantiles
  A_tmp <- A[[1]]
  A_dens <- density(A_tmp[A_tmp < quantile(A_tmp, 0.75)])
  q75_plot <- ggplot(data.frame(A = A_dens$x, density = A_dens$y), aes(A)) + geom_line(aes(y = density)) + ggtitle("Density of A < 75th Quantile")
  A_dens <- density(A_tmp[A_tmp < quantile(A_tmp, 0.9)])
  q90_plot <- ggplot(data.frame(A = A_dens$x, density = A_dens$y), aes(A)) + geom_line(aes(y = density)) + ggtitle("Density of A < 90th Quantile")
  if (length(A) > 1) {
    for(i in 2:length(A)) {
      A_tmp <- A[[i]]
      A_dens <- density(A_tmp[A_tmp < quantile(A_tmp, 0.75)])
      q75_plot <- q75_plot + geom_line(data = data.frame(A = A_dens$x, density = A_dens$y), aes(y = density))
      A_dens <- density(A_tmp[A_tmp < quantile(A_tmp, 0.9)])
      q90_plot <- q90_plot + geom_line(data = data.frame(A = A_dens$x, density = A_dens$y), aes(y = density))
    }
  }
  # make histograms
  A_tmp <- A[[1]]
  log_hist <- ggplot(data = data.frame(A = log2(A_tmp)), aes(A)) + geom_histogram() + ggtitle("log2(A) Histogram") + labs(x = "Log2(A)")
  A_tmp <- A_tmp[A_tmp < quantile(A_tmp, 0.75)]
  raw_hist <- ggplot(data = data.frame(A = (A_tmp)), aes(A)) + geom_histogram() + ggtitle("A < 75th Quantile Histogram")
  
  # combine plots
  final_plot <- gridExtra::grid.arrange(q75_plot, q90_plot, log_hist, raw_hist, ncol = 2)
  return(final_plot)
}