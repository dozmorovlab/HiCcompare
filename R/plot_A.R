# function to plot distribution of A to get input values for weight_M
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