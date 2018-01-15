#' Function to visualize p-values from HiCcompare results
#' 
#' @param hic.table A hic.table object that has been
#'     normalized and has had differences detected.
#' @param alpha The alpha level at which you will 
#'     call a p-value significant. If this is set to
#'     a numeric value then any p-values >= alpha will
#'     be set to 1 for the visualization in the heatmap.
#'     Defaults to NA for visualization of all p-values.
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

visualize_pvals <- function(hic.table, alpha = NA) {
  # check that hic.table is entered
  if (!is(hic.table, "data.table")) {
    stop("Enter a hic.table object")
  }
  # check that hic.table has necessary columns
  if (ncol(hic.table) < 16) {
    stop('Enter a hic.table that you have already performed normalization and difference detection on')
  }
  if (!is.na(alpha)) {
    if (!is.numeric(alpha)) {
      stop("Alpha must be a numeric value between 0 and 1")
    }
    if (alpha < 0 | alpha > 1) {
      stop("Alpha must be a numeric value between 0 and 1")
    }
  }
  # convert to full matrix
  m <- sparse2full(hic.table, hic.table = TRUE, column.name = 'p.value') 
  # get fold change
  fc <- sparse2full(hic.table, hic.table = TRUE, column.name = 'adj.M')
  # remove non significant values from matrix if alpha is set to a value
  if (!is.na(alpha)) {
    m[m >= alpha] <- 1
  }
  # plot heatmap
  pheatmap::pheatmap(-log10(m) * sign(fc), cluster_rows = FALSE, 
                     cluster_cols = FALSE, show_rownames = FALSE, 
                     show_colnames = FALSE)
  # # get TADs
  # m2 <- sparse2full(hic.table, hic.table = TRUE, column.name = 'IF1') 
  # tads <- HiCseg::HiCseg_linkC_R(size_mat = nrow(m2), nb_change_max = round(nrow(m2) / 5) + 1, distrib = "B", mat_data = m2, model = "Dplus")
  # image(1:nrow(m2), 1:nrow(m2), -log10(m)*sign(fc), xlab="", ylab="")
  # t_hat=c(1,tads$t_hat[tads$t_hat!=0]+1)
  # for (i in 1:(length(t_hat)-1)) {
  #   lines(c(t_hat[i],t_hat[i]),c(t_hat[i],(t_hat[(i+1)]-1)))
  #   lines(c(t_hat[(i+1)]-1,t_hat[(i+1)]-1),c(t_hat[i],t_hat[(i+1)]-1))
  #   lines(c(t_hat[i],t_hat[(i+1)]-1),c(t_hat[i],t_hat[i]))
  #   lines(c(t_hat[i],t_hat[(i+1)]-1),c(t_hat[(i+1)]-1,t_hat[(i+1)]-1))
  # }
  # 
  # library(ComplexHeatmap)
  # mat <- -log10(m) * sign(fc)
  # Heatmap(mat, name = "pval", cluster_rows = FALSE, cluster_columns = FALSE)
  # decorate_heatmap_body("pval", {
  #   # i = t_hat[i]
  #   # x = i/ncol(mat)
  #   # grid.lines(c(x,x), c(0,1), gp=gpar(lwd=2))
  #   for (i in 1:(length(t_hat)-1)) {
  #     i = t_hat[i]
  #     x = i/ncol(mat)
  #     y = t_hat[(i+1)]-1
  #     y = y/ ncol(mat)
  #     grid.lines(x = c(x,x), y = c(x,y))
  #     grid.lines(x = c(y,y), y = c(x,y))
  #     grid.lines(x = c(x,y), y = c(x,x))
  #     grid.lines(x = c(x,y), y = c(y,y))
  #   }
  # })
  
}


#' #' Alternate heatmap for dealing with rankings
#' #' @param which_rank The column name for the ranking that you
#' #'    want to plot.
#' #' @param only_toprank Logical, Should only the top ranks be 
#' #'    plotted? If TRUE then the proportion option will be used
#' #'    to determine how many ranks will be considered top.
#' #' @param proportion The proportion of ranks to be considered
#' #'    as top ranks. Defaults to 0.05.
#' visualize_differences <- function(hic.table, which_rank = 'rnkMean', only_toprank = TRUE, proportion = 0.05) {
#'   # check that hic.table is entered
#'   if (!is(hic.table, "data.table")) {
#'     stop("Enter a hic.table object")
#'   }
#'   # check that hic.table has necessary columns
#'   if (ncol(hic.table) < 16) {
#'     stop('Enter a hic.table that you have already performed normalization and difference detection on')
#'   }
#'   # sort hic.table for top rankings
#'   hic.table <- hic.table[order(get(which_rank)), ]
#'   idx <- 1:(proportion * nrow(hic.table))
#'   # make indicator for top ranks
#'   tRank <- rep(0, nrow(hic.table))
#'   tRank[idx] <- 1
#'   hic.table[, topRank := tRank]
#'   # convert to full matrix
#'   m <- sparse2full(hic.table, hic.table = TRUE, column.name = which_rank) 
#'   # normalize ranks by total length of ranks
#'   m <- m / nrow(hic.table)
#'   # get fold change
#'   fc <- sparse2full(hic.table, hic.table = TRUE, column.name = 'adj.M')
#'   if (only_toprank) {
#'     # top rank indicator matrix
#'     tr <- sparse2full(hic.table, hic.table = TRUE, column.name = 'topRank')
#'     # set non-top ranks to NA
#'     m <- m * tr
#'   }
#'   # make 0 values 1 since we are -log10 transforming m
#'   m[m == 0] <- 1
#'   
#'   # plot heatmap
#'   pheatmap::pheatmap(-log10(m) * sign(fc), cluster_rows = FALSE, 
#'                      cluster_cols = FALSE, show_rownames = FALSE, 
#'                      show_colnames = FALSE)
#'  
#'   
#' }
