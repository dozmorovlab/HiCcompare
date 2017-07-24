# Plotting functions for hic_loess

#' Visualize the MD plot before and after loess normalization
#'
#' @export
#' @param M The M component of the MD plot. Available from the
#'     hic.table object.
#' @param D The D component of the MD plot. The unit distance of
#'     the interaction. Available from the hic.table object.
#' @param mc The correction factor. Calculated by \code{hic_loess}.
#'     Available from the hic.table object after running \code{hic_loess}.
#'
#' @return An MD plot.
#' @examples
#' # Create hic.table object using included Hi-C data in sparse upper
#' # triangular matrix format
#' data('HMEC.chr22')
#' data('NHEK.chr22')
#' hic.table <- create.hic.table(HMEC.chr22, NHEK.chr22, chr = 'chr22')
#' # Plug hic.table into hic_loess()
#' result <- hic_loess(hic.table)
#' MD.plot1(result$M, result$D, result$mc)
#'
MD.plot1 <- function(M, D, mc) {
  MDplot <- data.frame(D = D, M = M, mc = mc)
  plot1 <- ggplot(MDplot, aes(x = D, y = M)) + geom_point() +
    geom_abline(slope = 0, intercept = 0) + stat_smooth(geom = "smooth", size = 1) +
    labs(title = "Before loess",  x = "Distance")
  plot2 <- ggplot(MDplot, aes(x = D, y = M - mc)) + geom_point() +
    geom_abline(slope = 0, intercept = 0) + labs(title = "After loess", x = "Distance") +
    stat_smooth(geom = "smooth", size = 1) + ylab("")
  gridExtra::grid.arrange(plot1, plot2, ncol = 2)
}


#' Visualize the MD plot.
#'
#' @export
#' @param M The M component of the MD plot. Available from the
#'     hic.table object.
#' @param D The D component of the MD plot. The unit distance of the
#'     interaction. Available from the hic.table object.
#' @param p.val An optional p-value vector to provide color to the plot
#'     based on the significance of the differences between the IFs.
#' @param diff.thresh A difference threshold used for calculating p-values.
#'     If set to a value will add dotted horizontal lines to the plot to
#'     display the threshold cutoffs. See `hic_loess` or `hic_diff`
#'     functions help for more information on this parameter.
#'
#' @return An MD plot.
#'
#' @examples
#' data('HMEC.chr22')
#' data('NHEK.chr22')
#' hic.table <- create.hic.table(HMEC.chr22, NHEK.chr22, chr='chr22')
#' # Plug hic.table into hic_loess()
#' result <- hic_loess(hic.table)
#' # perform difference detection
#' diff.result <- hic_diff(result, diff.thresh = 1)
#' MD.plot2(diff.result$M, diff.result$D, diff.result$p.value)
#'
MD.plot2 <- function(M, D, p.val = NA, diff.thresh = NA) {
  if (is.na(p.val[1])) {
    MAplot <- data.frame(M = M, D = D)
    plot1 <- ggplot(MAplot, aes(x = D, y = M)) + geom_point(colour = "black") +
      geom_abline(slope = 0, intercept = 0) + labs(title = "MD Plot") +
      theme(legend.position = "top", legend.title = element_blank(),
            legend.direction = "horizontal", legend.margin = margin(t = 0,
                                                                    b = 0, unit = "cm"))
  } else {
    p.colors <- ifelse(p.val < 0.01, "p < 0.01", ifelse(p.val < 0.05,
                                                        "p < 0.05", "p >= 0.05"))
    group.colors <- c("p < 0.01" = "#F8766D", "p < 0.05" = "#00BA38", "p >= 0.05" = "#619CFF")
    MAplot <- data.frame(M = M, D = D, p.colors = factor(p.colors))
    plot1 <- ggplot(MAplot, aes(x = D, y = M, color = p.colors)) +
      geom_point() + geom_abline(slope = 0, intercept = 0) + labs(title = "MD Plot") +
      theme(legend.position = "top", legend.title = element_blank(),
            legend.direction = "horizontal", legend.margin = margin(t = 0,
                                                                    b = 0, unit = "cm")) +
      scale_color_manual(values = group.colors)
    if (!is.na(diff.thresh)) {
      plot1 <- plot1 + geom_hline(yintercept = diff.thresh, linetype = "dashed") +
        geom_hline(yintercept = -diff.thresh, linetype = "dashed")
    }
  }
  return(plot1)
}
