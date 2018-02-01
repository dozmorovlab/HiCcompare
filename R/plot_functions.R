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
#' @param smooth Should smooth scatter plots be used? If set to FALSE
#'     ggplot scatter plots will be used. When option is TRUE, plots
#'     generate quicker. It is recommend to use the smooth scatter
#'     plots.
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
MD.plot1 <- function(M, D, mc, smooth = TRUE) {
  MDplot <- data.frame(D = D, M = M, mc = mc)
  if (smooth) {
    # using ggplot
    # plot1 <- ggplot(MDplot, aes(x = D, y = M)) +
    #   stat_density2d(aes(fill = ..density..^0.25), geom = "tile", contour = FALSE, n = 200) +
    #   scale_fill_continuous(low = "white", high = "dodgerblue4") +
    #   geom_point(alpha = 0.01, shape = 20) + theme_classic() +
    #   geom_abline(slope = 0, intercept = 0) + stat_smooth(geom = "smooth", size = 1, se = FALSE) +
    #   labs(title = "Before loess",  x = "Distance") +
    #   theme(legend.position = "none")
    #
    # plot2 <- ggplot(MDplot, aes(x = D, y = M - mc)) +
    #   stat_density2d(aes(fill = ..density..^0.25), geom = "tile", contour = FALSE, n = 200) +
    #   scale_fill_continuous(low = "white", high = "dodgerblue4") +
    #   geom_point(alpha = 0.01, shape = 20) + theme_classic() +
    #   geom_abline(slope = 0, intercept = 0) + labs(title = "After loess", x = "Distance") +
    #   stat_smooth(geom = "smooth", size = 1, se = FALSE) + ylab("") +
    #   theme(legend.position = "none")
    # gridExtra::grid.arrange(plot1, plot2, ncol = 2)

    # using smoothScatter
    par(mfrow = c(1,2), mai = c(0.8, 0.759, 0.5, 0.1))
    smoothScatter(D, M, xlab = 'Distance', ylab = 'M', main = 'Before loess')
    lines(D[order(D)], mc[order(D)], col = 'red')
    abline(h = 0)
    smoothScatter(D, M - mc, ylab = '', xlab = 'Distance', main = 'After loess')
    abline(h = 0)
    .reset_par()
  } else {
    plot1 <- ggplot(MDplot, aes(x = D, y = M)) + geom_point() +
      geom_abline(slope = 0, intercept = 0) + stat_smooth(geom = "smooth", size = 1, se = FALSE) +
      labs(title = "Before loess",  x = "Distance")
    plot2 <- ggplot(MDplot, aes(x = D, y = M - mc)) + geom_point() +
      geom_abline(slope = 0, intercept = 0) + labs(title = "After loess", x = "Distance") +
      stat_smooth(geom = "smooth", size = 1, se = FALSE) + ylab("")

    gridExtra::grid.arrange(plot1, plot2, ncol = 2)
  }

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
#' @param smooth Should smooth scatter plots be used? If set to FALSE
#'     ggplot scatter plots will be used. When option is TRUE, plots
#'     generate quicker. It is recommend to use the smooth scatter
#'     plots.
#' @return An MD plot.
#'
#' @examples
#' data('HMEC.chr22')
#' data('NHEK.chr22')
#' hic.table <- create.hic.table(HMEC.chr22, NHEK.chr22, chr='chr22')
#' # Plug hic.table into hic_loess()
#' result <- hic_loess(hic.table)
#' # perform difference detection
#' diff.result <- hic_compare(result)
#' MD.plot2(diff.result$M, diff.result$D, diff.result$p.value)
#'
MD.plot2 <- function(M, D, p.val = NA, diff.thresh = NA, smooth = TRUE) {
  # smooth scatter version
  if (smooth) {
    smoothScatter(D, M, xlab = 'Distance', ylab = 'M', main = 'MD Plot')
    abline(h = 0)
    if (!is.na(p.val[1])) {
      p0.01 <- which(p.val < 0.01)
      p0.05 <- which(p.val >= 0.01 & p.val < 0.05)
      points(D[p0.01], M[p0.01], col = "red", pch = 20)
      points(D[p0.05], M[p0.05], col = 'yellow', pch = 20)
      legend('bottomright', legend = c('P < 0.01', 'P < 0.05'), fill = c('red', 'yellow'), bty = 'n', horiz = TRUE)
    }
    if (!is.na(diff.thresh)) {
      abline(h = diff.thresh, lty = 2)
      abline(h = -diff.thresh, lty = 2)
    }
  } else {
    # ggplot version
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
}


# function to reset par
.reset_par <- function(){
  op <- structure(list(xlog = FALSE, ylog = FALSE, adj = 0.5, ann = TRUE,
                       ask = FALSE, bg = "transparent", bty = "o", cex = 1, cex.axis = 1,
                       cex.lab = 1, cex.main = 1.2, cex.sub = 1, col = "black",
                       col.axis = "black", col.lab = "black", col.main = "black",
                       col.sub = "black", crt = 0, err = 0L, family = "", fg = "black",
                       fig = c(0, 1, 0, 1), fin = c(6.99999895833333, 6.99999895833333
                       ), font = 1L, font.axis = 1L, font.lab = 1L, font.main = 2L,
                       font.sub = 1L, lab = c(5L, 5L, 7L), las = 0L, lend = "round",
                       lheight = 1, ljoin = "round", lmitre = 10, lty = "solid",
                       lwd = 1, mai = c(1.02, 0.82, 0.82, 0.42), mar = c(5.1, 4.1,
                                                                         4.1, 2.1), mex = 1, mfcol = c(1L, 1L), mfg = c(1L, 1L, 1L,
                                                                                                                        1L), mfrow = c(1L, 1L), mgp = c(3, 1, 0), mkh = 0.001, new = FALSE,
                       oma = c(0, 0, 0, 0), omd = c(0, 1, 0, 1), omi = c(0, 0, 0,
                                                                         0), pch = 1L, pin = c(5.75999895833333, 5.15999895833333),
                       plt = c(0.117142874574832, 0.939999991071427, 0.145714307397962,
                               0.882857125425167), ps = 12L, pty = "m", smo = 1, srt = 0,
                       tck = NA_real_, tcl = -0.5, usr = c(0.568, 1.432, 0.568,
                                                           1.432), xaxp = c(0.6, 1.4, 4), xaxs = "r", xaxt = "s", xpd = FALSE,
                       yaxp = c(0.6, 1.4, 4), yaxs = "r", yaxt = "s", ylbias = 0.2), .Names = c("xlog",
                                                                                                "ylog", "adj", "ann", "ask", "bg", "bty", "cex", "cex.axis",
                                                                                                "cex.lab", "cex.main", "cex.sub", "col", "col.axis", "col.lab",
                                                                                                "col.main", "col.sub", "crt", "err", "family", "fg", "fig", "fin",
                                                                                                "font", "font.axis", "font.lab", "font.main", "font.sub", "lab",
                                                                                                "las", "lend", "lheight", "ljoin", "lmitre", "lty", "lwd", "mai",
                                                                                                "mar", "mex", "mfcol", "mfg", "mfrow", "mgp", "mkh", "new", "oma",
                                                                                                "omd", "omi", "pch", "pin", "plt", "ps", "pty", "smo", "srt",
                                                                                                "tck", "tcl", "usr", "xaxp", "xaxs", "xaxt", "xpd", "yaxp", "yaxs",
                                                                                                "yaxt", "ylbias"))
  par(op)
}
