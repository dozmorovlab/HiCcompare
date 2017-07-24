#' Perform MA normalization on a hic.table object
#'
#' @export
#' @param hic.table A hic.table object
#' @param degree The degree for loess normalization
#' @param Plot logical, should the MA plot be output?
#' @param span The span for loess. If left as the default value of NA the span will be calculated automatically
#' @param loess.criterion The criterion for calculating the span for loess
#'
#' @details Performs loess normalization on the MA plot of the data.
#' @return An extended hic.table with adjusted IFs and M columns.
#'
#' @examples
#' # create hic.table
#' data("HMEC.chr22")
#' data("NHEK.chr22")
#' hic.table <- create.hic.table(HMEC.chr22, NHEK.chr22, chr= 'chr22')
#' # Plug hic.table into MA_norm()
#' MA_norm(hic.table)

MA_norm <- function(hic.table, degree = 2, Plot = FALSE, span = NA,
                    loess.criterion = "gcv") {
  A <- 0.5 * log2(hic.table$IF1 * hic.table$IF2)
  # perform loess on data
  # if (is.na(span)) {
  #   l <- .loess.as(x = hic.table$D, y = hic.table$M, degree = degree,
  #                  criterion = loess.criterion,
  #                  control = loess.control(surface = "interpolate",
  #                                          statistics = "approximate", trace.hat = "approximate"))
  #   # calculate gcv and AIC
  #   traceL <- l$trace.hat
  #   sigma2 <- sum(l$residuals^2)/(l$n - 1)
  #   aicc <- log(sigma2) + 1 + 2 * (2 * (traceL + 1))/(l$n - traceL -
  #                                                       2)
  #   gcv <- l$n * sigma2/(l$n - traceL)^2
  # } else {
  #   l <- loess(hic.table$M ~ A, degree = degree, span = span,
  #              control = loess.control(surface = "interpolate", statistics = "approximate",
  #                                      trace.hat = "approximate"))
  #   # calculate gcv and AIC
  #   traceL <- l$trace.hat
  #   sigma2 <- sum(l$residuals^2)/(l$n - 1)
  #   aicc <- log(sigma2) + 1 + 2 * (2 * (traceL + 1))/(l$n - traceL -
  #                                                       2)
  #   gcv <- l$n * sigma2/(l$n - traceL)^2
  # }
  # # print the span picked by gcv
  # message(paste("Span for loess: ", l$pars$span, sep = ""))
  # message(paste("GCV for loess: ", gcv, sep = ""))
  # message(paste("AIC for loess: ", aicc, sep = ""))

  l <- loess(hic.table$M ~ A)
  # get the correction factor
  mc <- predict(l, A)
  if (Plot) {
    MDplot <- data.frame(D = A, M = hic.table$M, mc = mc)
    plot1 <- ggplot(MDplot, aes(x = D, y = M)) + geom_point() +
      geom_abline(slope = 0, intercept = 0) + stat_smooth(geom = "smooth", size = 1) +
      labs(title = "Before loess",  x = "A")
    plot2 <- ggplot(MDplot, aes(x = D, y = M - mc)) + geom_point() +
      geom_abline(slope = 0, intercept = 0) + labs(title = "After loess", x = "A") +
      stat_smooth(geom = "smooth", size = 1) + ylab("")
    gridExtra::grid.arrange(plot1, plot2, ncol = 2)
  }
  # create mhat matrix using mc/2 which will be subtracted/added to the
  # original matrices to produce the loess normalized matrices
  mhat <- mc/2
  # normalize original matrices
  if (sum(hic.table$IF1 == 0) + sum(hic.table$IF2 == 0) > 0) {
    hic.table[, `:=`(adj.IF1, 2^(log2(IF1 + 1) + mhat))]
    hic.table[, `:=`(adj.IF2, 2^(log2(IF2 + 1) - mhat))]
    hic.table[, `:=`(adj.M, log2((adj.IF2)/(adj.IF1)))]
    hic.table[, `:=`(mc, mc)]
  } else {
    hic.table[, `:=`(adj.IF1, 2^(log2(IF1) + mhat))]
    hic.table[, `:=`(adj.IF2, 2^(log2(IF2) - mhat))]
    hic.table[, `:=`(adj.M, log2((adj.IF2)/(adj.IF1)))]
    hic.table[, `:=`(mc, mc)]
  }
  hic.table[, A := A]
  return(hic.table)
}
