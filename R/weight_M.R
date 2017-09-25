#' Function to weight M values on the MD plot by Average expression
#' 
#' @param hic.table A hic.table object or hic.table list. 
#'     You must run hic_loess on the hic.table before performing M weighting.
#' @param A.min The A value at which to start the powerlaw estimation of the
#'     density. This should be determined by looking at the plots of A from
#'     the plot_A function. See vignette for more details. Defaults to 1.
#' @param quant The quantile value at which to stop weighting M values.
#'     Typical values should range from 0.6 to 0.9. This can be determined
#'     by looking at the plots of the distribution of A from the plot_A
#'     function. See vignette for more details. Defaults to 0.75.
#' @param Plot.diagnostic Logical, should a histogram of the subset of 
#'     A between A.min and A < quant be plotted with the powerlaw fit
#'     overlaid. It is recommend to check the fit of the powerlaw to 
#'     the distribution. If a good fit is not achieved you may need
#'     to adjust A.min or quant parameters. Defaults to TRUE.
#' @param Plot.MD Logical, should the MD plot of the weighted M 
#'     values be plotted?  Defaults to TRUE.
#' @details Weighting of M is used as a way to account for the fact
#'     that IFs with low average values are less trustworthy compared
#'     to IFs with high average values. i.e. a comparison between 
#'     IF1 = 1 and IF2 = 2 results in an M value of 1, but a 
#'     comparison of IF1 = 100 and IF2 = 200 also resuts in an
#'     M value of 1. However the average expression of the first
#'     comparison is 1.5 while the second is 150. Naturally the
#'     first comparison is much more likely to have occured as a 
#'     fluke compared to a true biological difference between the
#'     conditions being tested. This function performs down-weighting
#'     of M values which arise from IFs exhibiting low average
#'     expression. For a more in-depth explanation with examples
#'     please see the vignette.
#' @return A hic.table with weighted adj.M values.
#' @examples 
#' # Create hic.table object using included Hi-C data in sparse upper
#' # triangular matrix format
#' data("HMEC.chr22")
#' data("NHEK.chr22")
#' hic.table <- create.hic.table(HMEC.chr22, NHEK.chr22, chr= 'chr22')
#' # Plug hic.table into hic_loess()
#' hic.table <- hic_loess(hic.table, Plot = TRUE)
#' # weight M-values
#' hic.table <- weight_M(hic.table, A.min = 75, quant = 0.9)
#'
weight_M <- function(hic.table, A.min = 1, quant = 0.75, Plot.diagnostic = TRUE, Plot.MD = TRUE) {
  # check for correct input
  if (is(hic.table, "list")) {
    if ( sapply(hic.table, ncol) %>% min() < 13) {
      stop("Make sure you run hic_loess() on your hic.table before inputting it into hic_diff()")
    }
  } else {
    if (ncol(hic.table) < 13) {
      stop("Make sure you run hic_loess() on your hic.table before inputting it into hic_diff()")
    }
  }
  if (A.min < 1) {
    stop("A.min should be a numeric value >= 1")
  }
  if (quant < 0 | quant > 1) {
    stop("quant should be a numeric value in [0,1]")
  }
  # check if list or single hic.table
  if (is.data.table(hic.table)) hic.table = list(hic.table)
  
  # Run weighting function
  hic.table <- lapply(hic.table, .internal_weight_M, A.min = A.min, quant = quant, Plot.diagnostic = Plot.diagnostic, Plot.MD = Plot.MD)
  
  # clean up if single hic.table
  if (length(hic.table) == 1) {
    hic.table <- hic.table[[1]]
  }
  return(hic.table)
}


# internal function to applied over list of hic.tables
.internal_weight_M <- function(hic.table, A.min = 1, quant = 0.75, Plot.diagnostic = TRUE, Plot.MD = TRUE) {
  # calculate A
  hic.table[, A := (adj.IF1 + adj.IF2) / 2]
  # get specified quantile for A
  A_quant <- quantile(hic.table$A, quant)
  # # get density of A between 1 and A_quant
  density_A <- density(hic.table$A[hic.table$A > A.min & hic.table$A < A_quant], cut = 0)
  # get PDF of density_A
  A_pdf <- approxfun(density_A)
  
  # fit powerlaw curve to PDF
  resid_plaw <- function(A, par) {
    resid <- A_pdf(A) - .plaw(A, par[1], par[2])
    SSR <- sum(resid^2, na.rm = TRUE)
    return(SSR)
  }
  # minimize the residuals between powerlaw and density of A 
  est <- optim(par = c(0,1), resid_plaw, A = hic.table$A)
  
  # plot histogram of A and powerlaw fit to it
  if (Plot.diagnostic) {
    A_seq <- seq(A.min, max(hic.table$A[hic.table$A > A.min & hic.table$A < A_quant]), by = 0.001)
    hist(hic.table$A[hic.table$A > A.min & hic.table$A < A_quant], freq = FALSE, main = 'Powerlaw fit to A', xlab = 'A')
    lines(.plaw(A_seq, est$par[1], est$par[2]) ~ A_seq)
  }
  
  
  # calculate weight - weight is 1 - (f(A) / maximum value f(A)) where f(A) is the powerlaw fit
  # 1 - because f of A is decreasing
  # divide by max(f(A)) to make the weights stay on (0, 1)
  # set weights for A < peak of A_pdf to 0 and for A > quant_A to 1
  weight_func <- function(A) {
    peak_A <- A.min
    peak_f <- max(.plaw(A.min, est$par[1], est$par[2]))
    message("M values corresponding to Average expression < ", peak_A, " being weighted to 0.")
    weight <-ifelse(A < peak_A, 0, ifelse(A > A_quant, 1, 
                                             1 - (.plaw(A, est$par[1], est$par[2]) / peak_f)))
    return(weight)
  }
  weight <- weight_func(hic.table$A)
  # check weights are between 0 and 1
  if (sum(weight < 0) > 1 | sum(weight > 1) > 1) {
    warning("Weights outside of [0,1] range")
  }
  # multiply weight by adj.M
  hic.table[, adj.M := weight * adj.M]
  if (Plot.MD) {
    MD.plot2(hic.table$adj.M, hic.table$D)
  }
  return(hic.table)
}


# power-law function
.plaw <- function(x, C, alpha) {
  C * x^(-alpha)
}








# # old version
# weight_M <- function(hic.table, quant = 0.75, Plot = TRUE) {
#   # calculate A
#   hic.table[, A := (adj.IF1 + adj.IF2) / 2]
#   # get specified quantile for A
#   A_quant <- quantile(hic.table$A, quant)
#   # get density of A between 1 and A_quant
#   density_A <- density(hic.table$A[hic.table$A > 1 & hic.table$A < A_quant])
#   # get PDF of density_A
#   A_pdf <- approxfun(density_A)
#   # calculate weight - weight is 1 - (pdf(A) / maximum value pdf(A))
#   # 1 - because pdf of A is decreasing
#   # divide by max(density_A$y) to make the weights stay on (0, 1)
#   # set weights for A < peak of A_pdf to 0 and for A > quant_A to 1
#   weight_func <- function(A) {
#     mode_of_A <- density_A$x[which.max(density_A$y)]
#     message("M values corresponding to Average expression < ", mode_of_A, " being weighted to 0.")
#     weight <-ifelse(A < mode_of_A, 0, ifelse(A > A_quant, 1, 
#                                              1 - (A_pdf(A) / max(density_A$y))))
#     return(weight)
#   }
#   weight <- weight_func(hic.table$A)
#   # multiply weight by adj.M
#   hic.table[, adj.M := weight * adj.M]
#   if (Plot) {
#     MD.plot2(hic.table$adj.M, hic.table$D)
#   }
#   return(hic.table)
# }
