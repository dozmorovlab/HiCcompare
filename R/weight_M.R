weight_M <- function(hic.table, A.min = 1, quant = 0.75, Plot = TRUE) {
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
  A_seq <- seq(A.min, max(hic.table$A[hic.table$A > A.min & hic.table$A < A_quant]), by = 0.001)
  hist(hic.table$A[hic.table$A > A.min & hic.table$A < A_quant], freq = FALSE, main = 'Powerlaw fit to A', xlab = 'A')
  lines(.plaw(A_seq, est$par[1], est$par[2]) ~ A_seq)
  
  
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
  if (Plot) {
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
