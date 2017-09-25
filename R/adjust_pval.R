# Adjust p-values based on Average IF between matrices


adjust_pval <- function(hic.table) {
  # calculate A
  hic.table[, A := (adj.IF1 + adj.IF2) / 2]
  # get powerlaw fit for A
  fit <- .fit_powerlaw(hic.table$A)
  message("Powerlaw Xmin = ", fit$xmin, "; Powerlaw alpha = ", fit$alpha)
  # calculate adjusted p-value 
  hic.table[, adj.p := p.value * (1 + .plaw(A, fit$xmin, fit$alpha))]
  # for any A < 1 set adj.p to 1
  hic.table$adj.p[hic.table$A < 1] <- 1
  # set p-values > 1 to 1
  hic.table$adj.p[hic.table$adj.p > 1] <- 1
  MD.plot2(hic.table$adj.M, hic.table$D, hic.table$adj.p)
  return(hic.table)
}






# power-law function
.plaw <- function(x, C, alpha) {
  C * x^(-alpha)
}


# fit power law by deleting most extreme values
.fit_powerlaw = function(x) {
  x = na.omit(x) # remove any NAs
  fit = igraph::fit_power_law(x) # get initial fit
  pval = fit$KS.p
  # loop until fit is good
  while(pval <= 0.05) {
    x = x[-which(x == max(x))]
    fit = igraph::fit_power_law(x)
    pval = fit$KS.p
  }
  return(fit)
}
