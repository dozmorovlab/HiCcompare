weight_M <- function(hic.table, quant = 0.75, Plot = TRUE) {
  # calculate A
  hic.table[, A := (adj.IF1 + adj.IF2) / 2]
  # get specified quantile for A
  A_quant <- quantile(hic.table$A, quant)
  # get density of A between 1 and A_quant
  density_A <- density(hic.table$A[hic.table$A > 1 & hic.table$A < A_quant])
  # get PDF of density_A
  A_pdf <- approxfun(density_A)
  # calculate weight - weight is 1 - (pdf(A) / maximum value pdf(A))
  # 1 - because pdf of A is decreasing
  # divide by max(density_A$y) to make the weights stay on (0, 1)
  # set weights for A < peak of A_pdf to 0 and for A > quant_A to 1
  weight_func <- function(A) {
    mode_of_A <- density_A$x[which.max(density_A$y)]
    message("M values corresponding to Average expression < ", mode_of_A, " being weighted to 0.")
    weight <-ifelse(A < mode_of_A, 0, ifelse(A > A_quant, 1, 
                                             1 - (A_pdf(A) / max(density_A$y))))
    return(weight)
  }
  weight <- weight_func(hic.table$A)
  # multiply weight by adj.M
  hic.table[, adj.M := weight * adj.M]
  if (Plot) {
    MD.plot2(hic.table$adj.M, hic.table$D)
  }
  return(hic.table)
}
