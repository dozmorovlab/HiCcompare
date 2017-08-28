# function to adjust p-values based on weighting for distance

p_adjust <- function(pval, D) {
  zscore <- qnorm(pval)
  scaled_D <- (D + 1) / max(D + 1)
  zscore <- zscore / scaled_D
  new_pval <- pnorm(zscore)
  new_pval <- p.adjust(new_pval, method = 'bonferroni')
  return(new_pval)
}

p_adjust <- function(pval, D) {
  maxD <- max(D + 1)
  normD <- (D + 1) / maxD
  zscore <- qnorm(pval)
}


#
# tab <- hic_tables[[1]]
#
# adj.p <- p_adjust(tab$p.value, tab$D)
#
# tab[, adj.p := adj.p]
# sum(adj.p < 0.05)
# sum(pval < 0.05)
# MD.plot2(tab$M, tab$D, tab$adj.p)



# calc_p_AC <- function(x, y) {
#   N1 <- sum(x, na.rm = TRUE)
#   N2 <- sum(y, na.rm = TRUE)
#   message('N1: ', N1, ' N2: ', N2)
#   # p.AC <- ifelse(x + y < 100,
#   #                (N2/N1)^y * ( choose(round(x + y),round(x)) * (1 / (1 + (N2/N1))^(x + y + 1))),
#   #                exp(y * (log(N2) - log(N1)) + .factorial.approx(x + y) - .factorial.approx(x) - .factorial.approx(y) - ((x + y + 1) * log(1 + (N2 / N1))) )
#   # )
#
#   p.AC <- exp((y * (log(N2) - log(N1))) + (.factorial.approx(x + y) - (.factorial.approx(x) + .factorial.approx(y) + (x + y + 1) * log(1 + (N2 / N1)))))
#   #p.AC <- (N2/N1)^y * ( choose(round(x + y),round(x)) * (1 / (1 + (N2/N1))^(x + y + 1)))
#   return(p.AC)
# }
#
#
# # stirling's approximation for factorials
# # returns log(x!)
# .factorial.approx = function(x) {
#   result = ifelse(x == 0, 0, x*log(x)-x+0.5*log(2*pi*x))
#   return(result)
# }
