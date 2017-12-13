# function to adjust p-values based on distance holm method
adjust_pval <- function(hic.table, theta, alpha = 0.05) {
  # get distance percentages
  d_percent <- ((hic.table$D + 1) / max(hic.table$D + 1)) * 100
  # calculate threshold for checking if p-value is < thresh to be a rejection of the null hypothesis
  threshhold <- alpha / (d_percent * theta)
  P <- hic.table$p.val
  P[P > threshhold] <- 1
  hic.table[, thresh := threshhold]
  hic.table[, p.adj := P]
  return(hic.table)
} 




