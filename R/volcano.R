volcano <- function(hic.table) {
  hic.table[, A := (adj.IF1 + adj.IF2) / 2]
  max_D <- max(hic.table$D)
  max_A <- max(log10(hic.table$A))
  hic.table[, volcano := -log10(p.value) * (log10(A) / log10(max_A)) * (1 - (D / max_D))]
  plot(hic.table$adj.M, hic.table$volcano)
  return(hic.table)
}