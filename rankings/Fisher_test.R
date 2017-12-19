# We compare two samples. For a given point, IFs may be different, or similar. We need the significance of the differences
# IF1 = 5, IF2 = 15. IF1 + IF2 = 20 - low total IF. Null hypothesis - IF1 = IF2 = 10
m1 <- matrix(c(5, 15, 10, 10), ncol = 2)
fisher.test(m1)$p.value
# IF1 = 50, IF2 = 150. IF1 + IF2 = 200 - high total IF. Null hypothesis - IF1 = IF2 = 100
m2 <- matrix(c(50, 150, 100, 100), ncol = 2)
fisher.test(m2)$p.value

# Test of proportions, as in Kal, A. J., A. J. van Zonneveld, V. Benes, M. van den Berg, M. G. Koerkamp, K. Albermann, N. Strack, et al. “Dynamics of Gene Expression Revealed by Comparison of Serial Analysis of Gene Expression Transcript Profiles from Yeast Grown on Two Different Carbon Sources.” Molecular Biology of the Cell 10, no. 6 (June 1999): 1859–72.
# https://www.ncbi.nlm.nih.gov/pubmed/10359602


# After normalization, we assume the samples have the same number of reads
N1 = 1000 # total number of reads in the first sample
N2 = 1000 # total number of reads in the second sample

kal_test <- function(IF1, IF2, N1, N2) {
  p1 = IF1 / N1 # proportion for IF in the first sample
  p2 = IF2 / N2 # proportion for IF in the second sample
  p0 = (IF1 + IF2) / (N1 + N2)
  # Z score calculation, as per A. J. Kal et al., “Dynamics of Gene Expression Revealed by Comparison of Serial Analysis of Gene Expression Transcript Profiles from Yeast Grown on Two Different Carbon Sources,” Molecular Biology of the Cell 10, no. 6 (June 1999): 1859–72.
  Z = (p1 - p2) / (sqrt( (p0 * (1 - p0) / N1) + (p0 * (1 - p0) / N2) ))
  return(2 * pnorm(-abs(Z))) # P-value
}

kal_test(IF1 = 5, IF2 = 15, N1 = N1, N2 = N2)
kal_test(IF1 = 50, IF2 = 150, N1 = N1, N2 = N2)
