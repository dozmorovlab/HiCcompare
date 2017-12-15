#' Detect differences between two jointly normalized Hi-C datasets.
#'
#' @export
#' @param hic.table A hic.table or list of hic.tables output from the
#'     \code{hic_loess} function. hic.table must be jointly normalized
#'     before being entered.
#' @param Plot Logical, should the MD plot showing before/after loess normalization
#'     be output?
#' @param Plot.smooth Logical, defaults to TRUE indicating the MD plot
#'     will be a smooth scatter plot. Set to FALSE for a scatter plot
#'     with discrete points.
#' @param parallel Logical, set to TRUE to utilize the \code{parallel} package's
#'     parallelized computing. Only works on unix operating systems. Only useful if
#'     entering a list of hic.tables.
#' @param BP_param Parameters for BiocParallel. Defaults to bpparam(), see help
#'     for BiocParallel for more information
#'     \url{http://bioconductor.org/packages/release/bioc/vignettes/BiocParallel/
#'     inst/doc/Introduction_To_BiocParallel.pdf}
#'
#'
#' @details  The function takes in a hic.table or a list of hic.table objects created
#'     with the \code{hic_loess} function. If you wish to perform difference
#'     detection on Hi-C data for multiple chromosomes use a list of hic.tables. The process
#'     can be parallelized using the \code{parallel}
#'     setting. The adjusted IF and adjusted M calculated from \code{hic_loess} are used for
#'     difference detection. Fisher's exact test is performed for each pair of IFs between
#'     the two datasets. For example if IF1 = 50 and IF2 = 100 we create a 2x2 table of the
#'     form \code{matrix(c(50, 100, 75, 75), ncol=2)}. The null hypothesis is that IF1 = IF2
#'     = Average expression of IF1 and IF2. We can then perform fisher's exact test on this
#'     2x2 table. FDR p-value adjustment for multiple testing is then performed on a per
#'     distance basis. i.e. at each distance the vector of p-values corresponding to the
#'     interactions occuring at that distance have the FDR multiple testing correction 
#'     applied to them. See methods of Stansfield & Dozmorov 2017 for more details.
#'
#' @return A hic.table with additional columns containing a p-value for the significance
#'     of the difference and the raw fold change between the IFs of the two datasets.
#'
#' @examples
#' # Create hic.table object using included Hi-C data in sparse upper triangular
#' # matrix format
#' data('HMEC.chr22')
#' data('NHEK.chr22')
#' hic.table <- create.hic.table(HMEC.chr22, NHEK.chr22, chr = 'chr22')
#' # Plug hic.table into hic_loess()
#' result <- hic_loess(hic.table, Plot = TRUE)
#' # perform difference detection
#' diff.result <- hic_compare(result, Plot = TRUE)
#'
hic_compare <- function(hic.table,
                     Plot = FALSE, Plot.smooth = TRUE,
                     parallel = FALSE, BP_param = bpparam()) {
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
  # check if single hic.table or list
  if (is.data.table(hic.table)) {
    hic.table <- list(hic.table)
  }
  
  # perform fisher's exact test on each element of the list
  if (parallel) {
    hic.table <- BiocParallel::bplapply(hic.table, .fisher, Plot = Plot) ### May need to change this to calc_z2/calc_z
  } else {
    hic.table <- lapply(hic.table, .fisher, Plot = Plot) ### May need to change this to calc_z2/calc_z
  }
  
  # adjust p-values
  if (parallel) {
    hic.table <- BiocParallel::bplapply(hic.table, .adjust_pval, Plot = Plot) ### May need to change this to calc_z2/calc_z
  } else {
    hic.table <- lapply(hic.table, .adjust_pval, Plot = Plot) ### May need to change this to calc_z2/calc_z
  }
  
  # clean up if single hic.table
  if (length(hic.table) == 1) {
    hic.table <- hic.table[[1]]
  }
  return(hic.table)
}


# background functions for hic_diff

# function to perform fisher's exact test on individual hic.table
.fisher <- function(hic.table, Plot) {
  # make adj.IFs into 2 column matrix
  IF_mat <- cbind(hic.table$adj.IF1, hic.table$adj.IF2) %>% as.matrix()
  # get p-values
  pval <- apply(IF_mat, 1, .get_fisher)
  hic.table[, p.value := pval]
  
  # if (Plot) MD.plot2(hic.table$adj.M, hic.table$D, hic.table$p.adj)
  return(hic.table)
}

# function to take row from IF_mat and produce
# a 2x2 table for input into fisher.test()
# then return the fisher exact p-value
.get_fisher <- function(x) {
  IF_mean <- mean(x) %>% round(., digits = 0)
  x <- round(x, digits = 0)
  m <- matrix(c(x[1], x[2], IF_mean, IF_mean), ncol = 2)
  pval <- fisher.test(m)$p.value
  return(pval)
}


.adjust_pval <- function(hic.table, Plot) {
  # apply distance wise FDR p-value correction
  # split table up for each distance
  temp_list <- S4Vectors::split(hic.table, hic.table$D)
  # combined top 15% of distances into single data.table
  all_dist <- sort(unique(hic.table$D))
  dist_85 <- ceiling(0.85 * length(all_dist))
  temp_list2 <- temp_list[1:dist_85]
  temp_list2[[dist_85+1]] <- data.table::rbindlist(temp_list[(dist_85+1):length(temp_list)])
  temp_list <- temp_list2
  rm("temp_list2")
  # adjust p-values
  temp_list <- lapply(temp_list, function(x) {
    x[, p.adj := p.adjust(p.value, method = 'fdr')]
    return(x)
  })
  # recombine into one table
  hic.table <- rbindlist(temp_list)
  
  if (Plot) MD.plot2(hic.table$adj.M, hic.table$D, hic.table$p.adj)
  return(hic.table)
}


########
# OLD METHODS
########

# # Permutation test function called from hic_loess function or hic_diff
# # function
# .perm.test <- function(data, iterations) {
#   n <- length(data)
#   numerator <- sapply(data, function(x) {
#     # to ignore any NAs in the data
#     if (!is.finite(x))
#       return(NA) else {
#         perm.data <- sample(data, size = iterations, replace = TRUE)
#         test.stat <- ifelse(abs(perm.data) >= abs(x), 1, 0)
#         return(sum(test.stat, na.rm = TRUE))
#       }
#   })
#   p.value <- (numerator + 1)/(iterations + 1)
#   return(p.value)
# }


# # function to calculate a difference threshold based on the
# # distribution of M will produce a difference threshold of 2 * SD(M)
# .calc.diff.thresh <- function(hic.table) {
#   sd_M <- sd(hic.table$adj.M)
#   diff.thresh <- 2 * sd_M
#   return(diff.thresh)
# }
# 
# # Fucntion to calculate p-values based on distance Called from within
# # hic_loess or hic_diff functions uses perm.test function
# .calc.pval <- function(hic.table, diff.thresh = NA, p.adj.method = "fdr",
#                        Plot = TRUE, Plot.smooth = TRUE, iterations = 10000) {
#   # set up vector of included distances
#   all_dist <- sort(unique(hic.table$D))
#   dist_85 <- ceiling(0.85 * length(all_dist))
#   temp <- vector("list", dist_85 + 1)
#   for (dist_idx in seq_len(dist_85)) {
#     temp[[dist_idx]] <- subset(hic.table, D == all_dist[dist_idx])
#     p.temp <- .perm.test(temp[[dist_idx]]$adj.M, iterations = iterations)
#     temp[[dist_idx]][, `:=`(p.value, p.temp)]
#     temp[[dist_idx]][, `:=`(p.adj, p.adjust(p.temp, method = p.adj.method))]
#     # method to check for significant calls when the actual difference
#     # between the two values is very small
#     if (!is.na(diff.thresh)) {
#       # M specifies the log2 fold change between IF1 and IF2. Want to call
#       # differences less than user set diff.thresh fold change not clinically
#       # significant
#       temp[[dist_idx]][, `:=`(p.value, ifelse(p.value < 0.05 & abs(adj.M) <
#                                                 diff.thresh, 0.5, p.value))]
#     }
#   }
#   # for permutation to work need to combine top distances together into
#   # one group
#   temp[[dist_idx + 1]] <- subset(hic.table, D > all_dist[dist_idx])
#   p.temp <- .perm.test(temp[[dist_idx + 1]]$adj.M, iterations = iterations)
#   temp[[dist_idx + 1]][, `:=`(p.value, p.temp)]
#   temp[[dist_idx + 1]][, `:=`(p.adj, p.adjust(p.temp, method = p.adj.method))]
#   # method to check for significant calls when the actual difference
#   # between the two values is very small
#   if (!is.na(diff.thresh)) {
#     temp[[dist_idx + 1]][, `:=`(p.value, ifelse(p.value < 0.05 & abs(adj.M) <
#                                               diff.thresh, 0.5, p.value))]
#   }
#   hic.table <- rbindlist(temp)
#   ## Dont need p.adj so remove it from hic.table before returning
#   hic.table[, `:=`(p.adj, NULL)]
#   # add fold change column
#   hic.table[, `:=`(fold.change, adj.IF2/adj.IF1)]
# 
#   if (Plot) {
#     mdplot <- MD.plot2(M = hic.table$adj.M, D = hic.table$D, p.val = hic.table$p.value,
#                        diff.thresh = diff.thresh, smooth = Plot.smooth)
#     if (!is.null(mdplot)) print(mdplot)
#     return(hic.table)
#   } else {
#     return(hic.table)
#   }
# }
# 
# 
# # Function to perform ranking on hic.table similar to LOLA
# .rank_table <- function(hic.table) {
#   
#   # add average expression to table
#   A <- (hic.table$adj.IF1 + hic.table$adj.IF2) / 2 
#   hic.table[, A := A]
#   
#   ## Rank by distance
#   distance_rank <- data.table::frank(hic.table$D, ties.method = "min")
#   hic.table[, rnkD := distance_rank]
#   
#   # Rank M values at each distance
#   # do i want to rank M at each distance or rank based on all M values for chromosome???
#   ## Rank over all M
#   ranks <- data.table::frank(-abs(hic.table$adj.M), ties.method = "min")
#   hic.table[, rnkM := ranks]
#   
#   # split table up for each distance
#   temp_list <- S4Vectors::split(hic.table, hic.table$D)
#   # combined top 15% of distances into single data.table
#   all_dist <- sort(unique(hic.table$D))
#   dist_85 <- ceiling(0.85 * length(all_dist))
#   temp_list2 <- temp_list[1:dist_85]
#   temp_list2[[dist_85+1]] <- data.table::rbindlist(temp_list[(dist_85+1):length(temp_list)])
#   temp_list <- temp_list2
#   rm("temp_list2")
#   # rank M by distance
#   temp_list <- lapply(temp_list, function(x) {
#     ranks <- data.table::frank(-abs(x$adj.M), ties.method = "min")
#     x[, rnkM_D := ranks]
#     return(x)
#   })
#   # recombine into one table
#   hic.table <- rbindlist(temp_list)
#   
#   # Rank by raw difference
#   diff_rank <- -abs(hic.table$adj.IF1 - hic.table$adj.IF2) %>% data.table::frank(., ties.method = "min")
#   hic.table[, rnkDiff := diff_rank]
#   
#   # Rank by distance
#   # D is equivalent to distance rank
#   # # Rank by p-value
#   # pval_rank <- data.table::frank(hic.table$p.value, ties.method = "min")
#   # hic.table[, rnkPV := pval_rank]
#   
#   # Rank by average expression
#   
#   rank_A <- data.table::frank(-hic.table$A, ties.method = "min")
#   hic.table[, rnkA := rank_A]
#   
#   # Get max rank
#   # max_rank <- hic.table %>% dplyr::select(rnkM, rnkM_D, rnkDiff, rnkPV, rnkD, rnkA) %>% as.matrix() %>% apply(., 1, max)
#   max_rank <- hic.table %>% dplyr::select(rnkM, rnkDiff, rnkA) %>% as.matrix() %>% apply(., 1, max)
#   hic.table[, rnkMax := max_rank]
#   
#   # Get mean rank
#   mean_rank <- hic.table %>% dplyr::select(rnkM, rnkDiff, rnkA) %>% as.matrix() %>% apply(., 1, mean)
#   hic.table[, rnkMean := mean_rank]
#   
#   # sort by max rank
#   hic.table <- hic.table[order(rnkMax, rnkMean),]
#   
#   return(hic.table)
# }
# 
# 
# # version where z scores calculated first then z scores with A < thresh set to 0
# .calc_z <- function(hic.table, quant, Plot) {
#   # add average expression to table
#   A <- (hic.table$adj.IF1 + hic.table$adj.IF2) / 2 
#   hic.table[, A := A]
#   threshold <- quantile((hic.table$A), quant, na.rm = TRUE)
#   Z1 <- (hic.table$adj.M - mean(hic.table$adj.M)) / sd(hic.table$adj.M)
#   # set z-scores where A < threshold to 0
#   Z1[hic.table$A < threshold] <- 0
#   hic.table[, Z := Z1]
#   hic.table[, p.val := 2*pnorm(abs(Z), lower.tail = FALSE)]
#   if (Plot) MD.plot2(hic.table$adj.M, hic.table$D, hic.table$p.val)
#   return(hic.table)
# }
# 
# 
# # version where M values with A < 0 removed before z score calculations
# .calc_z2 <- function(hic.table, quant, Plot) {
#   # add average expression to table
#   A <- (hic.table$adj.IF1 + hic.table$adj.IF2) / 2 
#   hic.table[, A := A]
#   threshold <- quantile((hic.table$A), quant, na.rm = TRUE)
#   new_M <- hic.table$adj.M
#   new_M[hic.table$A < threshold] <- NA
#   Z1 <- (new_M - mean(new_M, na.rm = TRUE)) / sd(new_M, na.rm = TRUE)
#   hic.table[, Z := Z1]
#   hic.table[, p.val := 2*pnorm(abs(Z), lower.tail = FALSE)]
#   if (Plot) MD.plot2(hic.table$adj.M, hic.table$D, hic.table$p.val)
#   return(hic.table)
# }


#### OLD Z SCORE CALCULATIONS
# # Z scores for M, Diff and distance weighting
# .calc_zscores <- function(hic.table) {
#   # calculate z scores
#   Zm <- (hic.table$adj.M - mean(hic.table$adj.M)) / sd(hic.table$adj.M)
#   hic.table[, raw_diff := adj.IF2 - adj.IF1]
#   Zd <- (hic.table$raw_diff - mean(hic.table$raw_diff)) / sd(hic.table$raw_diff)
#   Zmean <- (Zm + Zd) / 2
#   hic.table[, ':=' (Zm = Zm, Zd = Zd, Zmean = Zmean)]
#   # calculate distance weighting
#   dist_weight <- 1 - ((hic.table$D + 1)/max(hic.table$D + 1))
#   hic.table[, D_wt := dist_weight]
#   hic.table[, Zwt := Zmean * D_wt]
#   hic.table[, p.val := pnorm(Zwt)]
#   hic.table[, p.adj := p.adjust(p.val, method = 'fdr')]
#   MD.plot2(hic.table$adj.M, hic.table$D, hic.table$p.adj)
#   MD.plot2(hic.table$adj.M, hic.table$D, hic.table$p.val)
# }
# 
# 
# # Z scores for M, Diff with NO distance weighting
# .calc_zscores2 <- function(hic.table) {
#   # calculate z scores
#   Zm <- (hic.table$adj.M - mean(hic.table$adj.M)) / sd(hic.table$adj.M)
#   hic.table[, raw_diff := adj.IF2 - adj.IF1]
#   Zd <- (hic.table$raw_diff - mean(hic.table$raw_diff)) / sd(hic.table$raw_diff)
#   Zmean <- (Zm + Zd) / 2
#   hic.table[, ':=' (Zm = Zm, Zd = Zd, Zmean = Zmean)]
#   # calculate distance weighting
#   # dist_weight <- 1 - ((hic.table$D + 1)/max(hic.table$D + 1))
#   # hic.table[, D_wt := dist_weight]
#   # hic.table[, Zwt := Zmean * D_wt]
#   hic.table[, p.val := pnorm(Zmean)]
#   hic.table[, p.adj := p.adjust(p.val, method = 'fdr')]
#   MD.plot2(hic.table$adj.M, hic.table$D, hic.table$p.adj)
#   MD.plot2(hic.table$adj.M, hic.table$D, hic.table$p.val)
# }
# 
# 
# # Z scores for M, Diff and distance weighting calculated by Distance
# .calc_zscores3 <- function(hic.table) {
#   hic.table[, raw_diff := adj.IF2 - adj.IF1]
#   # split table up for each distance
#   temp_list <- S4Vectors::split(hic.table, hic.table$D)
#   # combined top 15% of distances into single data.table
#   all_dist <- sort(unique(hic.table$D))
#   dist_85 <- ceiling(0.85 * length(all_dist))
#   temp_list2 <- temp_list[1:dist_85]
#   temp_list2[[dist_85+1]] <- data.table::rbindlist(temp_list[(dist_85+1):length(temp_list)])
#   temp_list <- temp_list2
#   rm("temp_list2")
#   # z score by distance
#   temp_list <- lapply(temp_list, function(x) {
#     Zm <- (hic.table$adj.M - mean(hic.table$adj.M)) / sd(hic.table$adj.M)
#     Zd <- (hic.table$raw_diff - mean(hic.table$raw_diff)) / sd(hic.table$raw_diff)
#     Zmean <- (Zm + Zd) / 2
#     x[, ':=' (Zm = Zm, Zd = Zd, Zmean = Zmean)]
#     return(x)
#   })
#   # recombine into one table
#   hic.table <- rbindlist(temp_list)
#   # calculate distance weighting
#   dist_weight <- 1 - ((hic.table$D + 1)/max(hic.table$D + 1))
#   hic.table[, D_wt := dist_weight]
#   hic.table[, Zwt := Zmean * D_wt]
#   hic.table[, p.val := pnorm(Zwt)]
#   hic.table[, p.adj := p.adjust(p.val, method = 'fdr')]
#   MD.plot2(hic.table$adj.M, hic.table$D, hic.table$p.adj)
#   MD.plot2(hic.table$adj.M, hic.table$D, hic.table$p.val)
# }
