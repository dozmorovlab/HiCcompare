#' Determine the A quantile cutoff to be used
#' 
#' @param hic.table A hic.table object
#' @param A.quantile Logical, Should the quantile of A be tested for a cut off (TRUE)
#'     or should a minimum value of A be tested for a cut off (FALSE).
#' @param SD The standard deviation of the fuzzing used to produce a Hi-C
#'     matrix from your data with few true differences.
#' @param numChanges The number of changes to add into the Hi-C matrix created.
#'    This should be proportional to the resolution of the data. High resolution
#'    data should use more changes i.e. 1MB resolution - 300 changes, 100KB resolution -
#'    1000 changes, etc.
#' @param FC The fold change of the changes added to the Hi-C matrix.
#' @param alpha The alpha level for hypothesis testing.
#' @param Plot logical, should MD plots for the normalization
#'     and difference detection be plotted?
#' 
#' @details This function will take your data and produce an additional
#'     Hi-C matrix using the IF1 vector. Random normal noise will be added
#'     to the vector to create a "fuzzed" matrix with few true differences.
#'     Then the specified number of true changes will be added at the specified
#'     fold change level to the matrices. The HiCcompare procedure is run on the
#'     data and a plot of the MCC based on either the A quantile filtered out
#'     or the A minimum value filtered out will be produced. This is to aid you
#'     in determining what value you should use when analyzing your data with
#'     the hic_compare() function. 
#' @return A plot of the MCC over either the A quantile filtered or the A minimum
#'     value filtered.
#'     
#' @examples 
#' data('HMEC.chr22')
#' data('NHEK.chr22')
#' hic.table <- create.hic.table(HMEC.chr22, NHEK.chr22, chr = 'chr22')
#' filter_params(hic.table)


filter_params <- function(hic.table, A.quantile = TRUE, filter_n = 100,
                          SD = 2, numChanges = 300, FC = 3, alpha = 0.05, Plot = FALSE) {
  # create new fuzzed table
  new.table <- randomize_IFs(hic.table, SD)
  # delete huge differences
  new.table <- new.table[abs(M) < 2, ]
  
  # add true changes to fuzzed table
  sample_space <- 1:nrow(new.table)
  changes <- sample(sample_space, numChanges)
  # set IFs to mean IF then multiply one by FC
  meanIF <- ((new.table[changes,]$IF1 + new.table[changes,]$IF2) / 2) %>% round() %>% as.integer()
  suppressWarnings(new.table[changes, IF1 := meanIF ])
  suppressWarnings(new.table[changes, IF2 := meanIF])
  midpoint <- floor(numChanges/2)
  newIF1 <- new.table[changes[1:midpoint],]$IF1 * FC %>% as.integer()
  newIF2 <- new.table[changes[(midpoint+1):numChanges],]$IF2 * FC %>% as.integer()
  new.table[changes[1:midpoint], IF1 :=  newIF1]
  new.table[changes[(midpoint+1):numChanges], IF2 :=  newIF2]
  new.table = new.table[, M := log2(IF2/IF1)]
  truth <- rep(0, nrow(new.table))
  truth[changes] <- 1
  new.table[, truth := truth]
  
  # normalize & detect differences
  new.table <- hic_loess(new.table, Plot = Plot)
  new.table <- hic_compare(new.table, Plot = Plot)
  
  # calculate MCC
  TP <- vector(length = 50)
  FP <- vector(length = 50)
  FN <- vector(length = 50)
  TN <- vector(length = 50)
  if (A.quantile) {
    A_seq <- seq(1, 100, by = 1)
    for (i in seq_along(A_seq)) {
      tmp.table <- hic_compare(new.table, ngroups = A_seq[i], filter_n = filter_n, adjust.dist = TRUE, p.method = 'fdr', Plot = FALSE)
      TP[i] <- sum(tmp.table$p.adj < alpha & tmp.table$truth == 1)
      FP[i] <- sum(tmp.table$p.adj < alpha & tmp.table$truth == 0)
      FN[i] <- sum(tmp.table$p.adj >= alpha & tmp.table$truth == 1)
      TN[i] <- sum(tmp.table$p.adj >= alpha & tmp.table$truth == 0)
    }
  } else {
    A_seq <- seq(1, 50, by = 1)
    for (i in seq_along(A_seq)) {
      tmp.table <- hic_compare(new.table, A.min = A_seq[i], adjust.dist = TRUE, p.method = 'fdr', Plot = FALSE)
      TP[i] <- sum(tmp.table$p.adj < alpha & tmp.table$truth == 1)
      FP[i] <- sum(tmp.table$p.adj < alpha & tmp.table$truth == 0)
      FN[i] <- sum(tmp.table$p.adj >= alpha & tmp.table$truth == 1)
      TN[i] <- sum(tmp.table$p.adj >= alpha & tmp.table$truth == 0)
    }
  }
  # Calculate MCC
  MCC <- ((TP * TN) - (FP * FN)) / 
    (sqrt((TP + FP)) * sqrt((TP + FN)) * sqrt((TN + FP)) *
       sqrt((TN + FN)))
  ## calculate FN & FP rates
  FPR <- FP / (FP + TP)
  FNR <- FN / (FN + TN)
  TPR <- TP / (TP + FP)
  # PLOT MCC
  if (A.quantile) {
    # plot(MCC ~ A_seq, type = 'l', col = 'red', main = 'MCC by A quantile filtered', ylab = 'MCC', xlab = 'A Quantile filtered', ylim = c(0,1)) # old method
    plot(MCC ~ A_seq, type = 'l', col = 'red', main = 'MCC by A groups filtered', ylab = 'MCC', xlab = 'A Quantile filtered', ylim = c(0,1)) # new thing
    lines(FPR ~ A_seq, col = 'blue')
    lines(FNR ~ A_seq, col= 'green')
    lines(TPR ~ A_seq, col = 'black')
    legend('topright', legend = c('MCC', 'FPR', 'FNR', 'TPR'), fill = c('red', 'blue', 'green', 'black'))
  } else {
    plot(MCC ~ A_seq, type = 'l', col = 'red', main = 'MCC by A quantile filtered', ylab = 'MCC', xlab = 'A minimum filtered', ylim = c(0,1))
    lines(FPR ~ A_seq, col = 'blue')
    lines(FNR ~ A_seq, col= 'green')
    lines(TPR ~ A_seq, col = 'black')
    legend('topright', legend = c('MCC', 'FPR', 'FNR', 'TPR'), fill = c('red', 'blue', 'green', 'black'))
  }
  
}


# old version
# filter_params <- function(hic.table, A.quantile = TRUE, SD = 2, numChanges = 300, FC = 3, alpha = 0.05, Plot = FALSE) {
#   # create new fuzzed table
#   new.table <- randomize_IFs(hic.table, SD)
#   # delete huge differences
#   new.table <- new.table[abs(M) < 2, ]
#   
#   # add true changes to fuzzed table
#   sample_space <- 1:nrow(new.table)
#   changes <- sample(sample_space, numChanges)
#   # set IFs to mean IF then multiply one by FC
#   meanIF <- ((new.table[changes,]$IF1 + new.table[changes,]$IF2) / 2) %>% round() %>% as.integer()
#   suppressWarnings(new.table[changes, IF1 := meanIF ])
#   suppressWarnings(new.table[changes, IF2 := meanIF])
#   midpoint <- floor(numChanges/2)
#   newIF1 <- new.table[changes[1:midpoint],]$IF1 * FC %>% as.integer()
#   newIF2 <- new.table[changes[(midpoint+1):numChanges],]$IF2 * FC %>% as.integer()
#   new.table[changes[1:midpoint], IF1 :=  newIF1]
#   new.table[changes[(midpoint+1):numChanges], IF2 :=  newIF2]
#   new.table = new.table[, M := log2(IF2/IF1)]
#   truth <- rep(0, nrow(new.table))
#   truth[changes] <- 1
#   new.table[, truth := truth]
#   
#   # normalize & detect differences
#   new.table <- hic_loess(new.table, Plot = Plot)
#   new.table <- hic_compare(new.table, Plot = Plot)
#   
#   # calculate MCC
#   TP <- vector(length = 50)
#   FP <- vector(length = 50)
#   FN <- vector(length = 50)
#   TN <- vector(length = 50)
#   if (A.quantile) {
#     A_seq <- seq(0.01, 0.5, by = 0.01)
#     for (i in seq_along(A_seq)) {
#       tmp.table <- hic_compare(new.table, A.quantile = A_seq[i], adjust.dist = TRUE, p.method = 'fdr', Plot = FALSE)
#       TP[i] <- sum(tmp.table$p.adj < alpha & tmp.table$truth == 1)
#       FP[i] <- sum(tmp.table$p.adj < alpha & tmp.table$truth == 0)
#       FN[i] <- sum(tmp.table$p.adj >= alpha & tmp.table$truth == 1)
#       TN[i] <- sum(tmp.table$p.adj >= alpha & tmp.table$truth == 0)
#     }
#   } else {
#     A_seq <- seq(1, 50, by = 1)
#     for (i in seq_along(A_seq)) {
#       tmp.table <- hic_compare(new.table, A.min = A_seq[i], adjust.dist = TRUE, p.method = 'fdr', Plot = FALSE)
#       TP[i] <- sum(tmp.table$p.adj < alpha & tmp.table$truth == 1)
#       FP[i] <- sum(tmp.table$p.adj < alpha & tmp.table$truth == 0)
#       FN[i] <- sum(tmp.table$p.adj >= alpha & tmp.table$truth == 1)
#       TN[i] <- sum(tmp.table$p.adj >= alpha & tmp.table$truth == 0)
#     }
#   }
#   # Calculate MCC
#   MCC <- ((TP * TN) - (FP * FN)) / 
#     (sqrt((TP + FP)) * sqrt((TP + FN)) * sqrt((TN + FP)) *
#        sqrt((TN + FN)))
#   ## calculate FN & FP rates
#   FPR <- FP / (FP + TP)
#   FNR <- FN / (FN + TN)
#   TPR <- TP / (TP + FP)
#   # PLOT MCC
#   if (A.quantile) {
#     plot(MCC ~ A_seq, type = 'l', col = 'red', main = 'MCC by A quantile filtered', ylab = 'MCC', xlab = 'A Quantile filtered', ylim = c(0,1))
#     lines(FPR ~ A_seq, col = 'blue')
#     lines(FNR ~ A_seq, col= 'green')
#     lines(TPR ~ A_seq, col = 'black')
#     legend('topright', legend = c('MCC', 'FPR', 'FNR', 'TPR'), fill = c('red', 'blue', 'green', 'black'))
#   } else {
#     plot(MCC ~ A_seq, type = 'l', col = 'red', main = 'MCC by A quantile filtered', ylab = 'MCC', xlab = 'A minimum filtered', ylim = c(0,1))
#     lines(FPR ~ A_seq, col = 'blue')
#     lines(FNR ~ A_seq, col= 'green')
#     lines(TPR ~ A_seq, col = 'black')
#     legend('topright', legend = c('MCC', 'FPR', 'FNR', 'TPR'), fill = c('red', 'blue', 'green', 'black'))
#   }
#   
# }


# function to add noise to IFs of one matrix 
randomize_IFs <- function(hic.table, SD) {
  # copy first IF vector
  newIF2 <- hic.table$IF1
  # add constant offset
  newIF2 <- newIF2 + 5
  # add random noise
  newIF2 <- newIF2 + rnorm(length(newIF2), 0, SD)
  # check for 0's and negatives
  newIF2[newIF2 <= 0] <- 1
  # create new hic.table with new IF vectors
  sparse1 <- cbind(hic.table$start1, hic.table$start2, hic.table$IF1)
  sparse2 <- cbind(hic.table$start1, hic.table$start2, newIF2)
  temp.table <- create.hic.table(sparse1, sparse2, chr = hic.table$chr1[1])
  return(temp.table)
}
