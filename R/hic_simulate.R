# to generate NB component of the cell value
.sim.NBvalue <- function(Bmean, r = 10) {
  value <- rnbinom(1, mu = Bmean, size = r)
  return(value)
}


# power-law function
.powerlaw <- function(distance, C, alpha) {
  ifelse(distance == 0, C, C * distance^(-alpha))
}

# function for generating P(IF = 0) based on linear equation
.prop.zero.linear <- function(distance, prop.zero.slope) {
  prop.zero.slope * distance
}

# simulate matrix function will create a two full contact maps.
# Matrices will have same signal component but different noise
# components
.sim.mat <- function(nrow = 100, ncol = 100, medianIF, sdIF, powerlaw.alpha,
                     sd.alpha, prop.zero.slope) {
  # check for invalid proportion of zeros slope
  if (.prop.zero.linear(prop.zero.slope, nrow - 1) > 1) {
    stop("prop.zero.slope is too large and will produce probabilities
         greater than 1 at the maximum distance in the matrix.")
  }
  cell1 <- matrix(nrow = nrow, ncol = ncol)
  cell2 <- matrix(nrow = nrow, ncol = ncol)
  col_num <- 1

  j <- unlist(Map(seq, 1:ncol(cell1), MoreArgs=list(ncol(cell1))))
  i <- rep(seq_len(nrow(cell1)), times=head(rev(seq_len(ncol(cell1))), nrow(cell1)))
  distance <- j - i + 1
  Bmean <- .powerlaw(distance, medianIF, powerlaw.alpha)
  noise.sd <- .powerlaw(distance, sdIF, sd.alpha)
  idx <- cbind(i, j)
  cell1[idx] <- round(Bmean) + round(rnorm(i, 0, noise.sd))
  cell2[idx] <- round(Bmean) + round(rnorm(i, 0, noise.sd))
  prob.zero <- .prop.zero.linear(distance, prop.zero.slope)
  u1 <- runif(length(prob.zero))
  u2 <- runif(length(prob.zero))
  zero_index1 <- ifelse(u1 < prob.zero, TRUE, FALSE)
  zero_index2 <- ifelse(u2 < prob.zero, TRUE, FALSE)
  cell1[idx[zero_index1,]] <- 0
  cell2[idx[zero_index2,]] <- 0
  idx2 <- idx[, c(2,1)]
  cell1[idx2] <- cell1[idx]
  cell2[idx2] <- cell2[idx]
  # for (i in 1:dim(cell1)[1]) {
  #   for (j in col_num:dim(cell1)[2]) {
  #     distance <- j - i + 1
  #     Bmean <- .powerlaw(distance, medianIF, powerlaw.alpha)
  #     noise.sd <- .powerlaw(distance, sdIF, sd.alpha)
  #     cell1[i, j] <- round(Bmean) + round(rnorm(1, 0, noise.sd))
  #     cell1[j, i] <- cell1[i, j]
  #     cell2[i, j] <- round(Bmean) + round(rnorm(1, 0, noise.sd))
  #     cell2[j, i] <- cell2[i, j]
  #     # add in proportion of 0's
  #     prob.zero <- .prop.zero.linear(distance, prop.zero.slope)
  #     u <- runif(2)
  #     if (u[1] < prob.zero) {
  #       cell1[i, j] <- 0
  #       cell1[j, i] <- 0
  #     }
  #     if (u[2] < prob.zero) {
  #       cell2[i, j] <- 0
  #       cell2[j, i] <- 0
  #     }
  #   }
  #   col_num <- col_num + 1
  # }
  # check for negative values
  cell1[cell1 < 0] <- 0
  cell2[cell2 < 0] <- 0
  return(list(cell1, cell2))
}

# default bias functions

.normal.bias <- function(distance) {
  (1 + exp(-((distance - 20)^2)/(2 * 30))) * 4
}


.no.bias <- function(distance) {
  1
}

# function to add bias to a matrix
.add.bias <- function(mat, slope, medianIF, powerlaw.alpha,
                      biasFunc = .normal.bias) {
  col_num <- 1
  for (i in 1:dim(mat)[1]) {
    for (j in col_num:dim(mat)[2]) {
      distance <- j - i + 1
      Bmean <- .powerlaw(distance, medianIF, powerlaw.alpha)
      new.value <- round(mat[i, j] * biasFunc(distance))
      new.value <- ifelse(new.value < 0, 0, new.value)
      mat[i, j] <- new.value
      mat[j, i] <- mat[i, j]
    }
    col_num <- col_num + 1
  }
  return(mat)
}


# add true differences to the matrices
.sim.differences <- function(mat1, mat2, fold.change = 2,
                             i.range, j.range) {
  for (n in 1:length(i.range)) {
    i <- i.range[n]
    j <- j.range[n]
    # get which direction the difference is
    diff_direction <- sign(mat1[i, j] - mat2[i, j])
    newIF <- as.integer(round(mat1[i, j] * fold.change^diff_direction))
    mat1[i, j] <- ifelse(newIF < 1, 1, newIF)
    mat1[j, i] <- mat1[i, j]
  }
  return(mat1)
}


# wrapper function for simulation studies will generate the matrices
# and perform HiCdiff analysis on them

#' Simulate a Hi-C matrix and perform hic_diff analysis on it
#'
#' @export
#' @param nrow Number of rows and columns of the full matrix
#' @param medianIF The starting value for a power law distribution
#'     for the interaction frequency of the matrix. Should use the median
#'     value of the IF at distance = 0. Typical values for 1MB data
#'     are around 50,000.
#'     For 500kb data typical values are 25,000. For 100kb data, 4,000.
#'     For 50kb data, 1,800.
#' @param sdIF The estimated starting value for a power law distriubtion
#'     for the standard deviaton of the IFs. Should use the SD of the IF at
#'     distance = 0. Typical value for 1MB data is 19,000.
#' @param powerlaw.alpha The exponential parameter for the power law
#'     distribution for the median IF. Typical values are 1.6 to 2.
#'     Defaults to 1.8.
#' @param sd.alpha The exponential parameter for the power law
#'     distribution for the SD of the IF. Typical values are 1.8 to 2.2.
#'     Defaults to 1.9.
#' @param prop.zero.slope The slope to be used for a linear function of
#'     the probability of zero in matrix = slope * distance
#' @param biasFunc A function used for adding bias to one of the simulated
#'     matrices. Should take an input of unit distance and generally have
#'     the form of 1 + Probability Density Function with unit distance as the
#'     random variable. Can also use a constant as a scaling factor
#'     to add a global offset to one of the matrices. The output of the bias
#'     function will be multiplied to the IFs of one matrix.
#'     Included are a normal kernel bias and a no bias function. If no function
#'     is entered, a normal kernel bias with an additional global scaling
#'     factor of 4 will be used. To use no bias set biasFunc = .no.bias, see
#'     examples section.
#' @param fold.change The fold change you want to introduce for true differences
#'     in the simulated matrices. Defaults to NA for no fold change added.
#' @param i.range The row numbers for the cells that you want to introduce true
#'     differences at. Must be same length as j.range.
#'     Defaults to NA for no changes added.
#' @param j.range The column numbers for the cells that you want to introduce
#'     true differences at.  Must be same length as
#'     Defaults to NA for no changes added.
#' @param Plot Logical, should the HiCdiff plots be output? Defaults to TRUE.
#' @param scale Logical, Should scaling be applied for the HiCdiff procedure?
#'     Defaults to TRUE.
#' @param alpha Type I error rate parameter. At what level should a significant
#'     difference be defined. Defaults to 0.05.
#' @param diff.thresh Parameter for hic_diff procedure. See ?hic_diff for more
#'     help. Defaults to NA.
#'
#' @return A list containing the true positive rate (TPR), the specificity (SPC),
#'     the p-values, the hic.table object, true differences - a data.table
#'     of the rows of the hic.table where a true difference was applied, the truth
#'     vector - a vector of 0's and 1's where 1 indicates a true
#'     difference was applied to that cell, sim.table - the hic.table object for
#'     the simulate matrices before hic_loess and hic_diff was run on it.
#'
#' @examples
#' # simulate two matrices with no fold changes introduced using default values
#' sim <- hic_simulate()
#'
#' # example of bias functions
#' ## the default function used
#' .normal.bias = function(distance) {
#'   (1 + exp(-((distance - 20)^2) / (2*30))) * 4
#' }
#'
#' ## an additional bias function
#' .no.bias = function(distance) {
#'   1
#' }
#'
#' # simulate matrices with 200 true differences using no bias
#' i.range = sample(1:100, replace=TRUE)
#' j.range = sample(1:100, replace=TRUE)
#' sim2 <- hic_simulate(nrow=100, biasFunc = .no.bias, fold.change = 5,
#'                      i.range = i.range, j.range = j.range)
#'
#'
hic_simulate <- function(nrow = 100, medianIF = 50000, sdIF = 14000,
                         powerlaw.alpha = 1.8,
                         sd.alpha = 1.9, prop.zero.slope = 0.001,
                         biasFunc = .normal.bias, fold.change = NA,
                         i.range = NA, j.range = NA, Plot = TRUE,
                         scale = TRUE, alpha = 0.05,
                         diff.thresh = NA) {

  if (is.na(fold.change) & (is.na(i.range[1]) | is.na(j.range[1]))) {
    i.range <- 1
    j.range <- 1
  }
  if (!is.na(fold.change) & (is.na(i.range[1]) | is.na(j.range[1]))) {
    stop("Error: Please enter values for i.range and j.range if
         you wish to produce a fold change in the simulated matrix")
  }
  ncol <- nrow
  # simulate matrices
  sims <- .sim.mat(nrow, ncol, medianIF, sdIF, powerlaw.alpha, sd.alpha,
                   prop.zero.slope)

  # if fold.change = NA no true differences will be added to the matrix
  if (!is.na(fold.change)) {
    # make sure no cells are duplicated
    temp.tab <- data.frame(i = i.range, j = j.range)
    temp.tab <- unique(temp.tab)
    i.range <- temp.tab$i
    j.range <- temp.tab$j
    diff <- .sim.differences(sims[[2]], sims[[1]], fold.change, i.range,
                             j.range)
    sims[[2]] <- diff
  }
  # add in sample specific bias to one matrix
  sims[[1]] <- .add.bias(sims[[1]], bias.slope, medianIF, powerlaw.alpha,
                         biasFunc = biasFunc)
  # perform HiCloess on simulated data convert matrix to sparse format
  colnames(sims[[1]]) <- 1:nrow
  colnames(sims[[2]]) <- 1:nrow
  sims[[1]] <- full2sparse(sims[[1]])
  sims[[2]] <- full2sparse(sims[[2]])
  backup.sim.table <- create.hic.table(sims[[1]], sims[[2]], chr = "ChrSim",
                                       scale = FALSE)
  sims <- create.hic.table(sims[[1]], sims[[2]], chr = "ChrSim",
                           scale = scale)
  normed <- hic_loess(sims, Plot = Plot, diff.thresh = diff.thresh,
                      check.differences = TRUE)
  pvals <- normed$p.value

  # get the true differences
  true.diffs <- data.table(i = c(i.range, j.range), j = c(j.range, i.range))
  true.diffs <- left_join(true.diffs, normed, by = c(i = "start1", j = "start2"))
  # remove duplicated rows
  true.diffs <- true.diffs[!duplicated(true.diffs), ]
  true.diffs <- as.data.table(na.omit(true.diffs))

  # calculate sensitivty and specificity
  true.pos <- sum(true.diffs$p.value < alpha)
  false.pos <- sum(normed$p.value < alpha, na.rm = TRUE) - true.pos
  false.neg <- nrow(true.diffs) - true.pos
  true.neg <- nrow(normed) - true.pos - false.pos - false.neg
  TPR <- true.pos/(true.pos + false.neg)
  SPC <- true.neg/(true.neg + false.pos)

  # make vectors to feed into ROC packages vector of p-value decision

  temp.true.diffs <- true.diffs[, `:=`(truth, 1)]
  temp.true.diffs <- temp.true.diffs[, c("i", "j", "truth"), with = FALSE]
  truth <- left_join(normed, temp.true.diffs, by = c(start1 = "i", start2 = "j"))
  truth$truth[is.na(truth$truth)] <- 0
  truth <- truth$truth

  results <- list(TPR = TPR, SPC = SPC, pvals = normed$p.value, hic.table = normed,
                  true.diff = true.diffs, truth = truth, sim.table = backup.sim.table)
  if (!is.na(fold.change)) {
    message("True Positives: ", true.pos, " Total added differences: ",
                nrow(true.diffs), " True Negatives: ", true.neg, sep = "")
    message("TPR: ", TPR, sep = "")
    message("SPC: ", SPC, sep = "")
  }
  return(results)
}


# function to get TPR/SPC and generate truth matrix results of HiCdiff
# on non-hicloess normalized data enter sim.table as a hic.table for
# the pre-normalized simulated matrices. Normalization can be done
# using some other method than hic_loess

#' Compare other normalization methods on simulated data
#'
#' @export
#' @param sim.table the sim.table object output from hic_simulate
#' @param i.range The row numbers for the cells that you want to introduce
#'     true differences at. Must be same length as j.range.
#' @param j.range The column numbers for the cells that you want to introduce
#'     true differences at.  Must be same length as i.range.
#' @param Plot Logical, should the HiCdiff plots be output? Defaults to TRUE.
#' @param alpha Type I error rate parameter. At what level should a significant
#'     difference be defined. Defaults to 0.05.
#' @param diff.thresh Parameter for hic_diff procedure. see ?hic_diff for
#'     more details.
#'
#' @return A list containing the true positive rate (TPR), the specificity (SPC),
#'     the p-values, the hic.table object, true differences - a data.table
#'     of the rows of the hic.table where a true difference was applied, the truth
#'     vector - a vector of 0's and 1's where 1 indicates a true
#'     difference was applied to that cell.
#' @examples
#' i.range = sample(1:100, replace=TRUE)
#' j.range = sample(1:100, replace=TRUE)
#' sim <- hic_simulate(i.range = i.range, j.range = j.range, fold.change = 2)
#' mat1 <- sim$sim.table[, c('start1', 'start2', 'IF1'), with=FALSE]
#' mat2 <- sim$sim.table[, c('start1', 'start2', 'IF2'), with=FALSE]
#' mat1 <- sparse2full(mat1) %>% KRnorm
#' mat2 <- sparse2full(mat2) %>% KRnorm
#' colnames(mat1) <- 1:ncol(mat1)
#' colnames(mat2) <-1:ncol(mat2)
#' mat1 <- full2sparse(mat1)
#' mat2 <- full2sparse(mat2)
#' new.tab <- create.hic.table(mat1, mat2, chr= 'chrsim')
#' sim2 <- sim.other.methods(new.tab, i.range = i.range , j.range = j.range)
#'
sim.other.methods <- function(sim.table, i.range, j.range, Plot = TRUE,
                              alpha = 0.05, diff.thresh = NA) {
  if (ncol(sim.table) < 11) {
    sim.table <- sim.table[, `:=`(adj.IF1, IF1)]
    sim.table <- sim.table[, `:=`(adj.IF2, IF2)]
    sim.table <- sim.table[, `:=`(adj.M, M)]
  }
  diffs <- hic_diff(sim.table, Plot = Plot, diff.thresh = diff.thresh)
  pvals <- diffs$p.value

  true.diffs <- data.table(i = c(i.range, j.range), j = c(j.range, i.range))
  true.diffs <- left_join(true.diffs, diffs, by = c(i = "start1", j = "start2"))
  # remove duplicated rows
  true.diffs <- true.diffs[!duplicated(true.diffs), ]
  true.diffs <- as.data.table(na.omit(true.diffs))

  # calculate sensitivty and specificity
  true.pos <- sum(true.diffs$p.value < alpha)
  false.pos <- sum(diffs$p.value < alpha, na.rm = TRUE) - true.pos  #### changed to < 0.1
  false.neg <- nrow(true.diffs) - true.pos
  true.neg <- nrow(diffs) - true.pos - false.pos - false.neg
  TPR <- true.pos/(true.pos + false.neg)
  SPC <- true.neg/(true.neg + false.pos)

  # make vectors to feed into ROC packages vector of p-value decision
  temp.true.diffs <- true.diffs[, `:=`(truth, 1)]
  temp.true.diffs <- temp.true.diffs[, c("i", "j", "truth"), with = FALSE]
  truth <- left_join(diffs, temp.true.diffs, by = c(start1 = "i", start2 = "j"))
  truth$truth[is.na(truth$truth)] <- 0
  truth <- truth$truth


  results <- list(TPR = TPR, SPC = SPC, pvals = diffs$p.value, hic.table = diffs,
                  true.diff = true.diffs, truth = truth)
  message("True Positives: ", true.pos, " Total added differences: ",
              nrow(true.diffs), " True Negatives: ", true.neg, sep = "")
  message("TPR: ", TPR, sep = "")
  message("SPC: ", SPC, sep = "")
  return(results)
}
