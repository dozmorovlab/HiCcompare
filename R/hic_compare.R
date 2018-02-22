#' Detect differences between two jointly normalized Hi-C datasets.
#'
#' @export
#' @param hic.table A hic.table or list of hic.tables output from the
#'     \code{hic_loess} function. hic.table must be jointly normalized
#'     before being entered.
#' @param A.min The required value of A in order for a differences to be
#'     considered. All Z-scores
#'     where the corresponding A value is < A.min will be set to 0. 
#'     Defaults to NA. If NA, then the 10th percentile of A will automatically
#'     be calculated and set as the A.min value.
#'     To better determine how to set A.min see the help for ?filter_params().
#' @param adjust.dist Logical, should the p-value adjustment be performed
#'     on a per distance basis. i.e. The p-values at distance 1 will be grouped
#'     and the p-value adjustment will be applied. This process is repeated for
#'     each distance. The highest 15% of distancecs are grouped together i.e. 
#'     if you matrix has a maximum distance of 100, then distances 85-100 will be
#'     pooled together for p-value adjustment.
#' @param p.method The method for p-value adjustment. See ?p.adjust() help for
#'     options and more information. Defaults to "fdr". Can be set to "none" for
#'     no p-value adjustments.
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
#'     difference detection. Difference detection is performed by converting adjusted M 
#'     values to Z-scores. Any M value with a corresponding average expression level (A;
#'     mean of IF1 and IF2)
#'     less than the specified A.quantile is not considered for Z-score calculation. 
#'     This throws out the
#'     untrustworthy interactions that tend to produce false positives.
#'     The Z-scores are assumed to follow a roughly standard normal
#'     distribution and p-values are obtained. P-value adjustment for multiple testing 
#'     is then performed on a per distance basis (or on all p-values, optionally). 
#'     i.e. at each distance the vector of p-values corresponding to the
#'     interactions occuring at that distance have the selected multiple testing correction 
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
hic_compare <- function(hic.table, A.min = NA, adjust.dist = TRUE, p.method = 'fdr',
                        Plot = FALSE, Plot.smooth = TRUE,
                        parallel = FALSE, BP_param = bpparam()) {
  # check for correct input
  if (is(hic.table, "list")) {
    if ( sapply(hic.table, ncol) %>% min() < 13) {
      stop("Make sure you run hic_loess() on your hic.table before inputting it into hic_compare()")
    }
  } else {
    if (ncol(hic.table) < 13) {
      stop("Make sure you run hic_loess() on your hic.table before inputting it into hic_compare()")
    }
  }
  # check if single hic.table or list
  if (is.data.table(hic.table)) {
    hic.table <- list(hic.table)
  }
  # check A.min
  if (!is.na(A.min)) {
    if (!is.numeric(A.min)) {
      stop("Enter either NA or a numeric value > 1 for A.min")
    } else {
      if (A.min < 1) {
        stop("Enter either NA or a numeric value > 1 for A.min")
      }
    }
  }
  # if A.min = NA then automatically set it to the 10th percentile
  if (is.na(A.min)) {
    # get A quantiles 
    A_q10 <- sapply(hic.table, function(x) quantile(x$A, 0.1))
    A.min <- mean(A_q10) %>% ceiling()
  }
  
  # # check A.quantile
  # if (is.na(A.min) & !is.numeric(A.quantile)) {
  #   stop('Enter a numeric value for A.quantile or A.min')
  # }
  # if (is.na(A.min) & is.numeric(A.quantile)) {
  #   if (A.quantile < 0 | A.quantile >= 1) {
  #     stop('Enter a value between 0 and 1 for A.quantile')
  #   }
  # }
  
  # calculate z-scores
  if (parallel) {
    hic.table <- BiocParallel::bplapply(hic.table, .calc_z2, A.min = A.min) 
  } else {
    hic.table <- lapply(hic.table, .calc_z2, A.min = A.min) 
  }
  
  # adjust p-values
  if (adjust.dist) {
    if (parallel) {
      hic.table <- BiocParallel::bplapply(hic.table, .adjust_pval, Plot = Plot, p.method = p.method, p.smooth = Plot.smooth)
    } else {
      hic.table <- lapply(hic.table, .adjust_pval, Plot = Plot, p.method = p.method, p.smooth = Plot.smooth) 
    }
  } else {
    hic.table <- lapply(hic.table, function(x) {
      x[, p.adj := p.adjust(x$p.value, method = p.method)]
      return(x)
    })
    if (Plot) lapply(hic.table, function(x) {
      MD.plot2(x$adj.M, x$D, x$p.adj, smooth = Plot.smooth)
    })
  }
  
  
  # clean up if single hic.table
  if (length(hic.table) == 1) {
    hic.table <- hic.table[[1]]
  }
  return(hic.table)
}




# version where M values with A < thresh removed before z score calculations
# this version makes M have a lower standard deviation and thus higher z-scores
.calc_z2 <- function(hic.table, A.min, Plot = TRUE) {
  # threshold <- quantile((hic.table$A), quant, na.rm = TRUE)
  new_M <- hic.table$adj.M
  # if (is.na(A.min)) {
  #   # set M to NA to be ignored if A < quantile or IF1 < 1 or IF2 < 1
  #   new_M[hic.table$A < threshold | hic.table$adj.IF1 < 1 | hic.table$adj.IF2 < 1] <- NA
  # } else {
    # set M to be NA to be ignored if A < A.min or IF1 < 1 or IF2 < 1
    new_M[hic.table$A < A.min | hic.table$adj.IF1 < 1 | hic.table$adj.IF2 < 1] <- NA
  # }
  
  Z1 <- (new_M - mean(new_M, na.rm = TRUE)) / sd(new_M, na.rm = TRUE)
  hic.table[, Z := Z1]
  p.val <- 2*pnorm(abs(hic.table$Z), lower.tail = FALSE)
  # set P-values with NA value to 1
  p.val[is.na(p.val)] <- 1
  hic.table[, p.value := p.val]
  
  
  # if (Plot) MD.plot2(hic.table$adj.M, hic.table$D, hic.table$p.value)
  return(hic.table)
}


# adjust p-values in distance dependent manner
.adjust_pval <- function(hic.table, Plot, p.method, p.smooth) {
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
    x[, p.adj := p.adjust(p.value, method = p.method)]
    return(x)
  })
  # recombine into one table
  hic.table <- rbindlist(temp_list)
  
  if (Plot) {
    # if plotting ggplot MD plot then need to save as object and print it
    if (!p.smooth) {
      p1 <- MD.plot2(hic.table$adj.M, hic.table$D, hic.table$p.adj, smooth = p.smooth)
      print(p1)
    } else {
      MD.plot2(hic.table$adj.M, hic.table$D, hic.table$p.adj, smooth = p.smooth)
    }
  }
  return(hic.table)
}




# # version where z scores calculated first then z scores with A < thresh set to 0
# # this version means that M values have wider standard deviation thus lower Z-scores
# .calc_z <- function(hic.table, quant, A.min, Plot  = TRUE) {
#   # add average expression to table
#   # A <- (hic.table$adj.IF1 + hic.table$adj.IF2) / 2
#   # hic.table[, A := A]
#   threshold <- quantile((hic.table$A), quant, na.rm = TRUE)
#   Z1 <- (hic.table$adj.M - mean(hic.table$adj.M)) / sd(hic.table$adj.M)
#   # set z-scores where A < threshold to 0
#   if (is.na(A.min)) {
#     Z1[hic.table$A < threshold] <- 0
#     # set z-score to 0 if either adj.IF1 or adj.IF2 < 1
#     Z1[hic.table$adj.IF1 < 1 | hic.table$adj.IF2 < 1] <- 0
#   } else {
#     # Set Z scores where A < A.min to 0
#     Z1[hic.table$A < A.min] <- 0
#   }
#   
#   hic.table[, Z := Z1]
#   hic.table[, p.value := 2*pnorm(abs(Z), lower.tail = FALSE)]
#   # if (Plot) MD.plot2(hic.table$adj.M, hic.table$D, hic.table$p.value)
#   return(hic.table)
# }

