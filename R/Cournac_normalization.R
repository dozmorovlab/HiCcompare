#' SCN normalization from Cournac 2012
#'
#' @export
#' @param a The matrix to be normalized. Any cols/rows that sum to
#'    0 will be removed before normalization.
#' @param max.iter maximum number of iterations to be performed
#' @details Performs Sequential Component Normalization as described by
#'     Cournac. Coded using details in the manuscript.
#'     Cournac A, Marie-Nelly H, Marbouty M, Koszul R, Mozziconacci J.
#'     Normalization of a chromosomal contact map.
#'     BMC Genomics. 2012;13: 436. doi:10.1186/1471-2164-13-436
#' @return An SCN normalized matrix
#'
#' @examples
#' m <- matrix(rpois(100, 5), 10, 10)
#' SCN(m)

SCN = function(a, max.iter = 10) {
  # remove 0 columns from matrix
  zeros = unique(which(colSums(a) == 0), which(rowSums(a) == 0))
  if (length(zeros) > 0) {
    a = a[-zeros, -zeros]
    message(paste0('Cols/Rows removed: '))
    message(paste(" ", zeros, sep = ""))
  }
  diff = 1
  iterations = 0
  while(iterations < max.iter) {
    a.old <- a
    a.new <- a

    n.rows <- dim(a)[1]; n.cols <- dim(a)[2]
    # Euclidean normalization of matrix
    a.col.norms <- apply(a,2, function(x) sqrt(sum(x^2)) )
    # a.new <- t(t(a) / a.col.norms) # old version of code. Not sure why matrices transposed
    a.new <- (a) / a.col.norms # new

    a.row.norms <- apply(a.new, 1, function(x) sqrt(sum(x^2)) )
    # a <- t(t(a.new) / a.row.norms) # old version
    a <- a.new / a.row.norms # new

    diff = a.old - a
    diff = c(diff)
    iterations = iterations + 1
  }
  return(a)
}
