#' Function to remove centromere columns and rows from a full Hi-C
#'     contact matrix
#'
#' @export
#' @param mat A full Hi-C matrix
#'
#' @return A list of (1) the column/row numbers of the centromere and
#'     (2) the Hi-c matrix with the centromere removed
#' @examples
#' m <- matrix(rpois(100, 5), 10, 10)
#' m[5,] <- 0
#' m[,5] <- 0
#' remove_centromere(m)

remove_centromere <- function(mat) {
  cnames <- NA
  if (length(colnames(mat)) > 0) {
    cnames <- colnames(mat)
  }
  # get row and col sums for matrix
  col.sum <- colSums(mat)
  row.sum <- rowSums(mat)
  centromere.col <- which(col.sum == 0)
  centromere.row <- which(row.sum == 0)
  centromere <- intersect(centromere.col, centromere.row)
  # break up into list by what cols are consecutive
  breaks <- c(0, which(diff(centromere) != 1), length(centromere))
  if (length(breaks) > 2) {
    consecutive.list <- sapply(seq(length(breaks) - 1),
                               function(i) centromere[(breaks[i] + 1):breaks[i + 1]])
    # find out which consecutive set of cols is the maximum; this should
    # represent the centromere
    lengths <- sapply(consecutive.list, length)
    centromere <- consecutive.list[[which(lengths == max(lengths))]]
  }
  # remove centromere rows/cols from matrix
  if (length(centromere) > 0) {
    new.mat <- mat[-centromere, -centromere]
    if (length(cnames) > 1) {
      colnames(new.mat) <- cnames[-centromere]
    }
  } else {
    new.mat <- mat
  }
  result <- list(centromere_loc = centromere, mat = new.mat)
  return(result)
}
