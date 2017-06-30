# Functions for matrix transformations
#' Transform a sparse upper triangular matrix to a full Hi-C contact matrix
#'
#' @export
#' @param sparse.mat A matrix in sparse upper triangular format.
#' @param hic.table Logical, is your sparse.mat a hic.table?
#' @param column.name Character, Required if hic.table set to TRUE; The column
#'     name of the hic.table that you want placed into the cells of the full matrix.
#'     i.e. IF1, or p.value.
#' @description sparse2full will transform a sparse upper triangular Hi-C matrix
#'     to a full Hi-C chromatin contact matrix. If you are entering a
#'     simple sparse matrix, i.e. there are only 3 columns leave
#'     hic.table = FALSE and column.name = NA. If you wish to transform  a Hi-C
#'     matrix in hic.table object format into a full matrix then set
#'     hic.table = TRUE. You will then need to specify the column name that you
#'     wish to be entered as the values for the cells in the full matrix using
#'     the column.name option.
#' @return A full Hi-C contact Matrix.
#' @examples
#' data('NHEK.chr22')
#' full.mat <- sparse2full(NHEK.chr22)
#'
# input sparse upper triangular matrix file containing 3 columns. cols
# 1, 2 are bin location, col 3 is frequency of interaction function
# will construct the full matrix from the sparse matrix notation
sparse2full <- function(sparse.mat, hic.table = FALSE, column.name = NA) {
  # get the min and max bins if the matrix is a simple sparse matrix
  if (ncol(sparse.mat) == 3) {
    sparse.mat <- as.matrix(sparse.mat)
    bins <- unique(c(sparse.mat[, 1], sparse.mat[, 2]))
    bins <- as.numeric(bins)
    bins <- bins[order(bins)]
    bin.size <- min(diff(bins))
    min.bin <- min(bins)
    max.bin <- max(bins)
    # set sequence of rows/cols for full matrix
    cols <- seq(min.bin, max.bin, by = bin.size)
    # set up matrix to be correct size and row/col names
    mat <- matrix(nrow = length(cols), ncol = length(cols))
    rownames(mat) <- cols
    colnames(mat) <- cols
  }
  # if the matrix is a hic.table
  if (ncol(sparse.mat) > 3 & !hic.table) {
    stop("Enter a sparse matrix with 3 columns or if you are trying to enter a
         hic.table please set hic.table = TRUE")
  }
  if (ncol(sparse.mat) > 3) {
    sparse.mat <- as.matrix(sparse.mat)
    bins <- unique(c(sparse.mat[, 2], sparse.mat[, 5]))
    bins <- as.numeric(bins)
    bins <- bins[order(bins)]
    bin.size <- min(diff(bins))
    min.bin <- min(bins)
    max.bin <- max(bins)
    # set sequence of rows/cols for full matrix
    cols <- seq(min.bin, max.bin, by = bin.size)
    # set up matrix to be correct size and row/col names
    mat <- matrix(nrow = length(cols), ncol = length(cols))
    rownames(mat) <- cols
    colnames(mat) <- cols
  }
  if (hic.table) {
    if (is.na(column.name)) {
      stop("Enter a value for column.name")
    }
    # reconstruct matrix using by converting bin names to matrix cell
    # locations
    # for (count in 1:dim(sparse.mat)[1]) {
    #   row_num <- which(as.numeric(sparse.mat[count, 2]) == cols)
    #   col_num <- which(as.numeric(sparse.mat[count, 5]) == cols)
    #   mat[row_num, col_num] <- as.numeric(sparse.mat[count, column.name])
    #   mat[col_num, row_num] <- as.numeric(sparse.mat[count, column.name])
    # }
    col.val <- which(colnames(sparse.mat) == column.name)
    sparse.mat <- sparse.mat[, c(2, 5, col.val)]
    sparse.mat <- apply(sparse.mat, 2, as.numeric)
    # match bin names to column/row number
    sparse.mat[,1] <- match(sparse.mat[,1], cols)
    sparse.mat[,2] <- match(sparse.mat[,2], cols)
    # populate matrix
    rcix <- sparse.mat[,c(1,2)]
    mat[rcix] <- sparse.mat[, 3]
    rcix <- rcix[, c(2,1)]
    mat[rcix] <- sparse.mat[, 3]

  } else {
    # reconstruct matrix using by converting bin names to matrix cell
    # locations
    # for (count in 1:dim(sparse.mat)[1]) {
    #   row_num <- which(sparse.mat[count, 1] == cols)
    #   col_num <- which(sparse.mat[count, 2] == cols)
    #   mat[row_num, col_num] <- sparse.mat[count, 3]
    #   mat[col_num, row_num] <- sparse.mat[count, 3]
    # }
    # match bin names to column/row number
    sparse.mat[,1] <- match(sparse.mat[,1], cols)
    sparse.mat[,2] <- match(sparse.mat[,2], cols)
    # populate matrix
    rcix <- sparse.mat[, 1:2]
    mat[rcix] <- sparse.mat[,3]
    rcix <- rcix[,c(2,1)]
    mat[rcix] <- sparse.mat[,3]
  }
  # replace NAs in matrix with 0
  mat[is.na(mat)] <- 0
  message(paste0("Matrix dimensions: ", dim(mat)[1], "x", dim(mat)[2]))
  return(mat)
}

#' Transfrom a full Hi-C contact matrix to a sparse upper triangular matrix
#'
#' @export
#' @param mat A matrix. Must have column names equal to the start location for
#'     each bin. i.e. for a 6x6 Hi-C matrix where the first region starts at 0 kb and
#'     the final region starts at 500KB and the resolution is 100kb, the column names
#'     of the matrix should be as follows:
#'     colnames(mat) = c(0, 100000, 200000, 300000, 400000, 500000)
#' @return A sparse upper triangular matrix.
#' @examples
#' m <- matrix(1:100, 10, 10)
#' colnames(m) <- 1:10
#' sparse <- full2sparse(m)
#' sparse
#'
# function to take full matrix and transform to upper triangular sparse
# matrix matrix should have column/row names corresponding to the start
# location of the region represented by that cell in the matrix
full2sparse <- function(mat) {
  if (is.null(colnames(mat))) {
    stop("Please set the column names of the matrix to the start location for
         each bin. See ?full2sparse")
  }
  regions <- as.numeric(colnames(mat))
  sparse <- as.data.table(expand.grid(regions, regions))
  up.triangle <- c(upper.tri(mat, diag = TRUE))
  sparse <- cbind(sparse, c(mat), up.triangle)
  colnames(sparse) <- c("region1", "region2", "IF", "up.triangle")
  sparse <- subset(sparse, up.triangle == TRUE & IF != 0)
  sparse <- sparse[, 1:3, with = FALSE]
  return(sparse)
}


#' Transform a .cool file to a sparse upper triangular matrix for input into
#'     hic_loess
#'
#' @export
#' @param cooler The plain text file from a .cool file loaded into an R
#'     data.frame object. See vignette for more details.
#'
#' @details The .cool format is linked a database of Hi-C experiments and
#'     allows access to many sets of Hi-C data which can be found
#'     at the Index of Coolers  \url{ftp://cooler.csail.mit.edu/coolers}. Once a .cool
#'     file is dumped into a contact matrix in plain text it can be
#'     read into R.  This function provides a method for converting the .cool matrix
#'     into a sparse upper triangular matrix ready to be entered into
#'     hic_loess.
#'
#' @return A Sparse upper triangular matrix or a list of sparse upper triangular matrices.
#'     If the .cool file contains data for more than one chromosome
#'     The function will split the data up into a list of matrices, one per chromosome.
#'
#' @examples
#' data('cooler')
#' sparse <- cooler2sparse(cooler)
#' head(sparse)

cooler2sparse <- function(cooler) {
  # convert to data.table
  cooler <- as.data.table(cooler)
  colnames(cooler) <- c("chr1", "start1", "end1", "chr2", "start2", "end2",
                        "IF")
  # subset to only include cis interactions since HiCloess doesn't work
  # with trans interaction if data for more than 1 chromosome entered
  if (is.factor(cooler$chr1)) {
    chroms <- levels(cooler$chr1)
  }
  if (is.character(cooler$chr1)) {
    chroms <- unique(cooler$chr1)
  }
  if (length(chroms) > 1) {
    cooler <- subset(cooler, chr1 == chr2)
    # split up data into list of sparse matrices for each chromosome
    sparse.list <- vector("list", length(chroms))
    for (i in seq_along(chroms)) {
      temp <- subset(cooler, chr1 == chroms[i])
      temp <- temp[, c("start1", "start2", "IF"), with = FALSE]
      colnames(temp) <- c("region1", "region2", "IF")
      temp <- subset(temp, IF != 0)
      sparse.list[[i]] <- temp
    }
    names(sparse.list) <- chroms
    return(sparse.list)
  } else {
    temp <- temp[, c("start1", "start2", "IF"), with = FALSE]
    colnames(temp) <- c("region1", "region2", "IF")
    temp <- subset(temp, IF != 0)
    return(temp)
  }
}
