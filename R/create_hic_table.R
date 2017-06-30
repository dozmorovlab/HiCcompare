#' Create hic.table object from a sparse upper triangular Hi-C matrix
#'
#' @export
#' @param sparse.mat1 Required, sparse upper triangular Hi-C matrix, 7 column
#'     BEDPE format of the upper triangle of the matrix, OR InteractionSet
#'     object with the genomic ranges of the interacting regions for the upper
#'     triangle of the Hi-C matrix and a single metadata column containing the
#'     interaction frequencies for each interacting pair for the first dataset
#'     you wish to jointly normalize.
#' @param sparse.mat2 Required, sparse upper triangular Hi-C matrix, 7 column
#'     BEDPE format of the upper triangle of the matrix, OR InteractionSet
#'     object with the genomic ranges of the interacting regions for the upper
#'     triangle of the Hi-C matrix and a single metadata column containing the
#'     interaction frequencies for each interacting pair for the second dataset
#'     you wish to jointly normalize.
#' @param chr The chromosome name for the matrices being entered i.e 'chr1' or
#'     'chrX'. Only needed if using sparse upper triangular matrix format.
#'     If using BEDPE format leave set to NA.
#' @param scale Logical, should scaling be applied to the matrices to adjust
#'     for total read counts. If TRUE the IFs of the second sparse
#'     matrix will be adjusted as follows: IF2_scaled = IF2 / (sum(IF2)/sum(IF1)).
#' @param include.zeros Logical, If set to TRUE the function will include pairwise
#'     interactions where one of the interaction frequencies is 0.
#' @param subset.dist Should the matrix be subset to only include interactions up
#'     to a user specified matrix unit distance? i.e. to only include
#'     the cells of the matrix which are at a unit distance less than or equal to
#'     100 set \code{subset.dist = 100}. Subsetting the matrix by distance
#'     will cut out any interactions occuring at a unit distance greater than the
#'     specified value. This could be used to speed up computation time or if
#'     there is only interest in the interactions occuring below a specific
#'     distance in the matrix. Warning: If you subset by distance do NOT to
#'     transform the subsetted hic.table into a full matrix using `sparse2full`.
#'     If you plan on trasforming the matrix to a full contact matrix
#'     use subset.index instead.
#' @param subset.index Should the matrix be subset by a user specified distance?
#'     Input as a vector of 4 numbers (i.start, i.end, j.start, j.end).
#'     i.e. to only include a subset of the matrix with row numbers 20 <= i <= 40
#'     and column numbers 30 <= j <= 50 set as \code{subset.index = c(20, 40, 30, 50)}.
#'     This can be used to speed up computation time if only a subset of the matrix is
#'     of interest. The indices used here correspond to the indices of the full
#'     Hi-C contact matrix. The `sparse2full` function can be used to view the full
#'     contact matrix and make a decision about subsetting based on index.
#'
#'@details This function is used to transform two sparse upper triangular Hi-C matrices
#'    into an object usable in the \code{hic_loess} function.
#'    Sparse upper triangular Hi-C matrix format is typical of the Hi-C data available
#'    from the Aiden Lab \url{http://www.aidenlab.org/}. If you have a full
#'    Hi-C contact matrix, first transform it to sparse upper triangular format using
#'    the \code{full2sparse} function. Sparse matrices should have 3 columns
#'    in the following order: Start location of region 1, Start location of region 2,
#'    Interaction Frequency. Matrices in 7 column BEDPE format should
#'    have 7 columns in the following order: Chromosome name of the first region,
#'    Start location of first region, End location of first region,
#'    Chromosome name of the second region, Start location of the second region,
#'    End location of the second region, Interaction Frequency. Please enter either
#'    two sparse matrices or two matrices in 7 column BEDPE format or two
#'    InteractionSet objects; do not mix and match.
#'
#' @return A hic.table object.
#' @examples
#' # Create hic.table object using included Hi-C data in sparse upper
#' # triangular matrix format
#' data('HMEC.chr22')
#' data('NHEK.chr22')
#' hic.table <- create.hic.table(HMEC.chr22, NHEK.chr22, chr = 'chr22')
#' # View result
#' hic.table


create.hic.table <- function(sparse.mat1, sparse.mat2, chr = NA, scale = TRUE,
                             include.zeros = FALSE, subset.dist = NA, subset.index = NA) {
  if (!is.na(subset.dist) & !is.na(subset.index[1])) {
    stop("Enter a value for only one of the subsetting options")
  }
  # code to let it accept a GInteractions input for sparse.mat1 & sparse.mat2
  if ( (is(sparse.mat1, "GInteractions") & !is(sparse.mat2, "GInteractions")) |
       (!is(sparse.mat1, "GInteractions") & is(sparse.mat2, "GInteractions")) ) {
    stop("Make sure the classes of the sparse matrices match")
  }
  if (is(sparse.mat1, "GInteractions") & is(sparse.mat2, "GInteractions")) {
    sparse.mat1 <- as.data.table(sparse.mat1)
    sparse.mat2 <- as.data.table(sparse.mat2)
    chr <- sparse.mat1$seqnames1[1]
    sparse.mat1 <- sparse.mat1[, c(2, 7, 11), with = FALSE]
    sparse.mat2 <- sparse.mat2[, c(2, 7, 11), with = FALSE]
  }
  # check if sparse matrices are 3 column sparse upper triangular format
  # or 7 column BEDPE format
  if (ncol(sparse.mat1) == 7 & ncol(sparse.mat2) == 7) {
    chr <- as.character(sparse.mat1[1, 1])
    sparse.mat1 <- sparse.mat1[, c(2, 5, 7)]
    sparse.mat2 <- sparse.mat2[, c(2, 5, 7)]
  }
  if (ncol(sparse.mat1) == 3 & ncol(sparse.mat2) == 3) {
    if (is.na(chr)) {
      stop("Enter a value for chr")
    }
    # convert matrices to data.tables
    hic.table1 <- data.table::as.data.table(sparse.mat1)
    colnames(hic.table1) <- c("region1", "region2", "IF1")
    hic.table2 <- data.table::as.data.table(sparse.mat2)
    colnames(hic.table2) <- c("region1", "region2", "IF2")
    # define bin.size
    bins <- unique(c(hic.table1$region1, hic.table1$region2, hic.table2$region1,
                     hic.table2$region2))
    bins <- bins[order(bins)]
    bin.size <- min(diff(bins))
    # merge tables
    if (include.zeros) {
      hic.table <- data.table::as.data.table(dplyr::full_join(hic.table1,
                                                              hic.table2, by = c(region1 = "region1",
                                                                                 region2 = "region2")))
      hic.table$IF1[is.na(hic.table$IF1)] <- 0
      hic.table$IF2[is.na(hic.table$IF2)] <- 0
    } else {
      hic.table <- data.table::as.data.table(dplyr::inner_join(hic.table1,
                                                               hic.table2, by = c(region1 = "region1",
                                                                                  region2 = "region2")))
    }
    # calculate matrix distance of all interaction
    if (scale) {
      hic.table[, `:=`(D, abs(region2 - region1)/bin.size)]
      scale.factor <- sum(hic.table$IF2)/sum(hic.table$IF1)
      hic.table[, `:=`(IF2, IF2/scale.factor)]
      if (include.zeros) {
        hic.table[, `:=`(M, log2((IF2 + 1)/(IF1 + 1)))]
      } else {
        hic.table[, `:=`(M, log2((IF2)/(IF1)))]
      }
    } else {
      hic.table[, `:=`(D, abs(region2 - region1)/bin.size)]
      if (include.zeros) {
        hic.table[, `:=`(M, log2((IF2 + 1)/(IF1 + 1)))]
      } else {
        hic.table[, `:=`(M, log2((IF2)/(IF1)))]
      }
    }
    if (!is.na(subset.dist[1])) {
      hic.table <- subset(hic.table, D <= subset.dist)
    }
    if (!is.na(sum(subset.index))) {
      if (length(subset.index) != 4) {
        stop("subset.index must be a vector of 4 numbers: (i.start, i.end, j.start, j.end)
             where they correspond to i and j indices of the matrix that you want to subset to")
      }
      # add i and j coordinates for each entry in table
      regions <- unique(c(hic.table1$region1, hic.table$region2))
      regions <- regions[order(regions)]
      hic.table[, `:=`(i, match(region1, regions))]
      hic.table[, `:=`(j, match(region2, regions))]
      hic.table <- subset(hic.table, i >= subset.index[1] & i <=
                            subset.index[2] & j >= subset.index[3] & j <= subset.index[4])
      }
    new.hic.table <- data.table(chr1 = chr, start1 = hic.table$region1,
                                end1 = hic.table$region1 + bin.size, chr2 = chr, start2 = hic.table$region2,
                                end2 = hic.table$region2 + bin.size, IF1 = hic.table$IF1, IF2 = hic.table$IF2,
                                D = hic.table$D, M = hic.table$M)
  } else {
    stop("Enter both sparse matrices in the same format; either 7 column BEDPE or 3 column sparse upper triangular matrix")
  }
  return(new.hic.table)
}

