#' Total sum normalization for a list of hic.table objects
#'
#' @export
#' @param hic.list A list of hic.table objects created by the
#'     create.hic.table function. When creating the hic.list
#'     to be entered into this function you must set the
#'     scale option to FALSE.
#' @details This function will scale the IFs based on the total sum
#'     of the counts for the genome instead of on a per
#'     chromosome basis as done in the create.hic.table function
#'     when the scale option is set to TRUE.
#'     The idea behind this function to
#'     preserve more local CNV differences while still accounting
#'     for technical biases which can cause the total read counts
#'     to differ between sequencing runs.
#'
#' @return A list of hic.table objects.
#'
#' @examples
#' data('HMEC.chr22')
#' data('NHEK.chr22')
#' data('HMEC.chr10')
#' data('NHEK.chr10')
#' hic.table1 <- create.hic.table(HMEC.chr22, NHEK.chr22,
#'     chr = 'chr22', scale = FALSE)
#' hic.table2 <- create.hic.table(HMEC.chr10, NHEK.chr10,
#'     chr = 'chr10', scale = FALSE)
#' hic.list <- list(hic.table1, hic.table2)
#' scaled_list <- total_sum(hic.list)


total_sum <- function(hic.list) {
  # check for valid input
  if (!is.list(hic.list)) {
    stop('Enter a list of hic.table objects created by the create.hic.table function')
  }
  if (length(hic.list) < 2) {
    stop('Enter a list of more than 1 hic.table objects. If you only want
         chromosome specific scaling just use the scale option in create.hic.table')
  }
  # check for integer values of IF1 and IF2
  temp <- data.table::rbindlist(hic.list)
  IF1_whole <- sapply(temp$IF1, .is.whole)
  IF2_whole <- sapply(temp$IF2, .is.whole)
  if (sum(IF1_whole) != length(IF1_whole) | sum(IF2_whole) != length(IF2_whole)) {
    stop('When creating the hic.table objects set scale = FALSE')
  }
  # calculate scale factor
  scale.factor <- sum(temp$IF2)/sum(temp$IF1)
  message('Scaling factor = ', scale.factor, ' Total counts IF1: ', sum(temp$IF1), ' Total counts IF2: ', sum(temp$IF2))
  # update table
  temp[, IF2 := IF2 / scale.factor]
  if (sum(temp$IF1 == 0) > 0 | sum(temp$IF2 == 0) > 0) {
    temp[, M := log2((IF2 + 1)/(IF1 + 1))]
  } else {
    temp[, M := log2(IF2 / IF1)]
  }
  # recreate list
  chr_list <- unique(temp$chr1)
  new_list <- vector(mode = "list", length = length(chr_list))
  for (i in seq_along(chr_list)) {
    new_list[[i]] <- subset(temp, chr1 == chr_list[[i]])
  }

  # check for CNV by comparing chromosome sums
  lapply(new_list, .check_CNV)

  return(new_list)
}



# function to check if a value is a whole number vs decimal
.is.whole <- function(a) {
  (is.numeric(a) && floor(a)==a) ||
    (is.complex(a) && floor(Re(a)) == Re(a) && floor(Im(a)) == Im(a))
}

# function to check for CNVs by comparing total sums of IFs
.check_CNV <- function(hic.table) {
  chr <- hic.table$chr1[1]
  res <- t.test(hic.table$IF1, hic.table$IF2)
  if (res$p.value < 0.05) {
    message(chr, ' may have a Copy Number Variation. Check the MD plots.')
  }
}
