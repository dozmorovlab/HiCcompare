#' Function to split hic.table into 2 subsets at the centromere
#'
#' @export
#' @param hic.table A hic.table or list of hic.table objects
#' @return A list of hic.tables in which the matrices have been split at the centromere
#' @examples
#' data('HMEC.chr22')
#' data('NHEK.chr22')
#' hic.table <- create.hic.table(HMEC.chr22, NHEK.chr22, chr = 'chr22')
#' split <- split_centromere(hic.table)

split_centromere <- function(hic.table) {
  # check if list or single hic.table
  if (is.data.table(hic.table)) hic.table = list(hic.table)
  # apply centromere splitting
  tmp <- lapply(hic.table, .split_cent)
  # new_tables <- lapply(seq_along(new_tables), function(x) new_tables[[x]])
  new_tables <- tmp[[1]]
  if (length(tmp) > 1) {
    for (i in 2:length(tmp)) {
      new_tables <- c(new_tables, tmp[[i]])
    }
  }
  return(new_tables)
}



# centromere_locations <- read.table("C:/VM_Shared/DNA_int_analysis/HiCdiff/data/hg19_centromere.bed", header = FALSE)

.split_cent <- function(hic.table) {
  # get chromosome for current table
  chrom <- hic.table$chr1[1]
  cent_start <- centromere_locations[which(centromere_locations[,1] == chrom), 2]
  cent_end <- centromere_locations[which(centromere_locations[,1] == chrom), 3]
  arm1 <- subset(hic.table, start1 <= cent_start & start2 <= cent_start)
  arm2 <- subset(hic.table, start1 > cent_start & start2 > cent_start)
  arm1[, ':='(chr1 = paste0(chr1, '.1'), chr2 = paste0(chr2, '.1'))]
  arm2[, ':='(chr1 = paste0(chr1, '.2'), chr2 = paste0(chr2, '.2'))]
  if (nrow(arm1) < 1) {
    result <- list(arm2)
  }
  if (nrow(arm2) < 1) {
    result <- list(arm1)
  }
  if (nrow(arm1) > 0 & nrow(arm2) > 0) {
    result <- list(arm1, arm2)
  }
  return(result)
}



# wrapper function



