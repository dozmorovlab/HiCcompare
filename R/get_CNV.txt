#' Function to get the locations of copy number variations

#' @param path2bam The path to the folder containing the .bam files for the data
#' @param out.file The prefix for the filename which will be exported.
#' @param bin.size The bin.size (resolution) of the HiC data you will be using in kbp. e.g. for 1MB
#'     resolution enter 1000, for 500kb resolution enter 500. The available bin sizes
#'     for hg19 1, 5, 10, 15, 30, 50, 100, 500, and 1000 kbp
#' @param genome Character string for the genome. Defaults to 'hg19'.
#' @param CNV.level The value require to call a region a CNV. Should be either
#'     1, or 2. This setting will determine what level of CNV will be returned
#'     in the resulting BED file. 1 indicates any copy number variation will be
#'     returned. 2 indicates only deletions and amplifications will be returned.
#'     Defaults to 2.
#'
#' @details This function utilizeds QDNAseq to detect CNV (copy number variations) in your
#'     Hi-C data. You must have the .bam files in order to use this function. Place the
#'     .bam files in a folder and set path2bam as the path to this folder. The function
#'     will then perform CNV detection and export a .txt file containing the CNV
#'     calls for each region of the genome. In addition the results will be returned
#'     from the function as a data.frame. It is reccomended to exclude the regions
#'     with extreme CNV when using HiCcompare. Regions can be excluded using the
#'     create.hic.table() function. For more information about the methods
#'     used in this function see QDNAseq here:
#'     \url{https://bioconductor.org/packages/release/bioc/html/QDNAseq.html}
#'
#' @return A data.frame containing the regions for CNVs at the level
#'     specified with the CNV.level option.
#'     A data.frame containing the calls for all regions will also be
#'     written to the disk using the out.file prefix.
#'     The first 4 columns of the data.frame contain location information
#'     for the genomic regions. All following columns contain the CNV calls
#'     for each .bam file. CNV calls are defined as follows: -2 = deletion,
#'     -1 = loss, 0 = normal, 1 = gain, 2 = amplification. Additionally,
#'     a .bed file will be written to disk containing the same values
#'     as the returned data.frame.
#'
#' @examples
#'\dontrun{
#' cnv_calls <- get_CNV(path2bam = 'C:/bamfiles')
#'}
# @import QDNAseq
# @export

get_CNV <- function(path2bam, out.file = 'CNV_calls', bin.size = 1000, genome = 'hg19',
                    CNV.level = 2) {
  if (CNV.level != 1 & CNV.level != 2) {
    stop('Enter a value of 1 or 2 for CNV.level')
  }
  # get annotations
  bins <- QDNAseq::getBinAnnotations(binSize = bin.size, genome = genome)
  readCounts <- QDNAseq::binReadCounts(bins, path = path2bam)
  # do their filtering steps
  readCountsFiltered <- QDNAseq::applyFilters(readCounts, residual = TRUE, blacklist = TRUE)
  readCountsFiltered <- QDNAseq::estimateCorrection(readCounts)
  # make the copynumber object
  copyNumbers <- QDNAseq::correctBins(readCountsFiltered)
  copyNumbersNormalized <- QDNAseq::normalizeBins(copyNumbers)
  copyNumbersSmooth <- QDNAseq::smoothOutlierBins(copyNumbersNormalized)
  copyNumbersSegmented <- QDNAseq::segmentBins(copyNumbersSmooth, transformFun="sqrt")
  copyNumbersCalled <- QDNAseq::callBins(copyNumbersSegmented)
  # export calls
  QDNAseq::exportBins(copyNumbersCalled, format="tsv", type = 'calls',
                      logTransform = FALSE, file = paste0(out.file, '.txt'))
  res <- read.table(paste0(out.file, '.txt'), header = TRUE)
  # return data.frame with results from exported file
  res$chromosome <- paste0('chr', res$chromosome)
  write.table(res, file = paste0(out.file, '.txt'), sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)

  # get which rows are at or above CNV.level
  tmp <- res[, 5:ncol(res)] %>% as.matrix()
  above_level <- apply(tmp, 2, function(x) ifelse(abs(x) >= CNV.level, TRUE, FALSE))
  above_level <- apply(above_level, 1, sum)
  above_level <- ifelse(above_level > 0, TRUE, FALSE)
  if (sum(above_level) > 0) {
    BED <- res[above_level, 2:4]
    write.table(BED, file = paste0(out.file, '.bed'), sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
  } else {
    message("No CNVs at or greater than chosen level")
    BED <- NA
  }
  return(BED)
}
