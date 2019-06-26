#' Read a .cool file into R and output the data in BEDPE format
#' 
#' @param path The path to a .cool file on your disk.
#' @details .cool files are HDF5 containers that store Hi-C data.
#'     Many public Hi-C datasets are available in .cool format
#'     on the mirnylab ftp site \url{ftp://cooler.csail.mit.edu/coolers}.
#'     To use these files in HiCcompare simply download the .cool file
#'     and read it into R using this function. This function will dump
#'     the contents of the file and format them into BEDPE format in R.
#'     The resulting object cant then be used in HiCcompare.
#' @return A list with two items. Item 1, "cis"
#'     contains the intra-chromosomal contact matrices,
#'     one per chromosome.
#'     Item 2, "trans" contains the inter-chromsomal
#'     contact matrix.
#'    
#' @examples 
#' \dontrun{
#'  dat <- cooler2bedpe(path = "path/to/cool/file.cool")
#' }
#' @export
#' @importFrom rhdf5 h5dump
#' 
cooler2bedpe <- function(path) {
  # dump contents of cooler file
  dump <- rhdf5::h5dump(path)
  # genomic coordinates are stored in $bins
  # IFs are stored in $pixels where bin1_id and bin2_id correspond to $bins
  # make ids
  ids <- data.frame(chr = dump$bins$chrom, start = dump$bins$start, end = dump$bins$end, id = seq(1, length(dump$bins$chrom), by = 1))
  # make sparse matrix
  mat <- data.frame(bin1 = dump$pixels$bin1_id, bin2 = dump$pixels$bin2_id, IF = dump$pixels$count)
  # use hicpro2bedpe to convert to useable format
  bedpe <- hicpro2bedpe(mat, ids)
  return(bedpe)
}
