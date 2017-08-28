#' BED file for hg38 blacklisted regions
#'
#' A BED file with regions which may be
#'     excluded from your Hi-C data analysis.
#'
#' @format A data.frame with 3 columns and 38 rows:
#'     \describe{
#'     \item{chr}{chromosome for the region}
#'     \item{start}{start location of the region}
#'     \item{end}{end location of the region}
#'  }
#' @source Data from
#'     \url{http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg38-human/}
#'
#' @return A data.frame

#'
"hg38_blacklist"
