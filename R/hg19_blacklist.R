#' BED file for hg19 blacklisted regions
#'
#' A BED file with regions which may be
#'     excluded from your Hi-C data analysis.
#'
#' @format A data.frame with 3 columns and 1649 rows:
#'     \describe{
#'     \item{chr}{chromosome for the region}
#'     \item{start}{start location of the region}
#'     \item{end}{end location of the region}
#'  }
#' @source Data from UCSC
#'     \url{http://genome.ucsc.edu/cgi-bin/hgFileUi?db=hg19&g=wgEncodeMapability}
#'
#' @return A data.frame

#'
"hg19_blacklist"
