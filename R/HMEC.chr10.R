#' Hi-C data from HMEC cell line - chromosome 10 at 500kb resolution
#'
#' A sparse upper triangular matrix containing the interacting regions and the
#'     corresponding Interaction Frequency (IF).
#'
#'
#'
#' @format A matrix with 35386 rows and 3 columns:
#'     \describe{
#'     \item{region1}{The first interacting region - corresponds to the row
#'          name of a full Hi-C contact matrix}
#'     \item{region2}{The second interacting region - corresponds to the
#'          column name of a full Hi-C contact matrix}
#'     \item{IF}{The Interaction Frequency - number of read counts found
#'          for the interaction between region1 and region2}
#'     }
#'
#' @source Data from the Aiden Lab. See their website at
#'     \url{http://www.aidenlab.org/}
#'      Or the the GEO link to download the data
#'     \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525}
#'
#' @return A matrix

#'
"HMEC.chr10"
