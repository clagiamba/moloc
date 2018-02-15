#' Bed file needed for genome-wide analysis, goes with data_genome
#' 
#'
#' A tab delimited text file with ProbeID and regions to consider
#' for the moloc analysis
#'
#' @format A data frame with 3 variables:
#' \describe{
#'   \item{ProbeID}{Name of Genes or methylation Probes, as listed in the eqtl/mqtl summary stats}
#'   \item{CHR}{CHR}
#'   \item{START}{START of region}
#'   \item{STOP}{STOP of region}
#'   ...
#' }
"bed"