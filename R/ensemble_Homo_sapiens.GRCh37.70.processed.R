#' Bed file with Ensembl genes, start position; stop position; for genome-wide analysis.
#' If using this file, it is recemmended to add +/- 50000 for methylation analysis, or 200000 for eQTL/GWAS analysis.
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
"ensemble_Homo_sapiens.GRCh37.70.processed"