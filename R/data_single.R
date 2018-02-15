#' List including three data frames with summary statistics for three traits
#' 
#'
#' A list of data frames containing summary statistics from three indendent datasets,
#' for example GWAS, expression QTL, methylation QTL
#'
#' @format A list of three data frames, each with 831 rows and 7 variables:
#' \describe{
#'   \item{SNP}{SNP names, common to all data frames}
#'   \item{BETA}{Effect size estimate from regression, if OR must give the logOR}
#'   \item{SE}{Standard Error}
#'   \item{z}{Z-score}
#'   \item{PVAL}{P-value}
#'   \item{N}{Sample size}
#'   \item{MAF}{Minor allele frequency}
#'   ...
#' }
"data_single"