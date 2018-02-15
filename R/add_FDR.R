##' FDR using PPA
##'
##' Add a column with FDR of colocalization results
##' 
##' @title add_fdr
##' @param column name of posterior probability to use for FDR
##' @return data frame with added column
##' @author Claudia Giambartolomei
#' @keywords internal
add_fdr <- function(df, snpcol = "PP4") {
    df$PEP <- 1 - df[,snpcol]
    df <- df[order(df$PEP, decreasing=F),]
    df$FDR <- cumsum(df$PEP)/seq(df$PEP)
    message("There are ", nrow(df[df$FDR<=0.05,]), " results at FDR 5%")
    return(df)
}