# These functions align alleles and betas

#' Find complement SNP
#'
#' @title complement_snp
#' @param x character
#' @return character, complement
#' @family align functions
#' @examples
#' complement_snp("A")
#' 
#' @author James Boocock
#' @keywords internal
complement_snp <- function(x){
    as = x =="A"
    ts = x == "T"
    gs = x == "G"
    cs = x == "C"
    ins = x == "I"
    dels = x == "D"
    x[as] = "T"
    x[ts] = "A"
    x[gs] = "C"
    x[cs] = "G"
    x[ins] = "NA"
    x[dels] = "NA"
    return(x)
}

#' Change indel format to I/D
#'
#' @title change_indels
#' @param data, a data frame with columns A1 and A2
#' @return data frame
#' @family align functions
#' 
#' @author Claudia Giambartolomei
#' @keywords internal
change_indels <- function(data) {
                data$A2[nchar(data$A1)>1] = "D"
                data$A1[nchar(data$A1)>1] = "I"
                data$A1[nchar(data$A2)>1] = "D"
                data$A2[nchar(data$A2)>1] = "I"
                return(data)
}

#' Match alleles to a reference and flip betas
#'
#' @title match_alleles
#' @param data, is a data.frame with columns A1 and A2, and reference alleles 
#' that need to match
#' @param A1.ref column names of ref alleles for A1 in the same data frame
#' @param A2.ref column names of ref alleles for A2 in the same data frame
#' @param flip logical, whether to flip the betas if the alleles are flipped
#' @return data frame
#' @family align functions
#' 
#' 
#' 
#' @author Claudia Giambartolomei, James Boocock
#' @keywords internal
match_alleles <- function(data, A1.ref="A1.ref", A2.ref="A2.ref",  A1.data = "A1", A2.data = "A2", BETA.data="BETA", flip = TRUE) {
   match_correct = data[,A1.ref] == data[,A1.data] & data[,A2.ref]== data[,A2.data]
   match_flip = data[,A1.ref] == data[,A2.data] & data[,A2.ref] == data[,A1.data]
   match_comp_one = data[,A1.ref] == complement_snp(data[,A1.data]) & data[,A2.ref]== complement_snp(data[,A2.data])
   match_comp_two = data[,A1.ref] == complement_snp(data[,A2.data]) & data[,A2.ref] == complement_snp(data[,A2.data])
   snp_allele_match = match_flip | match_correct | match_comp_one | match_comp_two
   message(sum(snp_allele_match), " SNPs out of ", length(snp_allele_match), " had the correct alleles, discarding SNPs without the correct alleles")
   if (flip) {
           if (any(which(match_flip)>0)) {
               data[match_flip, A1.data]=data[match_flip, A1.ref]
               data[match_flip, A2.data]=data[match_flip, A2.ref]
               data[match_flip, BETA.data]=-data[match_flip, BETA.data]
           }
   }
   removed = data[!snp_allele_match,]
   data = data[snp_allele_match,]
   return(list(data, removed))
  }