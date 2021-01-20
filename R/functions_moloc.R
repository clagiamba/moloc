#library(mvtnorm)

##' Estimate trait standard deviation given vectors of variance of coefficients,  MAF and sample size
##'
##' Estimate is based on var(beta-hat) = var(Y) / (n * var(X))
##' var(X) = 2*maf*(1-maf)
##' so we can estimate var(Y) by regressing n*var(X) against 1/var(beta)
##' 
##' @title Estimate trait variance, internal function
##' @param vbeta vector of variance of coefficients
##' @param maf vector of MAF (same length as vbeta)
##' @param n sample size
##' @return estimated standard deviation of Y
##' 
##' @author Chris Wallace
sdY.est <- function(vbeta, maf, n) {
  oneover <- 1/vbeta
  nvx <- 2 * n * maf * (1-maf)
  m <- lm(nvx ~ oneover - 1)
  if(coef(m)[["oneover"]] < 0)
    stop("Trying to estimate trait variance from betas, and getting negative estimate.  Something is wrong.  You can 'fix' this by supplying an estimate of trait standard deviation yourself, as sdY=<value> in the dataset list.")
  return(sqrt(coef(m)[["oneover"]]))
}

#' variance of MLE of beta for quantitative trait, assuming var(y)=0
#'
#' Internal function
#' @title Var.data
#' @param f minor allele freq
#' @param N sample number
#' @return variance of MLE beta
#' @author Claudia Giambartolomei
#' @keywords internal
Var.data <- function(f, N) {
  1 / (2 * N * f * (1 - f))
}

##' variance of MLE of beta for case-control
##'
##' Internal function
##' @title Var.data
##' @inheritParams Var.data
##' @param s ???
##' @return variance of MLE beta
##' @author Claudia Giambartolomei
Var.data.cc <- function(f, N, s) {
  1 / (2 * N * f * (1 - f) * s * (1 - s))
}

#' Internal function, approx.bf.p
#'
#' Calculate approximate Bayes Factors
#' @title Internal function, approx.bf.p
#' @param p p value
#' @param f MAF
#' @param type "quant" or "cc"
#' @param N sample size
#' @param s proportion of samples that are cases, ignored if type=="quant"
#' @param suffix suffix to append to column names of returned data.frame
#' @return data.frame containing lABF and intermediate calculations
#' @author Claudia Giambartolomei
#' @keywords internal
approx.bf.p <- function(p,f,type, N, s, suffix=NULL, prior_var) {
  if(type=="quant") {
    V <- Var.data(f, N)
  } else {
    V <- Var.data.cc(f, N, s)
  }
  z <- qnorm(0.5 * p, lower.tail = FALSE)
  r <- prior_var / (prior_var + V)
  lABF = 0.5 * (log(1-r) + (r * z^2))
  ret <- data.frame(V,z,r,lABF)
  if(!is.null(suffix))
    colnames(ret) <- paste(colnames(ret), suffix, sep=".")
  return(ret)  
}

##' Internal function, approx.bf.estimates
##'
##' Calculate approximate Bayes Factors using supplied variance of the regression coefficients
##' @title Internal function, approx.bf.estimates
##' @param z normal deviate associated with regression coefficient and its variance
##' @param V its variance
##' @param sdY standard deviation of the trait. If not supplied, will be estimated.
##' @inheritParams approx.bf.p
##' @return data.frame containing lABF and intermediate calculations
##' @author Vincent Plagnol, Chris Wallace
approx.bf.estimates <- function (z, V, type, suffix=NULL, sdY=1, prior_var) {
  prior_var <- if (type == "quant") { prior_var*sdY^2} else { prior_var }
  r <- prior_var/(prior_var + V)
  lABF = 0.5 * (log(1 - r) + (r * z^2))
  ret <- data.frame(V, z, r, lABF)
  if(!is.null(suffix))
    colnames(ret) <- paste(colnames(ret), suffix, sep = ".")
  return(ret)
}

##' Internal function, process each dataset list for coloc.abf
##'
##' @title process.dataset
##' @param d list
##' @param suffix "df1" or "df2"
##' @return data.frame with log(abf) or log(bf)
##' @author Chris Wallace
process.dataset <- function(d, suffix=NULL, type=NULL, prior_var=0.15^2) { # if have an idea of sdY provide it in the data
  # message('Processing dataset')

  nd <- names(d)
  if (is.null(type)) stop('The variable type must be set, otherwise the Bayes factors cannot be computed')
  if ("BETA" %in% nd && "SE" %in% nd && !("MAF" %in% nd & "N" %in% nd || "sdY" %in% nd)) stop('Must give either MAF and N, or sdY')
  if ("BETA" %in% nd && "SE" %in% nd && ("MAF" %in% nd || "sdY" %in% nd)) {
  #   d$varbeta <- d$SE^2
  #   if (type=="quant" & !("sdY" %in% nd)) {
  #    d$sdY <- sdY.est(d$varbeta, d$MAF, d$N)
  #    message("Mean estimated sdY from data: ", mean(d$sdY))
  #    }
       df <- approx.bf.estimates(z=d$BETA/sqrt(d$varbeta),
                              V=d$varbeta, type=type, suffix=NULL, sdY=d$sdY, prior_var)
      return(df)
  }
  if ("PVAL" %in% nd & "MAF" %in% nd & "N" %in% nd) {
    if (type=="cc" & !("s" %in% nd))
      stop("Must specify s if type=='cc' and you want to use approximate Bayes Factors")
      df <- data.frame(PVAL = d$PVAL,
                     MAF = d$MAF,
                     snp=as.character(d$SNP))    
    # colnames(df)[-3] <- paste(colnames(df)[-3], suffix, sep=".")
    if (any(df$MAF==0 & df$PVAL==0)) {
       stop("There are some p-value equal to zero, remove these!")
       df <- subset(df, df$MAF>0 & df$PVAL>0) # all p values and MAF > 0
    }
    abf <- approx.bf.p(p=df$PVAL, f=df$MAF, type=type, N=d$N, s=d$s, suffix=NULL, prior_var)
    df <- cbind(df, abf)
    return(df)  
  }

  stop("Must give, as a minimum, either (beta, varbeta, type) or (pvalues, MAF, N, type)")
}

#' Internal function, logsum
#'
#' This function calculates the log of the sum of the exponentiated
#' logs taking out the max, i.e. insuring that the sum is not Inf
#' @title logsum
#' @param x numeric vector
#' @return max(x) + log(sum(exp(x - max(x))))
#' @author Claudia Giambartolomei, Vincent Plagnol
#' @keywords internal
logsum <- function(x) {
  my.max <- max(x)                              #take out the maximum value in log form
  my.res <- my.max + log(sum(exp(x - my.max ))) 
  return(my.res)
}

#' Internal function, logdiff
#'
#' This function calculates the log of the difference of the exponentiated
#' logs taking out the max, i.e. insuring that the difference is not negative
#' @title logdiff
#' @param x numeric
#' @param y numeric
#' @return max(x) + log(exp(x - max(x,y)) - exp(y-max(x,y)))
#' @author Chris Wallace
#' @keywords internal
logdiff <- function(x,y) {
  my.max <- max(x,y)                              #take out the maximum value in log form
  if (x>y) {
    my.res <- my.max + log(exp(x - my.max ) - exp(y-my.max))
  }
  if (x<y) {
    my.res <- my.max + log(exp(y-my.max) - exp(x - my.max ))
  }
  # If both x and y are zero, return zero?
  if (sum(x,y)==0) my.res = 0
  if (x==y) my.res = 0
  return(my.res)
}

#' Internal function, get_combos
#'
#' Find all possible combinations for the traits
#' @title get_combos
#' @param x names of each data frame
#' @return character vector with all combinations of names
#' @author Claudia Giambartolomei
#' @keywords internal
get_combos <- function(x) {
  x2 <- unlist(lapply(seq_along(x), function(m) 
    sapply(combn(x, m, simplify=FALSE), paste0, collapse='')))
  x3 <- expand.grid(rep(list(x2), length(x)))
  x4 <- sapply(apply(x3, 1, unique), function(x) paste0(sort(x), collapse=','))
  unique(grep('.*([^,]).*\\1', x4, val=TRUE, invert=TRUE))
}

#' Adjusted Bayes factors for each SNP and each single association and combinations of
#' sharing/not sharing of causal variants across the datasets in the input data
#'
#' @title adjust_bfs
#' @param listData A list of data frames to be analyzed. Each data frame
#'     must contain columns named SNP (the SNP name consistent across all the data frames), 
#'     BETA, 
#'     SE
#' @param prior_var One or more numbers for the prior variance of the ABF.
#' @param compute_sdY Estimate standard deviation of Y. If false assume sdY is 1
#' @return An array containing the log adjusted Bayes Factors for each SNP and
#'         for each SNP and each configuration combination
#' @examples
#' ABF <- adjust_bfs(listData, prior_var=c(0.01, 0.1, 0.5))
#' to view all configuration combinations for a SNP called "bp38878088"
#' ABF["bp38878088",,]
#' 
#' @keywords internal
#' @author Claudia Giambartolomei, Jimmy Liu
adjust_bfs <- function(listData, prior_var="default"){
  # keep only common SNPs in all data: 
  listData <- lapply(listData, function(x) x[x$SNP %in% Reduce(intersect, Map("[[", listData, "SNP")), ])
  if (nrow(listData[[1]])==0) stop("There are no common SNPs in the datasets: check that SNP names are consistent")       
  listData <- lapply(listData, function(df){ df[order(df$SNP),]})
  n_files <- length(listData)
  d <- letters[1:n_files]
  configs_cases <- do.call(expand.grid, lapply(d, function(x) c("", x)))[-1,]
  configs <- do.call(paste0, configs_cases)
  d_coloc <- configs[nchar(configs)>1]
  # coloc_cases <- configs_cases[apply(nchar(as.matrix(configs_cases)), 1, sum)>1,]
   
  # get adjusted ABF for each config combo
  snps <- as.character(listData[[1]]$SNP)
  ABF <- array(,dim=c(length(snps),1,length(configs)))
  dimnames(ABF) <- list(snps,NULL,configs) 
  # grid of priors?
  # grid_priors <- matrix(0, nrow(configs_cases), ncol(configs_cases))
  # This is coloc package using averages of priors
     single <- list()
     ABF_single = lapply (listData, function(x) {
                   adj_i_average <- data.frame(matrix(ncol=0,nrow=length(snps)),stringsAsFactors=FALSE)
                   if (length(prior_var)==1 & prior_var[1]=="default") {
                      if (x$type[1]=="quant") prior_var=0.15^2 else prior_var=0.2^2
                   }
                   message("Use prior variances of ", paste(prior_var, collapse=" "))

                   for (w_i in prior_var) {
                       adj_bf <- process.dataset(d=x, type=unique(x$type), prior_var=w_i)$lABF
                       adj_i_average<- cbind(adj_i_average, adj_bf)
                   }
                   single[[length(single) +1]] <- adj_i_average
                   #adj_bf <- apply(adj_i_average, 1, function(x) logsum(x) -log(length(prior_var)))
               })
    names(ABF_single) <- d


    for (i in 1:length(d)) {
          config <- d[i]
          ABF[,1,config] <- apply(ABF_single[[config]], 1, function(x) logsum(x) -log(length(prior_var)))
          }
          
    for (i in 1:length(d_coloc)) {
          config <- d_coloc[i]
          # ABF[,1,config] <- apply(ABF[,,unlist(strsplit(config, split=""))], 1, FUN=sum)
          internal_sum = Reduce('+', ABF_single[unlist(strsplit(config, split=""))])
          ABF[,1,config] <- apply(internal_sum, 1, function(x) logsum(x) -log(length(prior_var)))
          
          }
    
    return(ABF)
    }
    

#' Likelihood frame and posterior probability for combinations of sharing/not sharing
#' of causal variants across the datasets in the input data
#'
#' @title config_coloc
#' @param ABF is an array containing the single and coloc logBFs
#' @param n_files is one number, i.e. the number of traits to be analyzed
#' @param priors are numbers, the prior for one variant to be associated with 1 trait and with each additional trait
#' @return A data frame containing the likelihoods and posteriors for each configuration combination
#' @examples
#' lkl <- config_coloc(ABF, n_files=3, priors=c(1e-04, 1e-06, 1e-07))
#' 
#' @keywords internal
#' @author Jimmy Liu, Claudia Giambartolomei
config_coloc <- function(ABF, n_files, priors){
    # n_files = max(nchar(ABF))
    # have as many priors as traits
    d <- letters[1:n_files]
    configs_cases <- do.call(expand.grid, lapply(d, function(x) c("", x)))[-1,]
    configs <- do.call(paste0, configs_cases)
    final_configs <- get_combos(d) # includes the non-coloc configs

    # store sum of bfs for each configurarion, prior, and prior*sum(bfs_config)
    config_lkl <- array(,dim=c(length(final_configs),4,1), dimnames=list(final_configs,NULL,NULL))
    nsnps <- nrow(ABF)
    
    # likelihoods for colocalized configurations
    for (config in final_configs) {
        if (config %in% configs) {

        prior <- priors[nchar(config)]
        config_lkl[config,1,] <- prior
        #single_bfs[config,,] <- prior * sum(trait_bfs)
        config_lkl[config,2,] <- logsum(ABF[,,config])
        # Take out the - log(nsnps) to get per locus lkl
        config_lkl[config,4,] <- logsum(ABF[,,config]) -log(nsnps)
        }
        }
    # likelihoods for configurations that can be derived from colocalized configurations
     for (config in final_configs) {
        if (!config %in% configs) {
         # print(config)
         # lH3.abf <-  logdiff(lH1.abf + lH2.abf, lH4.abf)
            # left side
            left_trait_bfs <- 0
            prior_num <- 1.0
            composite_bfs <- unlist(strsplit(config, split=","))
            for (i in composite_bfs) { # gather info on relevant single_bfs
               # print(i)
               left_trait_bfs <- left_trait_bfs + config_lkl[i,2,]
               prior_num <- prior_num * config_lkl[i,1,]
               }
            # relevant coloc configuration
            coloc_config <- gsub(",", "", config)
            # must sort otherwise it doesn't find the match
            coloc_config <- as.character(lapply(lapply(strsplit(coloc_config,NULL),sort),paste,collapse=""))
            # right_trait_bfs <- ( prior_num / single_bfs[coloc_config,1,] ) * single_bfs[coloc_config,3,]
            # right_trait_bfs <- prior_num * single_bfs[coloc_config,3,]
            right_trait_bfs <- config_lkl[coloc_config,2,]
            # print(left_trait_bfs)
            # print(right_trait_bfs)
            config_bf <- logdiff(left_trait_bfs, right_trait_bfs)
            config_lkl[config,2,] <- config_bf
            # Re-check this when using >2 traits: log(nsnps)^length(composite_bfs)
            config_lkl[config,4,] <- config_bf - log(nsnps^2)
            config_lkl[config,1,] <- prior_num
        }
        }
      config_lkl[,3,] <- log(config_lkl[,1,]) + config_lkl[,2,] # prior * Sum(BFs)

    # add config where nothing is associated
    config_ppas <-  as.data.frame(config_lkl)
    names(config_ppas)=c("prior", "sumbf", "loglkl", "logBF_locus")
    # add zero
    config_ppas <- rbind(config_ppas, zero=c(0.999, 0, 0, 0))

    all.abf <- config_ppas$loglkl
    my.denom.log.abf <- logsum(all.abf)
    config_ppas$PPA <- exp(all.abf - my.denom.log.abf)

    # print(config_ppas)
    return(list(config_ppas, nsnps))
}


#' Posterior probability that each SNP is THE causal variant for a shared signal
#'
#' @title snp_ppa
#' @param ABF is an array containing the single and coloc logBFs
#' @param n_files is one number, i.e. the number of traits to be analyzed
#' @param priors is the prior for one variant to be associated with 1 trait and with each additional trait
#' @return A data frame containing the likelihoods and posteriors for each configuration combination
#' @examples
#' snp <- snp_ppa(ABF, n_files=3, config_ppas=lkl[[1]])
#' 
#' @keywords internal
#' @author Jimmy Liu, Claudia Giambartolomei
snp_ppa <- function(ABF, n_files, config_ppas, save.SNP.info){
    d <- letters[1:n_files]
    configs_cases <- do.call(expand.grid, lapply(d, function(x) c("", x)))[-1,]
    configs <- do.call(paste0, configs_cases)
    final_configs <- get_combos(d) # includes the non-coloc configs
    SNP.PP.H4.df <- as.data.frame(sapply(as.data.frame(ABF), function(x) exp(x-logsum(x)) ))
    names(SNP.PP.H4.df) = dimnames(ABF)[[3]]
    rownames(SNP.PP.H4.df) = dimnames(ABF)[[1]]

    #best.t <- exp(rowSums(as.data.frame(ABF)[,grep(i, colnames(as.data.frame(ABF)))]) - my.denom.log.abf)

    # printout
    # print out colocalization posteriors for each trait
    coloc_ppas = numeric()
    best.snp.coloc = character()
    # for (i in d) {
    for (i in configs) {
        final_config_ls <- strsplit(final_configs, ",")
        trait_coloc_list <- Filter(function(x) any(nchar(grep(i, x, value = TRUE))>1), final_config_ls)
        trait_coloc <- sapply(trait_coloc_list,function(x) paste(x, collapse=","))
        trait_coloc <- c(trait_coloc, grep("a[a-z]c", configs, perl=T, value=T))
        trait_coloc <- trait_coloc[!duplicated(trait_coloc)]
        # print(trait_coloc)
        hh4 <- sum(config_ppas[rownames(config_ppas) %in% trait_coloc,"PPA"])

        # what is the denominator here, other coloc models or all the coloc and non coloc?
        if (sum(names(SNP.PP.H4.df) %in% trait_coloc)>1) {
            SNP.coloc <- rowSums(SNP.PP.H4.df[,names(SNP.PP.H4.df) %in% trait_coloc]) 
            } else {
            SNP.coloc <- SNP.PP.H4.df[,names(SNP.PP.H4.df) %in% trait_coloc]
            }
        names(SNP.coloc) <- rownames(SNP.PP.H4.df)
        s = names(which.max(SNP.coloc))
        message("Best SNP per trait ", i, ": ", s)
        best.snp.coloc <- c(best.snp.coloc, s)
        #if (names(which.max(SNP.pp4)) != causal.1) browser("***causal not recognized")
        coloc_ppas <- c(coloc_ppas, hh4)
        message("Probability that trait ", i, " colocalizes with at least one other trait = ", signif(hh4, digits = 2))
     }
     # SNP.pp4.any <- rowSums(ABF[,,dimnames(ABF)[[3]] %in% any_coloc])/denom
     # message("Best SNP for any trait colocalizing: ", names(which.max(SNP.pp4.any)))
     
     # names(coloc_ppas) <- d
     best_snp = cbind.data.frame(coloc_ppas, best.snp.coloc)
     rownames(best_snp) <- configs
     if (!save.SNP.info) {
         return(best_snp)
     }
     if (save.SNP.info) {
         return(list(best_snp, SNP.PP.H4.df))
     }
}


#' Bayesian multiple trait colocalization analysis using list of data.frames
#'
#' Runs \code{\link{adjust_bfs}}, \code{\link{config_coloc}}, \code{\link{snp_ppa}}
#' 
#' @title moloc_test
#' @param listData A list of data frames to be analyzed. Each data frame
#'     must contain columns named SNP (the SNP name consistent across all the data frames), 
#'     BETA, 
#'     SE
#' @param prior_var One or more numbers for the prior variance of the ABF.
#' @param priors are numbers, the prior for one variant to be associated with 1 trait and with each additional trait
#' @param compute_sdY Estimate standard deviation of Y. If false assume sdY is 1
#' 
#' @param ... parameters passed to \code{\link{adjust_bfs}}, \code{\link{config_coloc}}, \code{\link{snp_ppa}}
#' @return A list of three elements: 
#'         First is a data frame with 4 variables: priors ('prior'), likelihoods ('sumbf') and 
#'           posteriors ('loglkl' and 'PPA') for each configuration; 
#'         Second is a number, the number of SNPs in common in the region;
#'         Third is a data frame with 2 variables: 
#'           SNP with the best posterior for each scenario ('coloc_ppas'), 
#'           SNP name with the best posterior for each scenario ('best.snp.coloc').
#' @export
#' @author Claudia Giambartolomei
#' @examples
#' moloc <- moloc_test(listData) # uses default priors
#' 
#' @author Claudia Giambartolomei
moloc_test <- function(listData, prior_var=c(0.01, 0.1, 0.5), priors=c(1e-04, 1e-06, 1e-07), save.SNP.info = FALSE) {
    n_files <- length(listData)
    if(missing(priors)) {
      priors <- 10^-(seq(from=4, to=4+n_files-1, by=1))
      message("Use default priors: ", paste(priors, collapse=","))
    } else {
      if (length(priors)!=n_files) stop("Priors need to be the same length as the number of traits")
      priors <- as.numeric(priors)
    }
    
    
    for (i in 1:length(listData)) {
       nd <- colnames(listData[[i]])
       # Add type
       listData[[i]]$type <- ifelse("Ncases" %in% names(listData[[i]]), "cc", "quant")
       # Add s
       # listData[[i]]$s = ifelse("Ncases" %in% names(listData[[i]]), as.numeric(listData[[i]]$Ncases/listData[[i]]$N), 0.5)
       listData[[i]]$s = 0
       if ("Ncases" %in% names(listData[[i]])) s=as.numeric(listData[[i]]$Ncases/listData[[i]]$N)
       if (!("Ncases" %in% names(listData[[i]]))) s=0.5
       listData[[i]]$s = s
       if ("BETA" %in% nd && "SE" %in% nd && !("MAF" %in% nd & "N" %in% nd || "sdY" %in% nd)) stop('Must give either MAF and N, or sdY')
       if ("BETA" %in% nd && "SE" %in% nd && ("MAF" %in% nd || "sdY" %in% nd)) {
          listData[[i]]$varbeta <- listData[[i]]$SE^2
          if (unique(listData[[i]]$type)=="quant" & !("sdY" %in% nd)) {
           listData[[i]]$sdY <- sdY.est(listData[[i]]$varbeta, listData[[i]]$MAF, listData[[i]]$N)
           message("Mean estimated sdY from data", i, ": ", mean(listData[[i]]$sdY))
          }
       }
    }

    # Add MAF if missing from the first available matching data with maf
    # for (i in 1:length(listData)) {
    #   if (length(grep("MAF", names(listData[[i]])))>0) 
    #   MAF <- listData[[i]]$MAF
    #   names(MAF) <- listData[[i]]$SNP
    #   print(i)
    #   break
    #}
      
    for (i in 1:length(listData)) {
       if (length(grep("MAF", names(listData[[i]])))==0) {
       stop("Dataset ", i, " does not contain MAF; use matching MAF from other datasets.") 
    #   id <- match(names(MAF), listData[[i]]$SNP)
    #   listData[[i]]$MAF = MAF[id]
       }
    }
    ABF <- adjust_bfs(listData=listData, prior_var=prior_var)
    lkl <- config_coloc(ABF, n_files, priors)
    snp <- snp_ppa(ABF, n_files=n_files, config_ppas=lkl[[1]], save.SNP.info)
    nsnps <- lkl[[2]]
    lkl <- lkl[[1]]
    # Report per locus likelihoods (adjusted for nsnps)
    # lkl$sumbf <- lkl$sumbf - log(nsnps); for the no coloc it's -2* log(nsnps)
    lkl <- lkl[,c("prior", "sumbf", "logBF_locus", "PPA")]
    if (!save.SNP.info) {
        res <- list(lkl, nsnps, snp)
    }
    if (save.SNP.info) {
       snp_info = merge(data.frame(listData), data.frame(lABF=ABF), by.x="SNP", by.y="row.names")
       SNP.PP.H4.df = snp[[2]]
       names(SNP.PP.H4.df) = paste("SNP.PP.", names(SNP.PP.H4.df), sep="")
       SNP.PP.H4.df$SNP = rownames(SNP.PP.H4.df)
       snp_info = merge(snp_info, SNP.PP.H4.df, by="SNP")
       res <- list(lkl, nsnps, snp, snp_info)
       }
    names(res) <- c("priors_lkl_ppa", "nsnps", "best_snp")
    return(res)
}


######
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
