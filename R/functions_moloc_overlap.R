#' Internal function, stats_p
#'
#' Get beta, SE from p
#' @title Internal function, stats_p
#' @param p p value
#' @param f MAF
#' @param type "quant" or "cc"
#' @param N sample size
#' @param s proportion of samples that are cases, ignored if type=="quant"
#' @param suffix suffix to append to column names of returned data.frame
#' @return data.frame containing beta, SE
#' @author Claudia Giambartolomei
#' @keywords internal
stats_p <- function(p, n, maf) {
    z <- qnorm(0.5 * p, lower.tail = FALSE)
    var_mle <- 1/(2*maf*(1-maf) * ( n + z^2))
    # var_mle <- 1 / (2 * maf * (1 - maf) * n)
    SE <- sqrt(var_mle)
    BETA = z * SE
    # vars = 2 * maf * ( 1 - maf) * n * SE^2 * (n -1) + 2 * maf * ( 1- maf) * n * BETA^2 
    # vars = sqrt(median(vars/(n-1)))
    # Neff_est <- vars^2/(2*maf*(1-maf)*SE^2) - (BETA^2/SE^2) +1
    df <- cbind.data.frame(BETA, SE)
    return(df)
}
 
#' Internal function, stats_p_cc
#'
#' Internal function
#' @title stats_p_cc
#' @inheritParams Var.data
#' @param s proportion of samples that are cases
#' @param f minor allele freq
#' @param N sample number
#' @return variance of MLE beta
#' @examples
#' stats = stats_p_cc(listData[[1]]$PVAL, listData[[1]]$N, listData[[1]]$F, listData[[1]]$Ncases/listData[[1]]$N)
#' @author Claudia Giambartolomei
#' @keywords internal
stats_p_cc <- function(p, n, maf, s) {
    z <- qnorm(0.5 * p, lower.tail = FALSE)
    var_mle <- (s*(1-s))/(2*maf*(1-maf) * ( n + z^2))
    SE <- sqrt(var_mle)
    BETA = z * SE
    df <- cbind.data.frame(BETA, SE)
    return(df)
}

#' Internal function, stats_p_cc
#'
#' Internal function
#' @title stats_p_cc
#' @inheritParams stats_p
#' @inheritParams stats_p_cc
#' @param s proportion of samples that are cases
#' @param f minor allele freq
#' @param N sample number
#' @return variance of MLE beta
#' @keywords internal
get_stats <- function(x) {
      message("Var and Beta from pvalues")
      if (!all(c("MAF", "N") %in% names(x))) stop("Must provide freq and N")
      if (all(c("N", "Ncases") %in% names(x))) {
         message("This is a case-control")
         x <- cbind.data.frame(x, stats_p_cc(x$PVAL, x$N, x$MAF, x$Ncases/x$N))
      } else {
         message("This is not a case-control")
         x <- cbind.data.frame(x, stats_p(x$PVAL, x$N, x$MAF))
      }
}

#' Internal function, is_pos_def
#'
#' This function tests if a matrix is positive definite
#' @title is_pos_def
#' @param A matrix
#' @author Claudia Giambartolomei, Jimmy Liu
#' @keywords internal
is_pos_def <- function(A) {
    cholStatus <- try(u <- chol(A), silent = TRUE)
    cholError <- ifelse(class(cholStatus) == "try-error", FALSE, TRUE)
    return(cholError)
    }


#' Internal function, get_psd
#'
#' This function attempts to make a matrix positive definite
#' @title get_psd
#' @param A matrix
#' @author Claudia Giambartolomei, Jimmy Liu
#' @keywords internal
get_psd <- function(A) {
    d = diag(A)
    # detA = np.linalg.det(A)
    while (!is_pos_def(A)) {
        # while detA <= 0:
        A = 0.999 * A
        for (i in 1:length(A[1,])){
            A[i, i] = d[i]
            # detA = np.linalg.det(A)
            }
        }
    return(A)
    }


#' Internal function, make_sigma
#'
#' This function creates Wakefield's ABF for each SNP, taking account of the correlation between traits
#' @title make_sigma
#' @param cor_mat matrix, correlation between the traits for each SNP
#'     If assuming no overlap, cor_mat is identity matrix.
#' @return matrix of ABF adjusted for correlations
#' @author Jimmy Liu
#' @keywords internal
make_sigma <- function(cor_mat, v, w){
    sigma = diag(v+w)
    for(i in 1:(length(v)-1)) {
        for(j in (i+1):length(v)) {
            c = cor_mat[i, j]
            sigma[i, j] = c * sqrt((v[i] + w[i]) * (v[j] + w[j]))
            sigma[j, i] = c * sqrt((v[i] + w[i]) * (v[j] + w[j]))
        }
    }
    return(sigma)
}

#' Adjusted Bayes factors for each SNP and each single association and combinations of
#' sharing/not sharing of causal variants across the datasets in the input data
#'
#' @title adjust_bfs
#' @param listData A list of data frames to be analyzed. Each data frame
#'     must contain columns named SNP (the SNP name consistent across all the data frames), 
#'     BETA, 
#'     SE
#' @param overlap Logical, do the individuals in the datasets overlap?
#' @param prior_var One or more numbers for the prior variance of the ABF.
#' @param compute_sdY Estimate standard deviation of Y. If false assume sdY is 1
#' @return An array containing the log adjusted Bayes Factors for each SNP and
#'         for each SNP and each configuration combination
#' @examples
#' ABF <- adjust_bfs_overlap(listData, overlap=FALSE, prior_var=c(0.01, 0.1, 0.5))
#' ABF <- adjust_bfs(listData=listData, overlap=overlap, prior_var=prior_var, compute_sdY=compute_sdY, from_p=from_p)
#' to view all configuration combinations for a SNP called "bp38878088"
#' ABF["bp38878088",,]
#' 
#' @keywords internal
#' @author Jimmy Liu, Claudia Giambartolomei
adjust_bfs_overlap <- function(listData, prior_var=0.15^2, overlap = FALSE, compute_sdY = FALSE, from_p=FALSE){
  # keep only common SNPs in all data: 
  listData <- lapply(listData, function(x) x[x$SNP %in% Reduce(intersect, Map("[[", listData, "SNP")), ])
  if (nrow(listData[[1]])==0) stop("There are no common SNPs in the datasets: check that SNP names are consistent")       
  listData <- lapply(listData, function(df){ df[order(df$SNP),]})
  ##### This is not right
  if (from_p) {
    listData <- lapply(listData, function(x) {
         y = get_stats(x)
         return(y)
         })    # if (!all(c("BETA", "SE") %in% names(x)))
  }       
  
  #######
  n_files <- length(listData)
  d <- letters[1:n_files]
  configs_cases <- do.call(expand.grid, lapply(d, function(x) c("", x)))[-1,]
  configs <- do.call(paste0, configs_cases)
  d_coloc <- configs[nchar(configs)>1]
  # coloc_cases <- configs_cases[apply(nchar(as.matrix(configs_cases)), 1, sum)>1,]
   
    # adjusted bfs
    # correlation matrix
    # if assuming no overlap, cor_mat is identity matrix
    cor_mat <- diag(1, n_files)
    if (overlap) { 
        # for each combination of traits, i.e. 
        #loop through unique combination
        for(i in 1:(n_files -1)) {
            for(j in (i+1):n_files) {
                beta1 = listData[[i]]$BETA
                beta2 = listData[[j]]$BETA
                cortest = cor.test(beta1,beta2) # pvalue extracted from cor.test is not exact...
                if (cortest$p.value >= 0.01) { 
                    cor_mat[i, j] = 0.0
                    cor_mat[j, i] = 0.0
                } else {
                    cor_mat[i, j] = cortest$estimate
                    cor_mat[j, i] = cortest$estimate
                } 
        }
       }
    }
    # configs - different adjusted ABF depending on configs
    # get adjusted ABF for each config combo
    snps <- as.character(listData[[1]]$SNP)
    ABF <- array(,dim=c(length(snps),1,length(configs)))
    dimnames(ABF) <- list(snps,NULL,configs) 
    # grid of priors?
    # grid_priors <- matrix(0, nrow(configs_cases), ncol(configs_cases))
    if (compute_sdY) {
        if (!all(unlist(lapply (listData, function(x) c("MAF") %in% names(x) ) ))) stop("Must provide a column MAF")
        sdY = mapply(sdY.est, lapply(lapply(listData, "[[", "SE"), '^',2), lapply(listData, "[[", "MAF"),  lapply(listData, "[[", "N"), SIMPLIFY = FALSE)
        varY = lapply(sdY, '^',2)
        names(varY) = d
        # varY_configs <- apply(configs_cases, 1, function(x) prod(unlist(varY[unlist(x)])) )
    }

    # create data frames of betas and var across all traits
    names(listData) <- d
    var = lapply(listData, "[[", "SE")
    var = lapply(var, '^',2)
    var = do.call(cbind, var)

    betas = lapply(listData, "[[", "BETA")
    betas = do.call(cbind, betas)

    # data_sets <- apply(configs_cases,1, function(x) do.call(cbind.data.frame, listData[unlist(x[x!=''])]))
    means <- rep(0, n_files)
    w_0 <- rep(0, n_files)
    sigma_h0 <- lapply(split(var, seq(NROW(var))), FUN=make_sigma, cor_mat=cor_mat, w=w_0)

    pdf_null <- c()
    for (j in 1:length(sigma_h0)) {
         x <- dmvnorm(betas[j,], means, sigma_h0[[j]])
         pdf_null <- c(pdf_null, x)
    }

    for (i in 1:length(d)) {
               config = d[i]
               adj_i_average <- data.frame(matrix(ncol=0,nrow=length(snps)),stringsAsFactors=FALSE)
               for (w_i in prior_var) {
                  # This needs to be in the correct order for each trait
                  w <- ifelse(d %in% as.matrix(config), w_i, 0)
                  # w <- as.numeric(ifelse(nchar(as.matrix(configs_cases[i,])), w_i, 0))
                  if (compute_sdY) {
                      w <- w * varY[[config]]
                  }

                  sigma_h1 <- lapply(split(var, seq(NROW(var))), FUN=make_sigma, cor_mat=cor_mat, w=w)
                  pdf_alt <- c()
                  for (j in 1:length(sigma_h0)) {
                      x <- dmvnorm(betas[j,], means, sigma_h1[[j]])
                      pdf_alt <- c(pdf_alt, x)
                  }
                  adj_bf <- pdf_alt / pdf_null

                  adj_bf <- mapply(function(x, y)
                          if (y ==0 | y<1e-300) res <- 1.0 else res <- x / y, 
                          pdf_alt, pdf_null)

                  adj_i_average<- cbind(adj_i_average, adj_bf)
               }  
               adj_bf <- apply(adj_i_average, 1, mean)

               if (any(log(adj_bf)== "Inf") | any(log(adj_bf)== "-Inf") ) stop("Division impossible for adj_bf")

          ABF[,1,config] <- log(adj_bf)
          }
          
    for (i in 1:length(d_coloc)) {
               config <- d_coloc[i]
          ABF[,1,config] <- apply(ABF[,,unlist(strsplit(config, split=""))], 1, FUN=sum)
          }
    return(ABF)
    }