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


#' Internal function, get_causalM
#'
#' Get causal matrix from binary matrix; v is the probability per SNP
#' @title get_combos
#' @param x matrix
#' @return matrix with causal configurations
#' @author Claudia Giambartolomei
#' @keywords internal
get_causalM <- function(mm, v) {
	H_causal <- matrix(0, nrow(mm), ncol(mm))

	# Build causal frame
	for (j in 1:ncol(mm)) {
	 for (i in 1:nrow(mm)) {
		H_causal[i,j] = ifelse(mm[i,j]> 0, v[j], 1-v[j])
	 }
	}
	return(H_causal)
}

#' Internal function, nested_lapply
#'
#' @title nested_lapply
#' @param
#' @return
#' @author Claudia Giambartolomei
#' @keywords internal
nested_lapply <- function(data, fun) {
    lapply(data, function(sublist) { lapply(sublist, fun) })
}

#' Internal function, nested_lapply_args
#'
#' @title nested_lapply_args
#' @param
#' @return
#' @author Claudia Giambartolomei
#' @keywords internal
nested_lapply_args <- function(data, fun, v) {
    lapply(data, function(sublist) { lapply(sublist, fun, v) })
}

#' Internal function, logsum
#'
#' @title logsum
#' @param
#' @return
#' @author Claudia Giambartolomei
#' @keywords internal
logsum <- function(x) {
  my.max <- max(x) ##take out the maximum value in log form
  my.res <- my.max + log(sum(exp(x - my.max ))) 
  return(my.res)
}

#' priors from causal probabilities
#'
#' @title get_priors
#' @param v, numeric vector of probabilities for each SNP included in the analysis
#' @param n_files, number of traits analyzed, for now can be either 2 or 3
#' @param priors, numeric vectors with three numbers;
#' @return probabilities for each hypotheses in moloc (5 for 2 traits; 15 for three traits).
#' @examples
#' v = runif(1:1000)
#' n_files = 3
#' get_priors(v, n_files)
#'
#' @export
#' @author Claudia Giambartolomei
get_priors = function(v, n_files) {

############################################
# Build configurations and group into lists
############################################
	nsnps <- length(v) # number of snps
	d <- letters[1:n_files]
	final_configs <- get_combos(d) 

	mylist = list()
	print("Constructing single trait configurations...")
	for (i in 1:nsnps) {
    	configuration<- matrix(0, nrow = n_files, ncol = nsnps)
    	configuration[1, i] = 1
    	mylist[[length(mylist)+1]] = configuration
 	}

        print("Constructing two traits configurations...")
	for (i in 1:nsnps) {
 		for (j in 1:nsnps) {
    		configuration<- matrix(0, nrow = n_files, ncol = nsnps)
    		configuration[1, i] = 1
    		configuration[2, j] = 1
    		mylist[[length(mylist)+1]] = configuration
 		}
	}

        print("Constructing three traits configurations...")
	for (i in 1:nsnps) {
 		for (j in 1:nsnps) {
  			for (k in 1:nsnps) {
   				configuration<- matrix(0, nrow = n_files, ncol = nsnps)
    			message("i ", i, "; j ", j, " ;k ", k) 

    			configuration[1, i] = 1
    			configuration[2, j] = 1
    			configuration[3, k] = 1
    			mylist[[length(mylist)+1]] = configuration
 			}
 		}
	}


	if (n_files==3) {
		three = data.frame(config=final_configs, numb_traits_sharing = c(1, 1, 1, 2, 1, 1, 1, 2, 1, 2, 2, 2, 2, 3), in_a = c(1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1), in_b = c(0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1), in_c = c(0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1), extra=c(NA,NA,NA,"0,1,1", NA, NA, NA, "1,0,1", NA, "1,1,0", NA, NA, NA, NA))
		# don't need all the configurations
		three = three[three$config %in% c("a", "ab", "a,b", "a,bc", "abc", "a,b,c"),]

		cond = three
		groupedll = list()
		for (i in 1:nrow(cond)) {
			subsetll = mylist[ unlist( lapply(mylist, FUN=function(x) rowSums(x)==cond$in_a[i] && rowSums(x)[2]==cond$in_b[i] && rowSums(x)[3]==cond$in_c[i] && max(colSums(x))==cond$numb_traits_sharing[i]) ) ]
			if (!is.na(cond$extra[i]) ) {
				z = as.numeric(unlist(strsplit(as.character(cond$extra[i]),",")))
				subsetll2 = subsetll [ unlist( lapply(subsetll, FUN=function(x)  any(apply(x, 2, function(x) identical(as.numeric(x),z)))))]
				subsetll = subsetll2
			}
			groupedll[[length(groupedll)+1]] = subsetll
		}
		names(groupedll) = cond$config
	}

	if (n_files==2) {
		two = data.frame(config=final_configs, numb_traits_sharing = c(1, 1, 1, 2), in_a = c(1, 1, 0, 1), in_b = c(0, 1, 1, 1), extra=c(NA,NA,NA,NA)) 

		cond = two
		groupedll = list()
		for (i in 1:nrow(cond)) {
			subsetll = mylist[ unlist( lapply(mylist, FUN=function(x) rowSums(x)==cond$in_a[i] && rowSums(x)[2]==cond$in_b[i] && max(colSums(x))==cond$numb_traits_sharing[i]) ) ]
			groupedll[[length(groupedll)+1]] = subsetll
		}
		names(groupedll) = cond$config
	}

	# names(groupedll) <- c("a", "ab", "a.b.c", "abc", "a.b", "ab.c")
    # a*3, ab*3, ab.c*3, a.b*3, not repeated: a.b.c, abc
	############################################
	# Fill in probabilities 
	############################################

	# ll is a *named* list of data frames with causal configurations for each hypotheses
	causal_M <- nested_lapply_args(groupedll, get_causalM, v=v)
	# fill in: names(causal_M)[as.logical(lapply(causal_M, function(x) length(x)==0))]
	causal_M$b=causal_M$a; causal_M$c=causal_M$a
	causal_M$ac=causal_M$ab; causal_M$bc=causal_M$ab
	causal_M$'a,c'=causal_M$'a,b'; causal_M$'b,c'=causal_M$'a,b'
        causal_M$'ac,b'=causal_M$'a,bc'; causal_M$'ab,c'=causal_M$'a,bc' # this is possible only if priors for a, b, c are the same!
	if (length(causal_M)!=14) stop("Not all hypotheses are computed")
	# Take the product of the SNPs and configuration
	causal_P <- nested_lapply(causal_M, prod)
	# Sum across configurations belonging to the same hypotheses
	h = c()
	for (i in 1:length(causal_P)) {
		c = Reduce("+", causal_P[[i]])
		h = c(h,c)
	}
	# names(h) = names(causal_M)[as.logical(lapply(causal_M, function(x) length(x)!=0))]
	names(h) = paste("pp_", names(causal_M), sep="")

	############################################
  	pp = h

  	print(signif(pp,3))
  	return(pp)
}
