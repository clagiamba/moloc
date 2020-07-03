#' Genome-wide co-localization with three traits
#'
#' @title coloc.genome
#' @param listData, a list of data frames in this order: gwas, and other QTLs
#' The data frames need to have columns: "SNP" or "CHR", "POS"; "BETA", "SE"; "N", "MAF" (to estimate sdY)
#' 										  if a case control: "Ncases"
#'                                        If eqtl/mqtl: "ProbeID"
#'                                        Optionally, "A1", "A2" if want to match alleles
#'                                        
#' @param bed, a data frame with "CHR", "START", "STOP", and 
#'             "ProbeID" matching either eqtl or mqtl
#' @param outfolder, path of folder where the results will be written to
#' @param prefix, character, name to give to the output file
#' @param save.SNP.info, logical, should info of each SNP be saved? 
#'        If this option is true, beware that it will takes up a lot of time and space.
#' @param cores, number. See foreach package.
#' @param have_alleles, logical. If TRUE, matches alleles to be the same as gwas data.
#' @param bychrpos, logical. If TRUE, uses "CHR", "POS" columns 
#'        as SNP id to match across data frames;
#'        If TRUE, make sure not to have another column called "SNP".
#' @param prior_var, numeric vector specifying any number of values;
#'        These will be used to average the ABF computation across different prior variances.
#' @param priors, numeric vectors with three numbers; 
#'        First number is the prior of associaion of a SNP with any of the traits;
#'        Second number is the prior of a SNP being associated with two traits;
#'        Third number is the prior of a SNP being associated with three traits.
#' @param min_nsnps, number. The minimum number of matching SNPs across 
#'        the three data frames to consider. 
#' @return results in outfolder: text file with colocalization results for each locus; 
#'         text file with removed snps; 
#'         if save.SNP.info = TRUE, a folder called "coloc.output.perSNP"
#'         will save a file per locus with SNP information such as ABF (one row per matching SNP)
#' @examples
#' library(mvtnorm)
#' library(foreach)
#' library(doParallel)
#' library(data.table)
#' 
#' Create bed file with combination of ProbeIDs
#' library(GenomicRanges)
#' DT <- as.data.table(data_genome[[3]])
#' methyl_table <- DT[, list(CHR= unique(CHR), START = min(POS), STOP = max(POS)), by = ProbeID]
#' bed1.gr <- GRanges(seqnames = bed$CHR,IRanges(start = bed$START, end= bed$STOP))
#' bed2.gr <- GRanges(seqnames = methyl_table$CHR,IRanges(start =methyl_table$START, end= methyl_table$STOP))
#' my.overlap <- findOverlaps(query = bed2.gr, subject = bed1.gr)
#' bed = cbind(bed[my.overlap@subjectHits,], methyl_table[my.overlap@queryHits,])
#' bed = bed[,c(1,6,7,8,9,10)]
#' Setting the prior_var to "default" uses coloc defaults
#' coloc.genome(listData = data_genome, bed, cores=1, have_alleles=TRUE, bychrpos=TRUE, prior_var=NULL, priors=c(1e-04, 1e-06, 1e-07), min_nsnps = 50, write=FALSE, outfolder = "test")
#' 
#' @export
#' @author Claudia Giambartolomei
coloc.genome <- function(listData, bed, prefix = "pref", save.SNP.info=FALSE, cores=20, have_alleles=TRUE, bychrpos=TRUE, prior_var="default", priors=c(1e-04, 1e-06, 1e-07), min_nsnps = 50, takelog = FALSE, write=TRUE, outfolder = "test", forcePVAL = FALSE, compare_coloc = TRUE){
# if there is a BETA in the data, coloc will use the BETA and SE, to force use of PVAL instead use forcePVAL = TRUE

  ########
  merge_results  <- function(a, b) {
    if(is.null(a) & is.null(b)) {
        return(NULL)
    } else if (is.null(a)){
        return(b)
    } else if (is.null(b)){
        return(a)
    } else {
        return(rbind(a,b))
    }
  }
  ########
  # You need the suggested package for this function    
  if (!requireNamespace("foreach", quietly = TRUE)) {
    stop("Pkg foreach needed for this function to work. Please install it.",
      call. = FALSE)
  }
  if (!requireNamespace("doParallel", quietly = TRUE)) {
    stop("Pkg doParallel needed for this function to work. Please install it.",
      call. = FALSE)
  }

  ########## 
  outfname = paste(outfolder, prefix, '_summary.tab', sep='')
  if (!file.exists(outfolder)) dir.create(outfolder)
  
  if (bychrpos) {
    # If there is a column called SNP rename it first
    for (i in 1:length(listData)) {
       s = grep("^SNP$", names(listData[[i]]))
       if (length(s)==1) {
       names(listData[[i]])[s] <- "SNP.input"
        }
       }
    chrpos = lapply(listData, function(x) x$chrpos = paste(x$CHR, x$POS, sep=":"))
    listData = Map(cbind, listData, SNP = chrpos)
    listData <- lapply(listData, function(df) {df[c("SNP")] <- lapply(df[c("SNP")], as.character); df})
  }
  
  for (i in 1:length(listData)) {
    haveBETA = "BETA" %in% names(listData[[i]])
    if (haveBETA &!forcePVAL) {
        if (length(listData[[i]]$SE[is.na(listData[[i]]$SE)])>0) {
            warning("There are missing SE in the data!")
            listData[[i]] = subset(listData[[i]], !is.na(SE))
        }
       if (length(listData[[i]]$BETA[listData[[i]]$BETA<0])>0) log=TRUE  else log=FALSE # if there are negative value, then it could be already logOR?
       if (!log) warning("Dataset ",  i, " looks like a case-control: log the betas! The betas will be logged if takelog is set to TRUE, otherwise log them beforehand")
       # before taking the log must remove the SNPs with beta = 0
       if (!log & takelog) (listData[[i]] = subset(listData[[i]], BETA != 0))
       if (!log & takelog) (listData[[i]]$BETA = log(listData[[i]]$BETA))
    }
  }
  
    # beta = any(unlist(lapply(listData, function(df) {length(df$BETA[df$BETA==0])>0})))  
    # pval = any(unlist(lapply(listData, function(df) {length(df$BETA[df$PVAL==1])>0})))
    # if (any(beta & pval)) {
    #   stop("There are some p-value equal to zero, remove these!")
    #}
     
  if (!all(c("CHR", "START", "STOP") %in% names(bed))) stop("Bed file is missing info")
  bed$CHR <- gsub("chr", "", bed$CHR)
  bed$CHR <- as.numeric(bed$CHR)
  
  for (i in 1:length(listData)) {
       nm = names(listData[[i]])
       cols = c("SNP", "CHR", "POS", "N", "MAF")
       if ("PVAL" %in% nm) cols = c("PVAL", cols)
       if (all(c("BETA", "SE") %in% nm) & !forcePVAL) cols = c("BETA", "SE", cols)       
       if (!(all(cols %in% nm) && (all(c("BETA", "SE") %in% nm) | "PVAL" %in% nm))) stop("One or more essential columns are missing")
       if ("ProbeID" %in% nm) cols = c("ProbeID", cols)
       if (all(c("CHR", "POS") %in% nm)) cols = c("CHR", "POS", cols)
       if ("sdY" %in% nm) cols = c("sdY", cols)
       if ("Ncases" %in% nm) cols = c("Ncases", cols)
       if (all(c("A1", "A2") %in% nm)) cols = c("A1", "A2", cols)
       listData[[i]] = listData[[i]][,cols]
       }
  
  ##################################################### now start the loop
  # Now go over all regions that overlap between eQTL table and input.data
  message("Running in parallel")
  registerDoParallel(cores=cores)

  message("Looping through ", nrow(bed), " regions")

  res.all <- data.frame()
  removed_snp_list = data.frame()

  res.all  <-  foreach(i=1:nrow(bed), .combine=merge_results) %dopar% {
  #for (i in 1:nrow(bed)) {
  res <- data.frame()
       data <- listData
       probeids = names(bed)[grep("ProbeID", names(bed))]
       pos.start = bed$START[i]
       pos.end = bed$STOP[i]
       chrom = bed$CHR[i]

       if (length(probeids)>0) {
          if (length(probeids)>1) id = paste(apply(bed[i, probeids], 2, as.character), collapse="_") else id=as.character(bed[i, probeids])
          for (p in 1:length(probeids)) {
              ProbeID = as.character(bed[i,probeids[p]])
              #message("Looping through Probe ID ", ProbeID)
              qtl_index = which(unlist(lapply(listData, function(x) any(bed[,probeids[p]] %in% x$ProbeID))))
              if (length(qtl_index)==0) stop("The ProbeID in the bed file does not match the ProbeIDs in the input data")
              if (length(qtl_index)>1) stop("ProbeID in bed file matches more than one input data")
              if (length(qtl_index)==1) {
                  qtl = listData[[qtl_index]]
                  region.eqtl = qtl[qtl$ProbeID == ProbeID,]
                   data[[qtl_index]] <- region.eqtl
               }
               }
         } else {
         	  id = paste(chrom, pos.start, pos.end, sep="_")
         }     

       message("Looping through ", id)  
       listRegion = list()
       for (j in 1:length(data)) {
         region = data[[j]]
         matches <- which(region$CHR==chrom & region$POS > pos.start & region$POS < pos.end )
         x <- region[matches, ]
         listRegion[[length(listRegion)+1]] <- x
       }

       matches_in_dfs = all(unlist(lapply(listRegion, function(x) nrow(x)>0)))
       if (matches_in_dfs) {

           # Out of each data frame:
           # 1. Remove duplicated IDs
           # 2. keep only common SNPs in all data
           # if (all(c("A1","A2") %in% names(data))) {
           # duplSNPs <- lapply(listRegion, function(x) remove_dupl(data=x, snpcol="SNP"))
           # try = lapply(duplSNPs, function(x) {
           try = list()
           for (d in 1:length(listRegion)) {
                  # x = remove_dupl(data=listRegion[[d]], snpcol="SNP")
                  x = listRegion[[d]]
                  dupl = duplicated(x$SNP)
                  message("Removing ", sum(dupl), " duplicate SNPs")
                  x = x[!dupl,]
                  try[[length(try)+1]] <- x
            }
           listRegion = try

           listRegion <- lapply(listRegion, function(x) x[x$SNP %in% Reduce(intersect, Map("[[", listRegion, "SNP")), ])

           if (nrow(listRegion[[1]])>0) {

           listRegion <- lapply(listRegion, function(df){
                                 df[order(df$SNP),]
                                 })
           # Check that the alleles match and flip based on first data
           # consider A/G T/C, match by strand (complement) 
           if (have_alleles){
               listRegion <- lapply(listRegion, function(df) {df[c("A1", "A2")] <- lapply(df[c("A1", "A2")], as.character); df})
               listRegion <- lapply(listRegion, function(df) {df[c("A1", "A2")] <- lapply(df[c("A1", "A2")], toupper); df})
               listRegion <- lapply(listRegion, change_indels)
               
               for (i in 2:length(listRegion)) {
               match_correct = (listRegion[[1]]["A1"] ==listRegion[[i]]["A1"]) & (listRegion[[1]]["A2"]== listRegion[[i]]["A2"])
               match_flip = (listRegion[[1]]["A1"] == listRegion[[i]]["A2"]) & (listRegion[[1]]["A2"] == listRegion[[i]]["A1"])
               match_comp_one = (listRegion[[1]]["A1"] == complement_snp(listRegion[[i]]["A1"])) & (listRegion[[1]]["A2"]== complement_snp(listRegion[[i]]["A2"]))
               match_comp_two = (listRegion[[1]]["A1"] == complement_snp(listRegion[[i]]["A2"])) & (listRegion[[1]]["A2"] == complement_snp(listRegion[[i]]["A2"]))
               snp_allele_match = match_flip | match_correct | match_comp_one | match_comp_two
               print(listRegion[[i]][!snp_allele_match,])

               if (any(which(match_flip)>0)) {
                 listRegion[[i]][which(match_flip), "A1"]=listRegion[[1]][which(match_flip), "A1"]
                 listRegion[[i]][which(match_flip), "A2"]=listRegion[[1]][which(match_flip), "A2"]
                 listRegion[[i]][which(match_flip), "BETA"]=-listRegion[[i]][which(match_flip), "BETA"]
               }

               if (sum(!snp_allele_match)>0) {
                  removed_snp_list = rbind(removed_snp_list, data.frame(Marker_removed=listRegion[[i]]["SNP"][!snp_allele_match], reason="Alleles do not match"))
                  listRegion[[i]] = listRegion[[i]][which(snp_allele_match),]
               }

             }
             
             listRegion <- lapply(listRegion, function(x) x[x$SNP %in% Reduce(intersect, Map("[[", listRegion, "SNP")), ])
             listRegion <- lapply(listRegion, function(df){
                                 df[order(df$SNP),]
                                 })
          } # end of have_alleles


           nsnps = nrow(listRegion[[1]])
           message("Number of SNPs matching ", nsnps)
           if (nsnps>min_nsnps) {
              if (save.SNP.info) {
                moloc = moloc_test(listRegion, prior_var = prior_var, priors = priors, save.SNP.info=TRUE)  
              } else {
                moloc = moloc_test(listRegion, prior_var = prior_var, priors = priors, save.SNP.info=FALSE) 
              }  
              n <- rownames(moloc$priors_lkl_ppa)
              ppa = signif(moloc$priors_lkl_ppa$PPA, digits=2) # configuration PPA using set priors
              names(ppa) = paste("PPA.",  n, sep="")
              logBF_locus = moloc$priors_lkl_ppa$logBF_locus # locus BF before priors
              names(logBF_locus) = paste("logBF_locus.",  n, sep="")
              bf = moloc$priors_lkl_ppa$sumbf # locus BF before priors
              names(bf) = paste("bf.",  n, sep="")
              
              if (!save.SNP.info) {
                  snp.info = moloc$best_snp
              }
              if (save.SNP.info) {
                  snp.info = moloc$best_snp[[1]]
              }
              best.n =  rownames(snp.info)
              best.snp.PPA= as.data.frame(signif(t(snp.info$coloc_ppas), digits=2))
              names(best.snp.PPA) = paste("best.snp.PPA.",  best.n, sep="")
              best.snp= as.data.frame(t(as.character(snp.info$best.snp.coloc)))
              names(best.snp) = paste("best.snp.",  best.n, sep="")
              if (any(is.na(ppa))) stop("Moloc gives missing values for ", bed[i,])

              nsnps = moloc$nsnps

              # Find min pval. beta, and best snps from the list of all the traits
              if ("PVAL" %in% Reduce(intersect, lapply(listData, colnames))) addpval=TRUE else addpval=FALSE
              if ("BETA" %in% Reduce(intersect, lapply(listData, colnames))) addbeta=TRUE else addbeta=FALSE

              if (addbeta) {
              best.snp.betas.df =c()
              for (d in 1:length(listRegion)) {
                 best.snp.betas =c()
                 for (s in (as.character(snp.info$best.snp.coloc))) {
                 b <- signif(listRegion[[d]]$BETA[listRegion[[d]]$SNP==s], digits=3)
                best.snp.betas = c(best.snp.betas, b)
                }
                best.snp.betas = paste(best.snp.betas, collapse=",")
                best.snp.betas.df = c(best.snp.betas.df,  best.snp.betas)
              }
              names(best.snp.betas.df) = paste("best.snp.betas.", paste("df", 1:length(listRegion), sep=""), sep="")

              best.snp.se.df =c()
              for (d in 1:length(listRegion)) {
                 best.snp.se =c()
                 for (s in (as.character(snp.info$best.snp.coloc))) {
                 b = signif(listRegion[[d]]$SE[listRegion[[d]]$SNP==s], digits=3)
                 best.snp.se = c(best.snp.se, b)
                 }
                best.snp.se = paste(best.snp.se, collapse=",")
                best.snp.se.df = c(best.snp.se.df,  best.snp.se)
              }
              names(best.snp.se.df) = paste("best.snp.se.", paste("df", 1:length(listRegion), sep=""), sep="")
              }
              
              # Find min pval and best snps from the list of all the traits
              if (addpval) {
                minpi = paste(lapply(listRegion, function(x) signif(min(x$PVAL), digits=3)), collapse=",")
                minpi.snp = paste(lapply(listRegion, function(x) as.character(x$SNP[which.min(x$PVAL)])), collapse=",")
                if (addbeta) {
                   minpi.snp.betas = paste(lapply(listRegion, function(x) x$BETA[which.min(x$PVAL)]), collapse=",")
                   minpi.snp.ses = paste(lapply(listRegion, function(x) x$SE[which.min(x$PVAL)]), collapse=",")
                   }
              }

              d <- letters[1:length(listRegion)]
              configs_cases <- do.call(expand.grid, lapply(d, function(x) c("", x)))[-1,]
              configs <- do.call(paste0, configs_cases)
              coloc_configs <- configs[nchar(configs)>1]

              if (any(as.numeric(ppa[paste("PPA.", coloc_configs, sep="")])>0.5)) message("Found a colocalized signal!!!")
              if (save.SNP.info) {
                coloc.out = paste(outfolder, "/", prefix, ".output.perSNP/", sep="")
                if (!file.exists(coloc.out)) dir.create(coloc.out)
                write.table(x=moloc[[4]], file=paste(coloc.out, id, '_results.tab', sep=''),row.names = FALSE, quote = FALSE, sep = '\t')

              }

              res = cbind(id, CHR = chrom, START=pos.start, STOP=pos.end, nsnps, t(logBF_locus), t(bf), t(ppa), best.snp, best.snp.PPA)
              if (addpval) res = cbind(res, minpi, minpi.snp)
              if (addpval & addbeta) res = cbind(res, minpi.snp.betas, best.snp.se)
              if (addbeta) res = cbind(res, best.snp.betas, best.snp.se)
              
              #res.all = rbind(res.all, res)
              

####***### COLOC OLD

if (compare_coloc) {
    print("Note: This computes COLOC using the first two datasets. When using two traits, PPA.ab and COLOC_ab are slightly diffeent becasue of approximations.")
    if (length(listData) > 2) warning("Using more than two traits, COLOC and moloc are not directly comparable")
    require(coloc)
    # useBETA = TRUE
    p1_coloc = 1e-04; p2_coloc = 1e-04; p12_coloc = 10^-5
    
    # Add type
    for (j in 1:length(listRegion)) {
        listRegion[[j]]$type <- ifelse("Ncases" %in% names(listRegion[[j]]), "cc", "quant")
        if (unique(listRegion[[j]]$type)=="quant") listRegion[[j]]$s = rep(0.5, length(listRegion[[j]]$N))
        if (unique(listRegion[[j]]$type)=="cc") listRegion[[j]]$s = listRegion[[j]]$Ncases/listRegion[[j]]$N
    }
    if (!haveBETA | forcePVAL) {
        dataset.biom = list(snp = listRegion[[1]]$SNP, pvalues = listRegion[[1]]$PVAL,
        N = listRegion[[1]]$N, s = listRegion[[1]]$s, type = unique(listRegion[[1]]$type), MAF=listRegion[[1]]$MAF)
        dataset.eqtl = list(snp = listRegion[[2]]$SNP, pvalues = listRegion[[2]]$PVAL,
        N = listRegion[[2]]$N, s = listRegion[[2]]$s, type = unique(listRegion[[2]]$type), MAF=listRegion[[2]]$MAF)
    }
    if (haveBETA & !forcePVAL) {
        dataset.biom = list(snp = listRegion[[1]]$SNP, beta = listRegion[[1]]$BETA, varbeta= (listRegion[[1]]$SE)^2, s=listRegion[[1]]$s, type = unique(listRegion[[1]]$type), MAF=listRegion[[1]]$MAF, N=listRegion[[1]]$N) #, sdY=unique(merged.data$sdY.biom))
        dataset.eqtl = list(snp = listRegion[[2]]$SNP, beta = listRegion[[2]]$BETA, varbeta= (listRegion[[2]]$SE)^2, s=listRegion[[2]]$s, type = unique(listRegion[[2]]$type), MAF=listRegion[[2]]$MAF, N=listRegion[[2]]$N)
    }

    ### COLOC OLD
    capture.output(coloc.res <- coloc.abf(dataset.biom, dataset.eqtl, p1 = p1_coloc, p2 = p2_coloc, p12 = p12_coloc))
    pp0       <- as.numeric(coloc.res$summary[2])
    pp1       <- as.numeric(coloc.res$summary[3])
    pp2       <- as.numeric(coloc.res$summary[4])
    pp3       <- as.numeric(coloc.res$summary[5])
    pp4       <- as.numeric(coloc.res$summary[6])
    
    coloc.res = data.frame(COLOC_zero=pp0, COLOC_a=pp1, COLOC_b=pp2, COLOC_a.b=pp3, COLOC_ab=pp4)
    res = cbind(res, coloc.res)
    
}

###########*******


           } # if nsnps>50 # after finding common alelles if have_alleles = TRUE
           } # if nrow(listRegion[[1]])>0  # first matching of SNPs across data frames
           } #  matches_in_dfs
           
           if(nrow(res)==0){
           return(NULL)
           }
           return(res)

}
  
  
  if (write) {
       write.table(x =  res.all , file = outfname, row.names = FALSE, quote = FALSE, sep = '\t')
       if (nrow(removed_snp_list)>0) write.table(x = removed_snp_list, file= paste(outfolder, prefix, '_removed.snps', sep=''), row.names = FALSE, quote = FALSE, sep = '\t')
       #if (nrow(res.all)==0) 
       #message("SNPs in common do not reach ", min_nsnps, " in any region: consider lowering this threshold")
    }
  
  return(res.all)

}
