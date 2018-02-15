# Remove duplicated IDs (duplicated snps and indels, keep only snps if have allele info)
remove_dupl = function(data=data, snpcol = "SNP") {
    removed_list <- data.frame()
    n_occur <- data.frame(table(data[,snpcol]))
    dupl = data[data[,snpcol] %in% n_occur$Var1[n_occur$Freq > 1],]
         while (nrow(dupl)>0) {
          #removed_list <- rbind(removed_list, data.frame(Marker_removed = dupl$SNPID, reason = "Duplicated SNPs"))
          if (length(grep("A1|A2", names(data), ignore.case=T))==2) {
             dupl <- transform(dupl, n=nchar(as.character(dupl$A1)) + nchar(as.character(dupl$A2)))
             dupl=dupl[order(dupl$n, decreasing=T),]
          } else {
             dupl=dupl[order(dupl$MAF, decreasing=T),]
          }
          toremove = rownames(dupl[ !duplicated(dupl[,snpcol]), ])
          if (length(toremove)>0) {
            removed_list <- rbind(removed_list, data.frame(Marker_removed = dupl[,snpcol][!duplicated(dupl[,snpcol])], reason = "Duplicated SNPs"))
          data = data[!(rownames(data) %in% toremove),]
          message("To remove: ", length(toremove), " duplicated SNP names")
          }  else {
          removed_list <- data.frame(Marker_removed = NA, reason = "Duplicated SNPs")
          }
          dupl = data[data[,snpcol] %in% n_occur$Var1[n_occur$Freq > 1],]
          }
    return(list(data, removed_list))
}