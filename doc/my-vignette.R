## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----three-tables-single, echo=FALSE, results='asis'---------------------
library(knitr)
options(scipen = 1, digits = 2)
## load single locus data (in a list) and bed file
data_single=get(load(file="../data/data_single.rda"))
t1 = knitr::kable(head(data_single[[1]], 5),  format='html', table.attr='cellpadding="1"', output = FALSE)

t2 = knitr::kable(head(data_single[[2]], 5),  format='html', table.attr='cellpadding="1"', output = FALSE)

t3 = knitr::kable(head(data_single[[3]], 5),  format='html', table.attr='cellpadding="1"', output = FALSE)
cat(c('<table><tr valign="top"><td>', t1, '</td><td>', t2, '</td><td>', t3, '</td><tr></table>'),
    sep = '')
## 

## ----message=FALSE, warning=F--------------------------------------------
options(scipen = 1, digits = 2)
library(moloc)
moloc <- moloc_test(data_single, prior_var=c(0.01, 0.1, 0.5), priors=c(1e-04, 1e-06, 1e-07))
# Posteriors
print(moloc[[1]])
# Number of SNPs analyzed
print(moloc[[2]])
# Posterior of the most likely SNP co-localizing with another trait
print(moloc[[3]])

## ----three-tables-genome, echo=FALSE, results='asis'---------------------
library(knitr)
options(scipen = 1, digits = 2)
## load single locus data (in a list) and bed file
library(moloc)
data_genome=get(load(file="../data/data_genome.rda"))
t1 = knitr::kable(head(data_genome[[1]], 5),  format='html', table.attr='cellpadding="1"', output = FALSE)
t2 = knitr::kable(head(data_genome[[2]], 5),  format='html', table.attr='cellpadding="1"', output = FALSE)
t3 = knitr::kable(head(data_genome[[3]], 5),  format='html', table.attr='cellpadding="1"', output = FALSE)
cat(c('<table><tr valign="top"><td>', t1, '</td><td>', t2, '</td><td>', t3, '</td><tr></table>'),
    sep = '')
## 

## ----echo=FALSE----------------------------------------------------------
knitr::kable(bed[1:2,])
## 

## ----message=FALSE, warning=F--------------------------------------------
options(scipen = 1, digits = 2)
library(moloc)
library(foreach)
library(doParallel)
res <- coloc.genome(data_genome, bed, prefix = "pref", save.SNP.info=FALSE, cores=20, have_alleles=TRUE, bychrpos=TRUE, prior_var="default", priors=c(1e-04, 1e-06, 1e-07), min_nsnps = 50, takelog = FALSE, write=TRUE, outfolder = "test", forcePVAL = FALSE)
# Posteriors
print(res[[1]])
# Bed file of loci analyzed
print(res[[2]])

