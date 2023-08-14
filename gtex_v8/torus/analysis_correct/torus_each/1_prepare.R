#
library(tidyverse)
library(data.table)

options(scipen=18)


rm(list=ls())

### passing argument
## args=commandArgs(trailingOnly=T)
## if (length(args)>0){
##    condition <- args[1]
##  }else{
##    condition <- "Bcell_CTRL"
## }

outdir <- "./torus_input/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)
 

###
###
x <- fread("../torus_input/zzz_torus.annot.gz", header=T, data.table=F)

## x2 <- x[,-c(1,2)]
## nsnp <- colSums(x2)
 
colSel <- colnames(x)[-c(1,2)]
colSel2 <- gsub("_d$", "", colSel)
write.table(colSel2, "annot_ls.txt", quote=F, row.names=F, col.names=F)


for (ii in colSel){
    ##
    x2 <- x[,c("SNP", "peaking_d", ii)]
    ##
    ii2 <- gsub("_d$", "", ii)
    opfn <- paste(outdir, ii2, ".torus.gz", sep="")
    fwrite(x2, opfn, quote=F, sep=" ")
    cat(ii2, "\n")
}






