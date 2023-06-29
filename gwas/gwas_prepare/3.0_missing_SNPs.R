##
library(tidyverse)
library(data.table)

rm(list=ls())

####
####

traits <- read.table("traits_of_interest.txt")$V1
traits <- sort(traits)


###
###
snp_WBL <- read.table("WBL_snplist.txt")$V1
for (ii in traits){
   ##
   outdir2 <- paste0("./gwas_imputefile/", ii, "/")
   if ( !file.exists(outdir2)) dir.create(outdir2, showWarnings=F, recursive=T)
    
   ###
   fn <- paste("./gwas_0/", ii, "_gwas.txt.gz", sep="")
   summ <- fread(fn, header=T, data.table=F)
   x <- summ$id_b38
   ###
   snp_miss <- setdiff(snp_WBL, x)
   ##
   opfn <- paste(outdir2, "zzz_allmissing.txt", sep="")
   write.table(snp_miss, file=opfn, row.names=F, col.names=F, quote=F)

   cat(ii, length(snp_miss), "\n")
}    


####
