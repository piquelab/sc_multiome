##
library(tidyverse)
library(data.table)
options(scipen=13)

rm(list=ls())

## combine files

traits <- read.table("traits_of_interest.txt")$V1
traits <- sort(traits)

for (ii in traits){
   ##
   fn1 <- paste("./gwas_0/", ii, "_gwas.txt.gz", sep="")
   x1 <- fread(fn1, header=T, data.table=F)
   x1 <- x1%>%mutate(id_b38_2=id_b38) 
   ##
   fn2 <- paste("./gwas_imputefile/", ii, "/gwas_concate.txt.gz", sep="")
   x2 <- fread(fn2, header=F, data.table=F)
   names(x2) <- names(x1) 

   ### combine and output 
   comb <- rbind(x1, x2)
   opfn <- gzfile(paste("./gwas_imputefile/", ii, "_impute_gwas.txt.gz", sep=""))
   write.table(comb, file=opfn, row.names=F, col.names=T, quote=F)
   cat(ii, nrow(comb), "\n")
}    
    
    
