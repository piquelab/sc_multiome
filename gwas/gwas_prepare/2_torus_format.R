##
library(tidyverse)
library(data.table)
options(scipen=13)

rm(list=ls())


### parsing argument
args=commandArgs(trailingOnly=T)
if ( length(args)>0){
   trait <- args[1]
}else{
   trait <- "CARDIoGRAM_C4D_CAD"
}   

### create directory 
outdir <- "./gwas_torusfile/"
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)


######################################################
### Generate torus files with LD block information ###
######################################################


### read summary data
fn <- paste("./gwas_0/", trait, "_gwas.txt.gz", sep="")  
summ <- fread(fn, header=T, data.table=F)


### LD block files
bed <- fread("eur_ld.hg38.bed", header=T, fill=T, data.table=F)
bed <- bed%>%drop_na(start, stop)



### add locus using LD block file
nblock <- nrow(bed)
res_torus <- lapply(1:nblock, function(i){
   ##
   cat("Locus", i, "\n")
   ## 
   chr_i <- bed$chr[i]
   s0 <- as.numeric(bed$start[i])
   s1 <- as.numeric(bed$stop[i])

   ## SNPs in each block   
   summ2 <- summ%>%dplyr::filter(chr==chr_i, pos>=s0, pos<s1)
   summ2 <- summ2%>%mutate(Locus=paste("Loc", i, sep=""))%>%dplyr::select(id_b38, Locus, zscore)
   summ2     
})
 
##
res_torus <- do.call(rbind, res_torus)

### output 
opfn <- gzfile(paste0(outdir, trait, "_torus.txt.gz"), "w")
write.table(res_torus, opfn, quote=F, row.names=F, col.names=F)


###END
