###
library(tidyverse)
library(data.table)

rm(list=ls())

options(scipen=200)


### passing argument
## args=commandArgs(trailingOnly=T)
## if (length(args)>0){
##    condition <- args[1]
##  }else{
##    condition <- "Bcell_CTRL"
## }

outdir <- "./torus_input/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)



################################
### asthma gwas summary data ###
################################

### add chr_pos information once and don't need run it again. 
fn <- "./gwas_input/Asthma_torus_zval.txt.gz"
summ <- fread(fn, header=T, data.table=F, fill=T)
snpSel <- summ$panel_variant_id

 
###
summ2 <- summ%>%dplyr::select(SNP=panel_variant_id, Locus, zscore)
gfn <- gzfile("./gwas_input/Asthma.torus.zval.gz")
write.table(summ2, gfn, row.names=F, col.names=F, quote=F)



####
#### annotation 

### snp list
traits <- read.table("traits_ls.txt")$V1
snpList <- lapply(traits, function(ii){
   ##
    fn <- paste("./gwas_input/", ii, ".torus.zval.gz", sep="")
    xx <- fread(fn, header=F, data.table=F)
    snp <- as.character(xx[,1])
    cat(ii, length(snp), "\n")
    snp
})
snpList <- do.call(c, snpList)
snpList <- unique(snpList)
opfn <- "./gwas_input/snpList.txt"
write.table(snpList, file=opfn, row.names=F, col.names=F, quote=F)



###
###
snpSel <- fread("./gwas_input/snpList.txt", header=F)$V1

fn <- "/nfs/rprdata/julong/sc_multiome/genetic_variant_annot/4_SNPAnnot.outs/zzz_th0.1_multiConditions_torus.annot.gz"
annot <- fread(fn, header=T, data.table=F)   
annot2 <- annot%>%filter(SNP%in%snpSel)

x <- annot2


##
nsnp <- colSums(x[,-1])
colSel <- names(nsnp)[nsnp>5000]
 
x2 <- x[,c("SNP", colSel)]

gfn <- gzfile("./torus_input/zzz_torus.annot.gz")
write.table(x2, gfn, row.names=F, col.names=T, quote=F)






## End, Mar-30-2023











       
 
