##
library(tidyverse)
library(data.table)
library(annotables)
library(INTACT, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")

###
### intersect SNPs
snp_gwas <- unique(read.table("../gwas/snpList_uniq.txt")$V1)
snp_eqtl <- unique(read.table("../eQTL_results/snpList.txt")$V1)
shared <- intersect(snp_gwas, snp_eqtl) ## 9,428,357


###
outdir <- "./INTACT_output/"
if ( !file.exists(outdir)) dir.create(outdir, recursive=T)



#######################
### INTACT analysis ###
#######################

traits <- sort(read.table("../traits_of_interest.txt")$V1)
 
## autosome <- as.character(1:22)
## annot <- grch38%>%dplyr::filter(chr%in%autosome)%>%dplyr::select(ensgene, symbol)


####
###
for ( ii in traits){

##    
fn <- paste("./enloc_output/", ii, ".enloc.gene.out", sep="")
enloc <- read.table(fn, header=T)


## add zscore    
fn <- paste("../../twas_analysis_SMR/analysis_correct/1_SMR_output/eQTL_dap_conditions/", ii, "_topPIP_twas.txt.gz", sep="")    
twas <- fread(fn, header=T, data.table=F)
twas <- twas%>%mutate(gene=gsub("\\..*", "", gene))    
twas2 <- twas%>%dplyr::select(Gene=gene, zscore=zscore_gwas)
 
## combine enloc and z-score of twas
DF <- enloc%>%inner_join(twas2, by="Gene")
    

## intact analysis
res_intact <- intact(GLCP=DF$GLCP, z_vec=DF$zscore)
DF$PCG <- res_intact
### FDR 
DF2 <- DF%>%mutate(LFDR=1-PCG)%>%arrange(LFDR)
x <- DF2$LFDR
FDR <- cumsum(x)/1:length(x)
DF2$FDR <- FDR
### 
opfn <- paste(outdir, ii, "_intact.txt", sep="")
write.table(DF2, opfn, row.names=F, quote=F)

###
cat(ii, nrow(DF2), "\n")    

}

## gene1 <- DF2%>%dplyr::filter(FDR<0.1)%>%pull(Gene)

## fn <- "/nfs/rprdata/julong/sc-atac/twas_analysis_2022-10-16/enloc_analysis/ALOFT_intact.txt"
## gene2 <- read.table(fn, header=T)%>%filter(FDR<0.1)%>%pull(Gene)

## fn <- "/nfs/rprdata/julong/sc-atac/twas_analysis_2022-10-16/enloc_GTEx_V8/GTEx_v8_intact.txt"
## gene3 <- read.table(fn, header=T)%>%filter(FDR<0.1)%>%pull(Gene)

## x12 <- intersect(gene1, gene2)
## x13 <- intersect(gene1, gene3)


