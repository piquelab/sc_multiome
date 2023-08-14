

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




#########################################
### SNP map file and annotation files ###
#########################################

###
fn <- "../eQTL_results/snpList.txt"
snpSel <- read.table(fn)$V1

###
### multiple condition response motifs annotation
fn <- "/nfs/rprdata/julong/sc_multiome/genetic_variant_annot/4_SNPAnnot_correct/zzz_th0.1_multiConditions_torus.annot.gz"
x <- fread(fn, header=T, data.table=F)
x2 <- x%>%filter(SNP%in%snpSel)

nsnp <- colSums(x2[,-1])
colSel <- names(nsnp)[nsnp>5000]

x3 <- x2[,c("SNP", colSel)]
opfn <- paste(outdir, "zzz_torus.annot.gz", sep="")
fwrite(x3, opfn, quote=F, sep=" ")


###
### snp map files
## snp.map <- data.frame(SNP=x2$SNP)%>%
##    separate(SNP, into=c("chr", "pos", "ref", "alt", "hg"), sep="_")%>%
##    mutate(chr2=as.character(gsub("chr", "", chr)))

## ##
## snp.map2 <- cbind(x2$SNP, snp.map[,c("chr2", "pos")])
## opfn2 <- paste(outdir, "zzz_snp.map.gz", sep="")
## fwrite(snp.map2, opfn2, quote=F, sep=" ")


###
### union annotation
fn <- "/nfs/rprdata/julong/sc_multiome/genetic_variant_annot/4_SNPAnnot_correct/zzz_th0.1_union_torus.annot.gz"
x <- fread(fn, header=T, data.table=F)
x2 <- x%>%filter(SNP%in%snpSel)

opfn3 <- paste(outdir, "zzz2_union_torus.annot.gz", sep="")
fwrite(x2, file=opfn3, quote=F, sep=" ")





## fn <- paste(outdir, "zzz_combine2_torus.annot.gz", sep="")
## x <- fread(fn, header=F)
## names(x) <- names(anno)
## gfn <- gzfile(paste(outdir, "zzz_combine2_torus.annot.gz", sep=""))
## write.table(x, gfn, quote=F, row.names=F, col.names=T)


###############
### summary ###
###############

## ### eqtl
## fn <- "./eQTL_results/snpList.txt"
## snp <- fread(fn, header=F)
## snp <- unique(snp$V1)

## ## gwas
## fn <- "../../gwas/torus/gwas_input/Asthma_torus_zval.txt.gz"
## summ <- fread(fn, header=T, data.table=F)
## snp <- unique(summ$panel_variant_id)


## ## fn <- "./torus_input/zzz_multiConditions_torus.annot.gz"
## ## x <- fread(fn, header=T, data.table=F)
## fn <- "../../genetic_variant_annot/gtex_v8_snpinfor.txt.gz"
## x <- fread(fn, data.table=F)

## ## 
## shared <- intersect(snp, x$V3)


















## res <- fread(fn, header=T, data.table=F, stringsAsFactors=F)






