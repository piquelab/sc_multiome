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


####################
### torus format ###
####################

### eQTL summary results
fn <- "./eQTL_results/Whole_Blood.allpairs.txt.gz"
res <- fread(fn, header=T, data.table=F, stringsAsFactors=F)
zscore <- res$slope/res$slope_se
summ <- data.frame(SNP=res$variant_id,
   gene=gsub("\\..*", "", res$gene_id), beta=res$slope, "t-stat"=zscore, "p-value"=res$pval_nominal)

##
##summ2 <- summ%>%mutate(gene=gsub("\\..*", "", gene))
 
opfn <- gzfile(paste(outdir, "Whole_Blood.eQTL.txt.gz", sep=""))
write.table(summ, opfn, quote=F, row.names=F)





##########################
### gene map for torus ###
##########################

###
### annotation files
## anno <- read.table("/wsu/home/groups/piquelab/data/gencode/Gencode_human/release_31/GRCh37_mapping/gencode.v31lift37.annotation.gff3.gz", header=F, stringsAsFactors=F) ##grch37

###
### gene annotation files from grch38

fn <- "/wsu/home/groups/piquelab/data/gencode/Gencode_human/release_31/gencode.v31.annotation.gff3.gz"
## fn <- "/wsu/home/groups/piquelab/data/gencode/Gencode_human/release_42/gencode.v42.chr_patch_hapl_scaff.annotation.gff3.gz"
anno <- read.table(fn, header=F, stringsAsFactors=F)
 
anno2 <- anno%>%
   mutate(gene_id=gsub("ID=|;.*", "", V9), ID=gsub("\\..*", "", gene_id))%>%
   dplyr::filter(V3=="gene")%>%
   dplyr::rename("Chr"="V1", "min"="V4", "max"="V5", "strand"="V7")%>%
   mutate(Chr=gsub("chr", "", Chr),
          start=ifelse(strand=="+", min, max),
          end=ifelse(strand=="-", max, min))

###
### gene used for GTEx 
geneList <- read.table("./eQTL_results/geneList.txt")$V1 ### 20,315 gene
geneDf <- data.frame(gene_id2=unique(geneList))%>%mutate(ID=gsub("\\..*", "", gene_id2))
 
anno3 <- anno2%>%inner_join(geneDf, by="ID")%>%dplyr::filter(!grepl("PAR_Y", gene_id)) ## 20207


###
### gene map file
gene.map <- data.frame(gene=anno3$gene_id2, chr=anno3$Chr, start=anno3$start, start2=anno3$start)
###
opfn <- gzfile(paste(outdir, "zzz_gene.map.gz", sep=""))
write.table(gene.map, opfn, sep="\t", quote=F, row.names=F, col.names=F)




#########################################
### SNP map file and annotation files ###
#########################################

fn <- "./eQTL_results/snpList.txt"
snpSel <- read.table(fn)$V1


fn <- "/nfs/rprdata/julong/sc_multiome/genetic_variant_annot/4_SNPAnnot.outs/zzz_th0.1_multiConditions_torus.annot.gz"
x <- fread(fn, header=T, data.table=F)
x2 <- x%>%filter(SNP%in%snpSel)

nsnp <- colSums(x2[,-1])
colSel <- names(nsnp)[nsnp>5000]

x3 <- x2[,c("SNP", colSel)]
opfn <- gzfile("./torus_input/zzz_torus.annot.gz")
write.table(x3, opfn, quote=F, row.names=F, col.names=T)


###
### snp map files
snp.map <- data.frame(SNP=x2$SNP)%>%
   separate(SNP, into=c("chr", "pos", "ref", "alt", "hg"), sep="_")%>%
   mutate(chr2=as.character(gsub("chr", "", chr)))

##
snp.map2 <- cbind(x2$SNP, snp.map[,c("chr2", "pos")])
opfn2 <- gzfile(paste(outdir, "zzz_snp.map.gz", sep=""))
write.table(snp.map2, opfn2,  quote=F, row.names=F, col.names=F)


###
###
fn <- "./torus_input/zzz_torus.annot.gz"
x <- fread(fn, header=T, data.table=F)


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






