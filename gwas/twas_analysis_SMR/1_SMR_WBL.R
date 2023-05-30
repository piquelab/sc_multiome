##
library(tidyverse)
library(data.table)

rm(list=ls())


### passing argument
## args=commandArgs(trailingOnly=T)
## if (length(args)>0){
##    condition <- args[1]
##    dir <- args[2]
##  }else{
##    condition <- "Bcell_CTRL"
##    dir <- "SCAIP_results"
## }




###
### 0, fastQTL results
fast <- fread("./eQTL_results/Whole_Blood.fastqtl.gz", header=F, data.table=F) ##9,719,090 SNPs
names(fast) <- c("gene", "id_b38", "DTSS", "pval", "beta", "se")
fast <- fast%>%mutate(gene_SNP=paste(gene, id_b38, sep="_"))
pval_eqtl <- fast$pval
names(pval_eqtl) <- fast$gene_SNP



###
### 1. PIP results  add grch38 position   
fn <- "./eQTL_results/WBL_PIP_allSNPs.txt.gz"
res <- fread(fn, header=F, data.table=F)
names(res) <- c("gene", "id_b38", "PIP")


###
outdir <- "./1_SMR_output/eQTL_conditions/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


### SMR analysis 
traits <- read.table("./traits_ls.txt")$V1

for (trait in traits){
##    
fn <- paste0("./gwas_data/gwas_impute/", trait, "_gwas.txt.gz")
summ <- fread(fn, header=T, data.table=F)
names(summ) <- c("id_b38", "chr", "pos", "zscore", "pvalue", "id_b38_2")
    
pval <- summ$pvalue
names(pval) <- as.character(summ$id_b38) ##summ$panel_variant_id


res$pval_gwas <- pval[as.character(res$id_b38)]
res2 <- res%>%drop_na(pval_gwas)


###
### high PIP
res2_gene <- res2%>%group_by(gene)%>%slice_max(order_by=PIP, n=1)%>%ungroup()%>%as.data.frame()
res2_gene <- res2_gene%>%mutate(gene_SNP=paste(gene, id_b38, sep="_"))

###
## if ties, we pick the SNP with minimum p-value
res2_gene$pval_eqtl <- pval_eqtl[as.character(res2_gene$gene_SNP)]
res3_gene <- res2_gene%>%group_by(gene)%>%slice_min(order=pval_eqtl, n=1)%>%
    ungroup()%>%as.data.frame()

### if ties, we pick the SNP with minimum of gwas 
res3_gene <- res3_gene%>%group_by(gene)%>%
    slice_min(order_by=pval_gwas, n=1, with_ties=F)%>%ungroup()%>%
    mutate(gene2=gsub("\\..*", "", gene))%>%
    as.data.frame() ### 19689 

gfn <- gzfile(paste(outdir, trait, "_topPIP_twas.txt.gz", sep=""))
write.table(res3_gene, gfn, row.names=F, col.names=T, quote=F)

cat(trait, nrow(res3_gene), "\n")
    
}    
### End


###
### correct pvalue, *2

## traits <- read.table("./traits_ls.txt")$V1
## dap <- c("eQTL_base", "eQTL_conditions", "eQTL_union")

## for (ii in dap[-1]){
## ## 
## for (trait in traits[-1]){
##    ###
##    fn <- paste("./1_SMR_output/", ii, "/", trait, "_topPIP_twas.txt.gz", sep="")
##    summ <- fread(fn, header=T, data.table=F)
##    summ <- summ%>%mutate(pval_gwas=pval_gwas*2)
##    ###
##    gfn  <- gzfile(paste("./1_SMR_output/", ii, "/", trait, "_topPIP_twas.txt.gz", sep=""))
##    write.table(summ, gfn, row.names=F, col.names=T, quote=F)
##    cat(ii, trait, "\n")
## } ### traits
    
## } ### dap    
    
    
