##
library(tidyverse)
library(data.table)

rm(list=ls())


### passing argument
args=commandArgs(trailingOnly=T)
if (length(args)>0){
   trait <- args[1]
 }else{
   trait <- "CARDIoGRAM_C4D_CAD"
}

## traits <- read.table("traits_of_interest.txt")$V1
## traits <- sort(traits)



###
### 0, fastQTL results
fast <- fread("./eQTL_results/Whole_Blood.fastqtl.gz", header=F, data.table=F) ##9,719,090 SNPs
names(fast) <- c("gene", "id_b38", "DTSS", "pval", "beta", "se")
fast <- fast%>%mutate(gene_SNP=paste(gene, id_b38, sep="_"))
fast <- fast%>%filter(!grepl("chrX", id_b38))

## pval_eqtl <- fast$pval
## names(pval_eqtl) <- fast$gene_SNP


###
### fine-mapping,  dap results from multiple condition annotation, default
fn <- "./eQTL_results/WBL_PIP_allSNPs.txt.gz"
dap <- fread(fn, header=F, data.table=F)
names(dap) <- c("gene", "id_b38", "PIP")
dap <- dap%>%mutate(gene_SNP=paste(gene, id_b38, sep="_")) 
dap <- dap%>%filter(!grepl("chrX", id_b38))
PIP <- dap$PIP
names(PIP) <- dap$gene_SNP


###
### fine-mapping, dap results from baseline annotation
fn <- "./eQTL_results/WBL_PIP_allSNPs_base.txt.gz"
dap2 <- fread(fn, header=F, data.table=F)
names(dap2) <- c("gene", "id_b38", "PIP")
dap2 <- dap2%>%mutate(gene_SNP=paste(gene, id_b38, sep="_"))
dap2 <- dap2%>%filter(!grepl("chrX", id_b38))
PIP2 <- dap2$PIP
names(PIP2) <- dap2$gene_SNP


###
### filtering gene_SNP
fast <- fast%>%filter(gene_SNP%in%dap$gene_SNP)



###
### gwas data
fn <- paste0("./gwas_data/", trait, "_impute_gwas.txt.gz")
summ <- fread(fn, header=T, data.table=F)
summ <- summ%>%dplyr::select(id_b38, chr, pos, zscore_gwas=zscore, pval_gwas=pval, id_b38_2)

## pval_gwas <- summ$pval
## names(pval_gwas) <- as.character(summ$id_b38) ##summ$panel_variant_id


###
### combine fastQTL and fine-mapping together 
res_comb <- fast%>%dplyr::select(gene, id_b38, gene_SNP, pval_eqtl=pval)%>%
    mutate(PIP_conditions=PIP[gene_SNP], PIP_base=PIP2[gene_SNP])




###########################################################################
### TWAS minimum p value to select the reprentative SNPs for each gene ####
###########################################################################

###
outdir <- "./1_SMR_output/eQTL_fastQTL/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


res_gene <- res_comb%>%group_by(gene)%>%slice_min(order_by=pval_eqtl, n=1)%>%ungroup()%>%as.data.frame() 
res_gene <- res_gene%>%left_join(summ, by="id_b38")

res2_gene <- res_gene%>%group_by(gene)%>%slice_min(order_by=pval_gwas, n=1, with_ties=F)%>%
    ungroup()%>%as.data.frame()%>%drop_na(pval_gwas)
 
gfn <- gzfile(paste(outdir, trait, "_minP_twas.txt.gz", sep=""))
write.table(res2_gene, gfn, row.names=F, col.names=T, quote=F)




 
#######################################################################
### TWAS top PIP from to select the reprentative SNPs for each gene ###
###       integration multiple conditions annotation                ###
#######################################################################

outdir <- "./1_SMR_output/eQTL_dap_conditions/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)
 
res_gene <- res_comb%>%drop_na(PIP_conditions)%>%
    group_by(gene)%>%slice_max(order_by=PIP_conditions, n=1)%>%ungroup()%>%as.data.frame()
res_gene <- res_gene%>%group_by(gene)%>%slice_min(order_by=pval_eqtl, n=1)%>%ungroup()%>%as.data.frame()
res_gene <- res_gene%>%left_join(summ, by="id_b38")

res2_gene <- res_gene%>%group_by(gene)%>%slice_min(order_by=pval_gwas, n=1, with_ties=F)%>%
    ungroup()%>%as.data.frame()%>%drop_na(pval_gwas)
     
gfn <- gzfile(paste(outdir, trait, "_topPIP_twas.txt.gz", sep=""))
write.table(res2_gene, gfn, row.names=F, col.names=T, quote=F)



#######################################################################
### TWAS top PIP from to select the reprentative SNPs for each gene ###
### integration baseline fine-mapping results                       ###
#######################################################################


outdir <- "./1_SMR_output/eQTL_dap_base/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)

 
res_gene <- res_comb%>%drop_na(PIP_base)%>%
    group_by(gene)%>%slice_max(order_by=PIP_base, n=1)%>%ungroup()%>%as.data.frame()
res_gene <- res_gene%>%group_by(gene)%>%slice_min(order_by=pval_eqtl, n=1)%>%ungroup()%>%as.data.frame()
res_gene <- res_gene%>%left_join(summ, by="id_b38")

res2_gene <- res_gene%>%group_by(gene)%>%slice_min(order_by=pval_gwas, n=1, with_ties=F)%>%
    ungroup()%>%as.data.frame()%>%drop_na(pval_gwas)
     
gfn <- gzfile(paste(outdir, trait, "_topPIP_twas.txt.gz", sep=""))
write.table(res2_gene, gfn, row.names=F, col.names=T, quote=F)



###END
    
    
