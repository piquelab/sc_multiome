##
library(Matrix)
library(tidyverse)
library(data.table)
## library(clusterProfiler)
## library(org.Hs.eg.db)
library(annotables)
library(ggplot2)
library(cowplot)
library(RColorBrewer)

rm(list=ls())


###
### calculate FDR


###
outdir <- "./5_summary.outs/Whole_Blood//"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


gene_all <- read.table("../geneList/geneList.txt")$V1

#########################################
### Calculate local FDR for each geen ###
#########################################

### use all gene to adjust for FDR

fn <- paste(outdir, "WBL_gene_lfdr.txt", sep="")
dap <- fread(fn, header=F, data.table=F)
names(dap) <- c("gene", "lfdr")
 
dap <- dap%>%mutate(lfdr=as.numeric(lfdr))%>%arrange(lfdr)
x <- dap$lfdr
FDR <- cumsum(x)/1:length(x)
dap$FDR <- FDR
dap <- dap%>%mutate(gene2=gsub("\\..*", "", gene))

### 
opfn2 <- paste(outdir, "WBL_gene_FDR.txt", sep="")
write.table(dap, opfn2, row.names=F, quote=F, sep="\t")


###
### extract protein conding genes

fn <- paste(outdir, "WBL_gene_FDR.txt", sep="")
dap <- read.table(fn, header=T)


autosome <- as.character(1:22) 
grch38_unq <- grch38%>%
    dplyr::filter(chr%in%autosome, grepl("protein", biotype))%>%
    distinct(ensgene, chr, .keep_all=T)%>%dplyr::select(gene=ensgene, chr, biotype)

dap2 <- dap%>%filter(gene2%in%grch38_unq$gene)
opfn <- paste(outdir, "WBL_gene_FDR_2.txt", sep="")
write.table(dap2, file=opfn, row.names=F, quote=F, sep="\t") ## protein coding and autosome genes


egene <- dap2%>%filter(FDR<0.1)%>%pull(gene)%>%unique()
opfn <- "./5_summary.outs/Whole_Blood/WBL_egene.txt"
write.table(egene, file=opfn, quote=F, row.names=F, col.names=F) ### 9,773 egene

# #
## egene <- read.table("./5_summary.outs/Whole_Blood/WBL_gene_FDR_2.txt",header=T)%>%
##     filter(FDR_2<0.1)%>%pull(gene)%>%unique()  ##9906

## opfn <- "./5_summary.outs/Whole_Blood/WBL_egene.txt"
## write.table(egene, file=opfn, quote=F, row.names=F, col.names=F)


## ### use all the gene to adjust FDR
## fn <- "./5_summary.outs/Whole_Blood/WBL_gene_FDR_2.txt"
## x <- read.table(fn, header=T)
## egene <- x%>%filter(FDR<0.1)%>%pull(gene)%>%unique()
## opfn <- "./5_summary.outs/Whole_Blood/WBL_egene.txt"
## write.table(egene, file=opfn, quote=F, row.names=F, col.names=F)




##############################################
### expected number of independent signals ###
##############################################

fn <- "./5_summary.outs/Whole_Blood/WBL_egene.txt"
egenes <- read.table(fn)$V1

## testing genes, 14,084
##
plotDF <- map_dfr(egenes, function(ens){
   ##    
   fn <- paste("./dap-g_outs/Whole_Blood/", ens, ".cluster.out",  sep="")
   cluster <- read.table(fn)
   df2 <- data.frame(gene=ens, neqtl=sum(cluster$V3))
   df2
})


p <- ggplot(plotDF, aes(x=as.numeric(neqtl)))+
    geom_histogram(color="grey30", fill="grey80", bins=50)+
    scale_x_continuous("Expected number of eQTLs per eGene", breaks=seq(0, 10, by=2))+
    ylab("Number of eGenes")+    
##    facet_wrap(~treat2, nrow=5)+ ##, scales="free")+
    ggtitle("WBL (9773 eGenes)")+
    theme_bw()+
    theme(##strip.text=element_text(size=12),
          plot.title=element_text(hjust=0.5, size=14),
          axis.text=element_text(size=12),
          axis.title=element_text(size=12))

###
figfn <- paste(outdir, "Figure1_hist_eqtls.png", sep="")
png(figfn, width=500, height=420, res=120)
print(p)
dev.off()

### End

fn <- "/nfs/rprdata/julong/sc-atac/genetic_analysis_GTExV8/DAP-G/5_summary.outs/Whole_Blood_union_all_gene_FDR_2.txt"
x <- read.table(fn, header=T)
egene2 <- x%>%dplyr::filter(FDR_2<0.1)%>%pull(gene)

