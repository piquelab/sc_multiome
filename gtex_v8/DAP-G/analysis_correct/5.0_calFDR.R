##
library(Matrix)
library(tidyverse)
library(data.table)
library(qvalue)
## library(clusterProfiler)
## library(org.Hs.eg.db)
## library(annotables)
## library(ggplot2)
## library(cowplot)
## library(RColorBrewer)

rm(list=ls())


###
### calculate FDR

args=commandArgs(trailingOnly=T)
if ( length(args)>0){
    ##
    geneFile <- args[1]
}else{
    geneFile <- "splitGene000"
}    


outdir <- "./5_summary.outs/Whole_Blood/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)




#########################################
### Calculate local FDR for each geen ###
#########################################

fn <- paste("../geneList/", geneFile, sep="")
geneList <- read.table(fn)$V1
dap <- NULL
for (ens in geneList){    
###
   dapfn <- paste("./dap-g_outs/Whole_Blood/", ens, ".model.out", sep="")
   if ( file.exists(dapfn)&file.size(dapfn)>0){
      ###    
      resTMP <- read.table(dapfn, fill=T, row.names=NULL, header=F)
      res2 <- resTMP%>%filter(V3==0)%>%mutate(gene=ens)%>%dplyr::select(gene, lfdr=V2)
      dap <- rbind(dap, res2)
  }    
}


### 
opfn2 <- paste(outdir, "calFDR_", geneFile, ".txt", sep="")
write.table(dap, opfn2, row.names=F, col.names=F, quote=F, sep="\t")
### End
