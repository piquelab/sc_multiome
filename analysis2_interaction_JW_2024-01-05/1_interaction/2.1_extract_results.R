###
library(tidyverse)
library(DESeq2)
library(biobroom)
library(parallel)

rm(list=ls())

outdir <- "./2_inter_ATAC.outs/dds_results/"
if ( ! file.exists(outdir)) dir.create(outdir, showWarnings=T, recursive=T) 

args <- commandArgs(trailingOnly=T)
if ( length(args)>0){
   ##
   ii <- as.numeric(args[1])
}else{
   ii <- 1
}   

contrast_df <- read.table("./2_inter_ATAC.outs/contrast_list.txt", header=T)

contrast_ii <- as.character(contrast_df[ii,])
contrast_nn <- contrast_df$con1[ii]

dds_full <- read_rds("./2_inter_ATAC.outs/2_full_ATAC.dds.rds")
res2 <- results(dds_full, contrast=contrast_ii)%>%tidy()
res2$condition <- contrast_nn

## output
opfn <- paste(outdir, ii, "_results.rds", sep="")
write_rds(res2, file=opfn)

### END



