##
library(tidyverse)
library(data.table)
library(annotables)
options(scipen=16)
library(openxlsx)

###
### The script for candidate gene information that are used for extracting genotypes 


fn <- "../1_candidate_genes.xlsx"
res <- read.xlsx(fn)

### get gene infor
geneSel <- unique(res$Gene)
anno <- grch38%>%filter(ensgene%in%geneSel, chr%in%as.character(1:22), !duplicated(ensgene))

### Start-End infor 
anno2 <- anno%>%mutate(start=ifelse(strand==1, start, end))%>%dplyr::select(ensgene, symbol, chr, s0=start)
anno2 <- anno2%>%mutate(start=s0-2e+06, end=s0+2e+06)%>%dplyr::select(-s0)

write.table(anno2, "gene_infor.txt", quote=F, row.names=F, col.names=F)

