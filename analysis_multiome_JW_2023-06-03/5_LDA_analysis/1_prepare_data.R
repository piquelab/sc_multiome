##
library(Matrix)
library(tidyverse)
library(data.table)
library(irlba)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(Signac)
library(SeuratWrappers)
library(SeuratObject)
library(GenomicRanges)



##
library(cowplot)
library(RColorBrewer)
library(ComplexHeatmap)
library(viridis)
library(circlize)
library(ggrepel)
library(ggrastr)
library(openxlsx)

rm(list=ls())


###
###


outdir <- "./1_input_data/"
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)

 


#################################
### multiome data object 
#################################

fn <- "../sc_multiome_data/1_processing/3_Clustering.outs/1_seurat.cluster.combined.mitofilt.rds"
sc <- read_rds(fn)


MCls_name <- c("0"="0_CD4Naive", "1"="1_TCM", "2"="2_NKcell", "3"="3_TEM",
    "4"="4_Bcell", "5"="5_CD8Naive", "6"="6_Monocyte", "7"="7_dnT",
    "8"="8_MAIT", "9"="9_Platelets", "10"="10_DC")

x <- sc@meta.data
x <- x%>%mutate(MCls=MCls_name[as.character(wsnn_res.0.1)],
                treat2=ifelse(treats%in%c("etOH", "water"), "control", treats))
sc <- AddMetaData(sc, x)



###
### Differential results

#### DEGs
fn <- "../sc_multiome_data/2_Differential/1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn_sampleID+new_treat.rds"
resDG <- read_rds(fn)%>%as.data.frame()%>%mutate(comb=paste(MCls, contrast, sep="_"))
DEG <- resDG%>%dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)%>%pull(gene)%>%unique()
###resDG <- resDG%>%filter(gene%in%DEG)


###
MCls <- unname((MCls_name[1:8]))
for (oneMCl in MCls){    
    ###
    sc2 <- subset(sc, subset=MCls==oneMCl) ##, features=DEG)
    ## x <- sc2[["SCT"]]$data    
    ##
    opfn <- paste(outdir, oneMCl, ".multiome.rds", sep="")
    write_rds(sc2, file=opfn)
    ###
    cat(oneMCl, "cells:", ncol(sc2), "genes:", nrow(sc2), "\n")    
}    





######################
### ATAC object 
######################

rm(list=ls())


outdir <- "./1_input_data/ATAC_outs/"
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)


fn <- "../sc_multiome_data/1_processing/5.1_reCallPeak.outs/3_scATAC.annot.mitofilt.rds"
sc <- read_rds(fn)


MCls_name <- c("0"="0_CD4Naive", "1"="1_TCM", "2"="2_NKcell", "3"="3_TEM",
    "4"="4_Bcell", "5"="5_CD8Naive", "6"="6_Monocyte", "7"="7_dnT",
    "8"="8_MAIT", "9"="9_Platelets", "10"="10_DC")

x <- sc@meta.data
x <- x%>%mutate(MCls=MCls_name[as.character(wsnn_res.0.1)],
                treat2=ifelse(treats%in%c("etOH", "water"), "control", treats))
sc <- AddMetaData(sc, x)
DefaultAssay(sc) <- "ATAC"

###
### normalization 
sc2 <- RunTFIDF(sc)
    
###



###
MCls <- unname((MCls_name[1:8]))
for (oneMCl in MCls){    
    ###
    sc0 <- subset(sc2, subset=MCls==oneMCl)
    ## x <- sc2[["SCT"]]$data    
    ##
    opfn <- paste(outdir, oneMCl, ".atac.rds", sep="")
    write_rds(sc0, file=opfn)
    ###
    cat(oneMCl, "cells:", ncol(sc2), "features:", nrow(sc2), "\n")    
}    
