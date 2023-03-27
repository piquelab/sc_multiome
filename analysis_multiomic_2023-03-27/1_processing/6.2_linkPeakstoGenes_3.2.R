###
library(tidyverse)
## library(parallel)
library(data.table)
## library(purrr)
##library(GenomicRanges)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(Signac)
library(SeuratWrappers)
library(SeuratObject)
##library(EnsDb.Hsapiens.v75)
## library(cicero, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
## library(monocle3)
###
##library(EnsDb.Hsapiens.v75)
##library(EnsDb.Hsapiens.v86)
##library(BSgenome.Hsapiens.UCSC.hg38)
###
##library(ggplot2)
##library(cowplot)
##library(RColorBrewer)
##library(viridis)
##library(ggrastr)
##theme_set(theme_grey())


rm(list=ls())

#outdir <- "./3_Clustering.outs/"
#if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


##############################################
##  sys_args ##
##############################################
### passing argument
## args=commandArgs(trailingOnly=T)

## #condition <- args[1]

## if (length(args)>0){
##     condition <- args[1]
## ###dir <- args[2]
## }else{
##     condition <- "caffeine.Bcell.AL-002"
## ###dir <- "SCAIP_results"
## }

## condition <- gsub("\r", "", condition)

## print(condition)








##############################################
##  df/corr files ##
##############################################
dat.ind.B <- read_rds("./6.2_linkPeakstoGenes_2.5.outs/df_ind_Bcell.rds") %>% as.data.frame
dat.ind.NK <- read_rds("./6.2_linkPeakstoGenes_2.5.outs/df_ind_NKcell.rds") %>% as.data.frame
dat.ind.mono <- read_rds("./6.2_linkPeakstoGenes_2.5.outs/df_ind_Monocyte.rds") %>% as.data.frame
dat.ind.TCM <- read_rds("./6.2_linkPeakstoGenes_2.5.outs/df_ind_TCM.rds") %>% as.data.frame
dat.ind.TEM <- read_rds("./6.2_linkPeakstoGenes_2.5.outs/df_ind_TEM.rds") %>% as.data.frame
dat.ind.CD4 <- read_rds("./6.2_linkPeakstoGenes_2.5.outs/df_ind_CD4Naive.rds") %>% as.data.frame
dat.ind.CD8 <- read_rds("./6.2_linkPeakstoGenes_2.5.outs/df_ind_CD8Naive.rds") %>% as.data.frame
dat.ind.DC <- read_rds("./6.2_linkPeakstoGenes_2.5.outs/df_ind_DC.rds") %>% as.data.frame
dat.ind.MAIT <- read_rds("./6.2_linkPeakstoGenes_2.5.outs/df_ind_MAIT.rds") %>% as.data.frame
dat.ind.platelet <- read_rds("./6.2_linkPeakstoGenes_2.5.outs/df_ind_Platelet.rds") %>% as.data.frame
dat.ind.dnT <- read_rds("./6.2_linkPeakstoGenes_2.5.outs/df_ind_dnT.rds") %>% as.data.frame

head(dat.ind.B)

nrow(dat.ind.B)
nrow(dat.ind.NK)
nrow(dat.ind.mono)
nrow(dat.ind.TCM)
nrow(dat.ind.TEM)
nrow(dat.ind.CD4)
nrow(dat.ind.CD8)
nrow(dat.ind.DC)
nrow(dat.ind.MAIT)
nrow(dat.ind.platelet)
nrow(dat.ind.dnT)


dat.ind <- rbind(dat.ind.B, dat.ind.NK, dat.ind.mono, dat.ind.TCM, dat.ind.TEM, dat.ind.CD4, dat.ind.CD8, dat.ind.DC, dat.ind.MAIT, dat.ind.platelet, dat.ind.dnT)

head(dat.ind)
nrow(dat.ind)

table(dat.ind$celltype, dat.ind$treatment)

dat.ind.cor.sig <- dat.ind %>% dplyr::filter(pvalue<0.05)
dat.ind.zcor.sig <- dat.ind %>% dplyr::filter(z_pvalue<0.05)





## treats ##
dat.treats.B <- read_rds("./6.2_linkPeakstoGenes_2.5.outs/df_treats_Bcell.rds") %>% as.data.frame
dat.treats.NK <- read_rds("./6.2_linkPeakstoGenes_2.5.outs/df_treats_NKcell.rds") %>% as.data.frame
dat.treats.mono <- read_rds("./6.2_linkPeakstoGenes_2.5.outs/df_treats_Monocyte.rds") %>% as.data.frame
dat.treats.TCM <- read_rds("./6.2_linkPeakstoGenes_2.5.outs/df_treats_TCM.rds") %>% as.data.frame
dat.treats.TEM <- read_rds("./6.2_linkPeakstoGenes_2.5.outs/df_treats_TEM.rds") %>% as.data.frame
dat.treats.CD4 <- read_rds("./6.2_linkPeakstoGenes_2.5.outs/df_treats_CD4Naive.rds") %>% as.data.frame
dat.treats.CD8 <- read_rds("./6.2_linkPeakstoGenes_2.5.outs/df_treats_CD8Naive.rds") %>% as.data.frame
dat.treats.DC <- read_rds("./6.2_linkPeakstoGenes_2.5.outs/df_treats_DC.rds") %>% as.data.frame
dat.treats.MAIT <- read_rds("./6.2_linkPeakstoGenes_2.5.outs/df_treats_MAIT.rds") %>% as.data.frame
dat.treats.platelet <- read_rds("./6.2_linkPeakstoGenes_2.5.outs/df_treats_Platelet.rds") %>% as.data.frame
dat.treats.dnT <- read_rds("./6.2_linkPeakstoGenes_2.5.outs/df_treats_dnT.rds") %>% as.data.frame


nrow(dat.treats.B)
nrow(dat.treats.NK)
nrow(dat.treats.mono)
nrow(dat.treats.TCM)
nrow(dat.treats.TEM)
nrow(dat.treats.CD4)
nrow(dat.treats.CD8)
nrow(dat.treats.DC)
nrow(dat.treats.MAIT)
nrow(dat.treats.platelet)
nrow(dat.treats.dnT)


dat.treats <- rbind(dat.treats.B, dat.treats.NK, dat.treats.mono, dat.treats.TCM, dat.treats.TEM, dat.treats.CD4, dat.treats.CD8, dat.treats.DC, dat.treats.MAIT, dat.treats.platelet, dat.treats.dnT)

head(dat.treats)
nrow(dat.treats)

dat.treats.cor.sig <- dat.treats %>% dplyr::filter(pvalue<0.05)
dat.treats.zcor.sig <- dat.treats %>% dplyr::filter(z_pvalue<0.05)


## box-whisker plot ##
figfn <- "./6.2_linkPeakstoGenes_3.outs/Figure1.2_cor_boxplot_only_common.png"
png(figfn, width=600, height=600, res=125)
boxplot(as.numeric(dat.ind$cor),
        as.numeric(dat.treats$cor),
        main = "score_correlation",
        names=c("individuals", "treatments"))
dev.off()

t.test(as.numeric(dat.ind$cor),as.numeric(dat.treats$cor))




figfn <- "./6.2_linkPeakstoGenes_3.outs/Figure1.2_zcor_boxplot_only_common.png"
png(figfn, width=600, height=600, res=125)
boxplot(as.numeric(dat.ind$z_cor),
        as.numeric(dat.treats$z_cor),
        main = "zscore_correlation",
        names=c("individuals", "treatments"))
dev.off()

t.test(as.numeric(dat.ind$z_cor),as.numeric(dat.treats$z_cor))




figfn <- "./6.2_linkPeakstoGenes_3.outs/Figure1.2_cor.sig_boxplot_only_common.png"
png(figfn, width=600, height=600, res=125)
##print(boxplot.cor.sig)
boxplot(as.numeric(dat.ind.cor.sig$cor),
        as.numeric(dat.treats.cor.sig$cor),
        main = "sig_score_correlation",
        names=c("individuals", "treatments"))
dev.off()

t.test(as.numeric(dat.ind.cor.sig$cor),as.numeric(dat.treats.cor.sig$cor))




figfn <- "./6.2_linkPeakstoGenes_3.outs/Figure1.2_zcor.sig_boxplot_only_common.png"
png(figfn, width=600, height=600, res=125)
boxplot(as.numeric(dat.ind.zcor.sig$z_cor),
        as.numeric(dat.treats.zcor.sig$z_cor),
        main = "sig_zscore_correlation",
        names=c("individuals", "treatments"))
##print(boxplot.zcor.sig)
dev.off()

t.test(as.numeric(dat.ind.zcor.sig$z_cor),as.numeric(dat.treats.zcor.sig$z_cor))


end


















