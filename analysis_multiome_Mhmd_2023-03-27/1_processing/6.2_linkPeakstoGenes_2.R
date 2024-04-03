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
args=commandArgs(trailingOnly=T)

#condition <- args[1]

if (length(args)>0){
    condition <- args[1]
###dir <- args[2]
}else{
    condition <- "caffeine.Bcell.AL-002"
###dir <- "SCAIP_results"
}

condition <- gsub("\r", "", condition)

print(condition)


#pbmc <- read_rds("./6.2_linkPeakstoGenes.outs/pbmc_comb.rds")
pbmc <- read_rds("./6.1_linkPeakstoGenes_2.outs/pbmc_comb.rds")

pbmc@meta.data$comb <- gsub("-", "_", pbmc@meta.data$comb)

##pbmc
##head(pbmc@meta.data)

#colnames(pbmc@meta.data)
#head(pbmc@meta.data$comb)

#pbmc[["peaks"]]
#peak.count <- pbmc@assays$peaks@counts
#head(peak.count)
#peak.count[chr1-629069-629348, ]

##for (i in unique(pbmc@meta.data$comb)){
##    print(i)
##    }




df <- subset(x = pbmc, subset = comb == condition)
count <- df@assays$RNA@counts
##head(count)
genes <- rownames(count)
##head(genes)



## pbmc.1 <- subset(x = pbmc, subset = comb == "caffeine.NKcell.AL-045")
## #pbmc.1
## count.1 <- pbmc.1@assays$RNA@counts
## ##head(count.1)
## genes.1 <- rownames(count.1)
## ##head(genes.1)
## ##length(genes.1)
## ##count.1.peaks <- pbmc.1@assays$peaks@counts
## ##head(count.1.peaks)
## ##peaks.1 <- rownames(count.1.peaks)
## ##head(peaks.1)
## ##length(peaks.1)




## pbmc.2 <- subset(x = pbmc, subset = comb == "nicotine.NKcell.AL-045")
## count.2 <- pbmc.2@assays$RNA@counts
## ##head(count.2)
## genes.2 <- rownames(count.2)
## ##head(genes.2)
## ##length(genes.2)
## ##count.2.peaks <- pbmc.2@assays$peaks@counts
## ##head(count.2.peaks)
## ##peaks.2 <- rownames(count.2.peaks)
## ##head(peaks.2)
## ##length(peaks.2)


# link peaks to genes
## pbmc <- LinkPeaks(
##     object = df,
##     peak.assay = "peaks",
##     expression.assay = "SCT",
##     genes.use = genes,
##     pvalue_cutoff = 1,
##     ##distance = 100000
## )

pbmc <- LinkPeaks(
    object = df,
    peak.assay = "peaks",
    expression.assay = "SCT",
    genes.use = genes,
    pvalue_cutoff = 1,
    score_cutoff = 0,
    min.cells = 0
    ##distance = 100000
)


## pbmc.1 <- LinkPeaks(
##     object = pbmc.1,
##     peak.assay = "peaks",
##     expression.assay = "SCT",
##     genes.use = genes.1,
##     pvalue_cutoff = 1,
##     score_cutoff = 0,
##     min.cells = 0
##     ##distance = 100000
## )

## pbmc.2 <- LinkPeaks(
##     object = pbmc.2,
##     peak.assay = "peaks",
##     expression.assay = "SCT",
##     genes.use = genes.2,
##     pvalue_cutoff = 1,
##     score_cutoff = 0,
##     min.cells = 0
##     ##distance = 100000
## )




## head(Links(pbmc))
## length(Links(pbmc))
## max(Links(pbmc)$pvalue)

##all files saved in 6.2_linkPeakstoGenes.outs folder are 100000 dist and 0.05 pvalue - if I try different options later will make resp folders
##opfn <- paste0("./6.2_linkPeakstoGenes_2.outs/linkPeaks.", condition, "dist100000.rds")
#opfn <- paste0("./6.2_linkPeakstoGenes_2.outs/linkPeaks_", condition, "_dist100000.rds")
#write_rds(pbmc, opfn)




## #----links data----
dat <- Links(pbmc) %>% as.data.frame()
#head(dat)
#nrow(dat)
#max(dat$zscore)
#min(dat$zscore)
#quantile(dat$zscore)
#length(unique(dat$peak))
#nrow(dat %>% dplyr::filter(zscore > 20))

opfn <- paste0("./6.2_linkPeakstoGenes_2.2.outs/links_", condition, "_pval1_distdefault_score0.rds")
write_rds(dat, opfn)


#dat2 <- dat %>% arrange(desc(zscore)) %>% dplyr::filter(zscore > 20)
#head(dat2)
#nrow(dat2)
#unique(dat2$gene)




## #----testing score0----
## dat1 <- Links(pbmc.1) %>% as.data.frame()
## dat2 <- Links(pbmc.2) %>% as.data.frame()


## opfn <- paste0("./6.2_linkPeakstoGenes_2.outs/links_", "caffeine.NKcell.AL-045", "_pval1_distdefault_score0.rds")
## write_rds(dat1, opfn)

## opfn <- paste0("./6.2_linkPeakstoGenes_2.outs/links_", "nicotine.NKcell.AL-045", "_pval1_distdefault_score0.rds")
## write_rds(dat2, opfn)


## #dat1 <- readRDS(paste0("./6.2_linkPeakstoGenes_2.outs/links_", "caffeine.NKcell.AL-045", "_pval1_distdefault_score0.rds"))
## #dat2 <- readRDS(paste0("./6.2_linkPeakstoGenes_2.outs/links_", "caffeine.NKcell.AL-045", "_pval1_distdefault_score0.rds"))

## head(dat1)
## head(dat2)


## dat1 %>% dplyr::filter(peak=="chr1-629069-629348")
## dat2 %>% dplyr::filter(peak=="chr1-629069-629348")
## dat2 %>% dplyr::filter(gene=="FAM87B")

## common_peaks <- intersect(dat1$peak, dat2$peak)

## length(common_peaks)
## head(common_peaks)

## nrow(dat1)
## nrow(dat2)

## dat1_ori <- readRDS(paste0("./6.2_linkPeakstoGenes_2.outs/links_", "caffeine.NKcell.AL_045", "_pval1_distdefault.rds"))
## dat2_ori <- readRDS(paste0("./6.2_linkPeakstoGenes_2.outs/links_", "nicotine.NKcell.AL_045", "_pval1_distdefault.rds"))

#dat1_ori <- readRDS(paste0("./6.2_linkPeakstoGenes_2.2.outs/links_", "caffeine.Bcell.AL_002", "_pval1_distdefault_score0.rds"))
#dat2_ori <- readRDS(paste0("./6.2_linkPeakstoGenes_2.2.outs/links_", "nicotine.Bcell.AL_002", "_pval1_distdefault_score0.rds"))

#head(dat1_ori)
#head(dat2_ori)
#nrow(dat1_ori)
#nrow(dat2_ori)





## #########################################################
## ############## dat -  2.2.outs ##########################
## #########################################################
## dat1 <- readRDS(paste0("./6.2_linkPeakstoGenes_2.2.outs/links_", "caffeine.NKcell.AL_045", "_pval1_distdefault_score0.rds"))

## dat1 <- dat1 %>% mutate(peak.gene = paste0(peak, ".", gene))

## dat2 <- readRDS(paste0("./6.2_linkPeakstoGenes_2.2.outs/links_", "nicotine.NKcell.AL_045", "_pval1_distdefault_score0.rds"))

## dat2 <- dat2 %>% mutate(peak.gene = paste0(peak, ".", gene))

## head(dat1)
## head(dat2)

## nrow(dat1)
## nrow(dat2)


## dat <- full_join(dat1, dat2, by="peak.gene")
## dat[is.na(dat)] <- 0


## scorr <- as.numeric(cor.test(dat$score.x, dat$score.y, method = "spearman")$estimate)
## scorr
## pval <- as.numeric(cor.test(dat$score.x, dat$score.y, method = "spearman")$p.value)
## pval


## scorr_z <- as.numeric(cor.test(dat$zscore.x, dat$zscore.y, method = "spearman")$estimate)
## scorr_z
## pval_z <- as.numeric(cor.test(dat$zscore.x, dat$zscore.y, method = "spearman")$p.value)
## pval_z



## head(dat)
## nrow(dat)

