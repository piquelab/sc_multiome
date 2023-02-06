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


pbmc <- read_rds("./6.2_linkPeakstoGenes.outs/pbmc_comb.rds")

pbmc@meta.data$comb <- gsub("-", "_", pbmc@meta.data$comb)

##pbmc
##head(pbmc@meta.data)

#colnames(pbmc@meta.data)
#head(pbmc@meta.data$comb)



##for (i in unique(pbmc@meta.data$comb)){
##    print(i)
##    }



df <- subset(x = pbmc, subset = comb == condition)


count <- df@assays$RNA@counts
##head(count)

genes <- rownames(count)
## head(genes)


# link peaks to genes
pbmc <- LinkPeaks(
    object = df,
    peak.assay = "peaks",
    expression.assay = "SCT",
    genes.use = genes,
    pvalue_cutoff = 1,
    ##distance = 100000
)


## head(Links(pbmc))
## length(Links(pbmc))
## max(Links(pbmc)$pvalue)

##all files saved in 6.2_linkPeakstoGenes.outs folder are 100000 dist and 0.05 pvalue - if I try different options later will make resp folders
##opfn <- paste0("./6.2_linkPeakstoGenes_2.outs/linkPeaks.", condition, "dist100000.rds")
opfn <- paste0("./6.2_linkPeakstoGenes_2.outs/linkPeaks_", condition, "_dist100000.rds")
write_rds(pbmc, opfn)








## ##############################################
## ##  combined and macs2peaks ##
## ##############################################
## #combined.2 <- read_rds("./2_demux.outs/2_seurat.merge.SNG.combined.rds")

## pbmc.combined <- read_rds("./3_Clustering.outs/1_seurat.cluster.combined.mitofilt.rds")
## pbmc.combined

## macs2 <- read_rds("./5.1_reCallPeak.outs/3_scATAC.annot.mitofilt.ndup.rds")
## macs2

## pbmc.combined[["peaks"]] <- macs2[["ATAC"]]
## pbmc.combined
## pbmc.combined[["peaks"]]
## head(pbmc.combined@meta.data)








## ##############################################
## ##  link peaks to genes ##
## ##############################################
## pbmc <- pbmc.combined
## DefaultAssay(pbmc) <- "peaks"
## pbmc

## # first compute the GC content for each peak
## pbmc <- RegionStats(pbmc, genome = BSgenome.Hsapiens.UCSC.hg38)

## #
## annotations <- GetGRangesFromEnsDb(ensdb=EnsDb.Hsapiens.v86)
## #seqlevelsStyle(annotations) <- "1000"
## seqlevelsStyle(annotations) <- "UCSC"
## genome(annotations) <- "hg38"
## Annotation(pbmc) <- annotations

## #
## count <- pbmc@assays$RNA@counts
## head(count)

## genes <- rownames(count)
## head(genes)

## count.peaks <- pbmc@assays$peaks@counts
## head(count.peaks)
## nrow(count.peaks)

## which(rownames( pbmc@assays$peaks@counts)=="chr1-992444-992774")
## sp <- count.peaks.sp <- pbmc@assays$peaks@counts[70:72, ] %>% as.data.frame()
## head(sp)
## rownames(sp)

## pbmc[["peaks"]]



## head(pbmc@meta.data)
## colnames(pbmc@meta.data)
## ##table(pbmc@meta.data$treats)








## #######################################################
## ##----data for diff cell-types/treats/individuals----##
## #######################################################
## cell <- c("0"="CD4Naive", "1"="TCM", "2"="NKcell", "3"="TEM",
##           "4"="Bcell", "5"="CD8Naive", "6"="Monocyte", "7"="dnT",
##           "8"="MAIT","9"="Platelet", "10"="DC")


## pbmc@meta.data$MCls <- cell[as.character(pbmc@meta.data$wsnn_res.0.1)]


## pbmc@meta.data$comb <- paste0(pbmc@meta.data$treats, ".", pbmc@meta.data$MCls, ".", pbmc@meta.data$SNG.BEST.GUESS)

## head(pbmc@meta.data)
## colnames(pbmc@meta.data)

## opfn <- paste0("./6.1_linkPeakstoGenes_2.outs/pbmc_comb.rds")
## write_rds(pbmc, opfn) 






## #### all combinations ####
## treats <- unique(pbmc@meta.data$treats)
## treats
## MCls <- unique(pbmc@meta.data$MCls)
## MCls
## inds <- unique(pbmc@meta.data$SNG.BEST.GUESS)
## inds

## df <- subset(x = pbmc, subset = treats == "caffeine")
## df2 <- subset(x = df, subset = MCls == "NKcell")
## df3 <- subset(x = df2, subset = SNG.BEST.GUESS == "AL-045")

## #### subset ####
## treats <- unique(pbmc@meta.data$treats)
## ##treats <- c("caffeine", "nicotine")
## treats

## for (i in treats){
##     df <- subset(x = pbmc, subset = treats == i)
##     df2 <- LinkPeaks(
##         object = df,
##         peak.assay = "peaks",
##         expression.assay = "SCT",
##         genes.use = genes,
##         ##pvalue_cutoff = 0.1,
##         distance = 100000
##         )
##     opfn <- paste0("./6.1_linkPeakstoGenes_2.outs/1_seurat.cluster.combined.mitofilt.peaks.linked_", i, "_dist100000.rds")
##     write_rds(df2, opfn)
##     assign(paste0("pbmc.",i), df2)
## }
## ##head(pbmc.caffeine)




## MCls <- unique(pbmc@meta.data$MCls)
## MCls

## for (i in MCls){
##     df <- subset(x = pbmc, subset = MCls == i)
##     df2 <- LinkPeaks(
##         object = df,
##         peak.assay = "peaks",
##         expression.assay = "SCT",
##         genes.use = genes,
##         ##pvalue_cutoff = 0.1,
##         distance = 100000
##         )
##     opfn <- paste0("./6.1_linkPeakstoGenes_2.outs/1_seurat.cluster.combined.mitofilt.peaks.linked_", i, "_dist100000.rds")
##     write_rds(df2, opfn)
##     assign(paste0("pbmc.",i), df2)
##     }
## ##head(pbmc.caffeine)




## inds <- unique(pbmc@meta.data$SNG.BEST.GUESS)
## inds

## for (i in inds){
##     df <- subset(x = pbmc, subset = SNG.BEST.GUESS == i)
##     df2 <- LinkPeaks(
##         object = df,
##         peak.assay = "peaks",
##         expression.assay = "SCT",
##         genes.use = genes,
##         ##pvalue_cutoff = 0.1,
##         distance = 100000
##         )
##     opfn <- paste0("./6.1_linkPeakstoGenes_2.outs/1_seurat.cluster.combined.mitofilt.peaks.linked_", i, "_dist100000.rds")
##     write_rds(df2, opfn)
##     assign(paste0("pbmc.",i), df2)
## }
## ##head(pbmc.caffeine)






## end








## #----links data----
## dat <- Links(pbmc) %>% as.data.frame()

## head(dat)

## nrow(dat)

## max(dat$zscore)
## min(dat$zscore)
## quantile(dat$zscore)


## length(unique(dat$peak))

## nrow(dat %>% dplyr::filter(zscore > 20))


## opfn <- paste0("./6.1_linkPeakstoGenes.outs/links_dist100000.rds")
## write_rds(dat, opfn)


## dat2 <- dat %>% arrange(desc(zscore)) %>% dplyr::filter(zscore > 20)
## head(dat2)
## nrow(dat2)
## unique(dat2$gene)



