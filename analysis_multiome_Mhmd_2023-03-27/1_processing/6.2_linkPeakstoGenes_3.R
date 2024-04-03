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

## pbmc <- read_rds("./6.1_linkPeakstoGenes_2.outs/pbmc_comb.rds")
## head(pbmc@meta.data)
## colnames(pbmc@meta.data)

## pbmc[["RNA"]]
## pbmc[["ATAC"]]
## pbmc[["peaks"]]


## meta <- pbmc@meta.data %>% as.data.frame()
## head(meta)
## colnames(meta)


## meta.1 <- meta %>% dplyr::filter(comb=="caffeine.Bcell.AL-002")
## meta.1$nFeature_ATAC
## nrow(meta.1)


## meta.2 <- meta %>% dplyr::filter(comb=="caffeine.NKcell.AL-045")
## nrow(meta.1)
## nrow(meta.2)




#dat1 <- read_rds("./6.2_linkPeakstoGenes_2.outs/links_caffeine.Bcell.AL_002_pval1_distdefault.rds") %>% as.data.frame
#head(dat1)
#nrow(dat1)
#length(unique(dat1$peak))


dat1 <- read_rds("./6.2_linkPeakstoGenes_2.outs/links_caffeine.NKcell.AL_045_pval1_distdefault.rds") %>% as.data.frame
dat2 <- read_rds("./6.2_linkPeakstoGenes_2.outs/links_nicotine.NKcell.AL_045_pval1_distdefault.rds") %>% as.data.frame
dat3 <- read_rds("./6.2_linkPeakstoGenes_2.outs/links_vitA.NKcell.AL_045_pval1_distdefault.rds") %>% as.data.frame
dat4 <- read_rds("./6.2_linkPeakstoGenes_2.outs/links_vitD.NKcell.AL_045_pval1_distdefault.rds") %>% as.data.frame
dat5 <- read_rds("./6.2_linkPeakstoGenes_2.outs/links_vitE.NKcell.AL_045_pval1_distdefault.rds") %>% as.data.frame
dat6 <- read_rds("./6.2_linkPeakstoGenes_2.outs/links_zinc.NKcell.AL_045_pval1_distdefault.rds") %>% as.data.frame

dat1[1:6,]
dat2[1:6,]

nrow(dat1)
nrow(dat2)

nrow(dat3)
nrow(dat4)
nrow(dat5)
nrow(dat6)


length(unique(dat2$peak))

dat1.2 <- dat1 %>% dplyr::filter(peak %in% peak1)
head(dat1.2)
nrow(dat1.2)


dat2.2 <- dat2 %>% dplyr::filter(peak %in% peak2)
head(dat2.2)
nrow(dat2.2)

dat1.peak <- unique(dat1$peak)
dat2.peak <- unique(dat2$peak)
dat.peak <- intersect(dat1.peak, dat2.peak)
dat1.peak.set <- setdiff(dat1.peak, dat.peak)
dat2.peak.set <- setdiff(dat2.peak, dat.peak)


length(dat1.peak)
length(dat2.peak)
length(dat.peak)
length(dat1.peak.set)
length(dat2.peak.set)




end




#############
## checking
#############
pbmc <- read_rds("./6.1_linkPeakstoGenes_2.outs/pbmc_comb.rds")


####################################
pbmc.1 <- subset(x = pbmc, subset = comb == "caffeine.NKcell.AL-045")

pbmc.1@meta.data

colnames(pbmc.1@meta.data)

head(pbmc.1@meta.data)

max(pbmc.1@meta.data$percent.mt.RNA)
min(pbmc.1@meta.data$percent.mt.RNA)

quantile(pbmc.1@meta.data$percent.mt.RNA, 0.99)

##--------------------------------
count.1 <- pbmc.1@assays$peaks@counts
#head(count.1)

anno.1 <- data.frame(rn=rownames(count.1), rnz=rowSums(count.1))
head(anno.1)

nrow(anno.1)
#anno.1["chr1-633631-634104",]

annoSel.1 <- anno.1%>%dplyr::filter(rnz>0)
#head(annoSel.1)
nrow(annoSel.1)
#annoSel.1["chr1-633631-634104", ]

peak1 <- unique(annoSel.1$rn)
head(peak1)
length(peak1)


##--------------------------------
count.1.rna <- pbmc.1@assays$SCT@counts
head(count.1.rna)

anno.1.rna <- data.frame(rn=rownames(count.1.rna), rnz=rowSums(count.1.rna))
head(anno.1.rna)

nrow(anno.1.rna)
#anno.1.rna["ISG15",]

annoSel.1.rna <- anno.1.rna%>%dplyr::filter(rnz>0)

head(annoSel.1.rna)

nrow(annoSel.1.rna)





####################################
pbmc.2 <- subset(x = pbmc, subset = comb == "nicotine.NKcell.AL-045")


##--------------------------------
count.2 <- pbmc.2@assays$peaks@counts
#head(count.2)

anno.2 <- data.frame(rn=rownames(count.2), rnz=rowSums(count.2))
head(anno.2)
#nrow(anno.2)
#anno.2["chr1-633631-634104",]

annoSel.2 <- anno.2%>%dplyr::filter(rnz>0)
#head(annoSel.2)
nrow(annoSel.2)
#annoSel.2["chr1-633631-634104", ]

peak2 <- unique(annoSel.2$rn)
head(peak2)
length(peak2)


##--------------------------------
count.2.rna <- pbmc.2@assays$SCT@counts
#head(count.2.rna)

anno.2.rna <- data.frame(rn=rownames(count.2.rna), rnz=rowSums(count.2.rna))
head(anno.2.rna)
#nrow(anno.2.rna)
#anno.2.rna["ISG15",]

annoSel.2.rna <- anno.2.rna%>%dplyr::filter(rnz>0)
#head(annoSel.2.rna)
nrow(annoSel.2.rna)




#which(rownames(pbmc.1@assays$peaks@counts)=="chr1-629069-629348")
#sp <- count.peaks.sp <- pbmc.1@assays$peaks@counts[22:24, ] %>% as.data.frame()
#head(sp)
#rownames(sp)


dat.test <- readRDS(paste0("./6.2_linkPeakstoGenes_2.2.outs/links_", "caffeine.Bcell.AL_002", "_pval1_distdefault_score0.rds")) 

head(dat.test)

write.table(dat.test, paste0("./6.2_linkPeakstoGenes_2.2.outs/links_", "caffeine.Bcell.AL_002", "_pval1_distdefault_score0.csv"), sep="\t")

end


















####################
##
####################

#pbmc <- read_rds("./6.2_linkPeakstoGenes_2.outs/linkPeaks_zinc.TEM.AL_105_dist100000.rds")
#pbmc.2 <- read_rds("./6.2_linkPeakstoGenes_2.outs/linkPeaks_zinc.TEM.AL_058_dist100000.rds")
#pbmc.3 <- read_rds("./6.2_linkPeakstoGenes_2.outs/linkPeaks_caffeine.TEM.AL_105_dist100000.rds")
#pbmc.4 <- read_rds("./6.2_linkPeakstoGenes_2.outs/linkPeaks_caffeine.TEM.AL_058_dist100000.rds")

pbmc <- read_rds("./6.2_linkPeakstoGenes.outs/linkPeaks_zinc.TEM.AL_105_dist100000.rds")
pbmc.2 <- read_rds("./6.2_linkPeakstoGenes.outs/linkPeaks_zinc.TEM.AL_058_dist100000.rds")
pbmc.3 <- read_rds("./6.2_linkPeakstoGenes.outs/linkPeaks_caffeine.TEM.AL_105_dist100000.rds")
pbmc.4 <- read_rds("./6.2_linkPeakstoGenes.outs/linkPeaks_caffeine.TEM.AL_058_dist100000.rds")




##pbmc@meta.data$comb <- gsub("-", "_", pbmc@meta.data$comb)

pbmc
head(pbmc@meta.data)
#colnames(pbmc@meta.data)
#head(pbmc@meta.data$comb)



#----links data----
dat <- Links(pbmc) %>% as.data.frame()
head(dat)
nrow(dat)

max(dat$zscore)
min(dat$zscore)
quantile(dat$zscore)


length(unique(dat$peak))

nrow(dat %>% dplyr::filter(zscore > 20))


opfn <- paste0("./6.1_linkPeakstoGenes.outs/links_dist100000.rds")
write_rds(dat, opfn)


dat2 <- dat %>% arrange(desc(zscore)) %>% dplyr::filter(zscore > 20)
head(dat2)
nrow(dat2)
unique(dat2$gene)



#---------------------------
dat.4 <- Links(pbmc.4) %>% as.data.frame()
head(dat.4)
nrow(dat.4)




##for (i in unique(pbmc@meta.data$comb)){
##    print(i)
##    }




## df <- subset(x = pbmc, subset = comb == condition)
## count <- df@assays$RNA@counts
## ##head(count)
## genes <- rownames(count)
## ## head(genes)


## # link peaks to genes
## pbmc <- LinkPeaks(
##     object = df,
##     peak.assay = "peaks",
##     expression.assay = "SCT",
##     genes.use = genes,
##     ##pvalue_cutoff = 0.1,
##     distance = 100000
## )


## ## head(Links(pbmc))
## ## length(Links(pbmc))
## ## max(Links(pbmc)$pvalue)

## ##all files saved in 6.2_linkPeakstoGenes.outs folder are 100000 dist and 0.05 pvalue - if I try different options later will make resp folders
## ##opfn <- paste0("./6.2_linkPeakstoGenes_2.outs/linkPeaks.", condition, "dist100000.rds")
## opfn <- paste0("./6.2_linkPeakstoGenes_2.outs/linkPeaks_", condition, "_dist100000.rds")
## write_rds(pbmc, opfn)








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









