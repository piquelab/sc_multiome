###
library(tidyverse)
## library(parallel)
library(data.table)
## library(purrr)
library(GenomicRanges)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(Signac)
library(SeuratWrappers)
library(SeuratObject)
library(EnsDb.Hsapiens.v75)
## library(cicero, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
## library(monocle3)
###
library(EnsDb.Hsapiens.v75)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
###
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(viridis)
library(ggrastr)
theme_set(theme_grey())


rm(list=ls())

outdir <- "./3_Clustering.outs/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


##############################################
##  combined and macs2peaks ##
##############################################
#combined.2 <- read_rds("./2_demux.outs/2_seurat.merge.SNG.combined.rds")

pbmc.combined <- read_rds("./3_Clustering.outs/1_seurat.cluster.combined.mitofilt.rds")
pbmc.combined

macs2 <- read_rds("./5.1_reCallPeak.outs/3_scATAC.annot.mitofilt.ndup.rds")
macs2

pbmc.combined[["peaks"]] <- macs2[["ATAC"]]
pbmc.combined
pbmc.combined[["peaks"]]
head(pbmc.combined@meta.data)








##############################################
##  link peaks to genes ##
##############################################
pbmc <- pbmc.combined
DefaultAssay(pbmc) <- "peaks"
pbmc

# first compute the GC content for each peak
pbmc <- RegionStats(pbmc, genome = BSgenome.Hsapiens.UCSC.hg38)

#
annotations <- GetGRangesFromEnsDb(ensdb=EnsDb.Hsapiens.v86)
#seqlevelsStyle(annotations) <- "1000"
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"
Annotation(pbmc) <- annotations

#
count <- pbmc@assays$RNA@counts
head(count)

genes <- rownames(count)
head(genes)

count.peaks <- pbmc@assays$peaks@counts
head(count.peaks)
nrow(count.peaks)

which(rownames( pbmc@assays$peaks@counts)=="chr1-992444-992774")
sp <- count.peaks.sp <- pbmc@assays$peaks@counts[70:72, ] %>% as.data.frame()
head(sp)
rownames(sp)

pbmc[["peaks"]]



head(pbmc@meta.data)
colnames(pbmc@meta.data)
##table(pbmc@meta.data$treats)








##----data for diff cell-types/treats/individuals----##
cell <- c("0"="CD4Naive", "1"="TCM", "2"="NKcell", "3"="TEM",
          "4"="Bcell", "5"="CD8Naive", "6"="Monocyte", "7"="dnT",
          "8"="MAIT","9"="Platelet", "10"="DC")


pbmc@meta.data$MCls <- cell[as.character(pbmc@meta.data$wsnn_res.0.1)]

head(pbmc@meta.data)
colnames(pbmc@meta.data)



#### subset ####
treats <- unique(pbmc@meta.data$treats)
##treats <- c("caffeine", "nicotine")
treats

for (i in treats){
    df <- subset(x = pbmc, subset = treats == i)
    df2 <- LinkPeaks(
        object = df,
        peak.assay = "peaks",
        expression.assay = "SCT",
        genes.use = genes,
        ##pvalue_cutoff = 0.1,
        distance = 100000
        )
    opfn <- paste0("./6.1_linkPeakstoGenes_2.outs/1_seurat.cluster.combined.mitofilt.peaks.linked_", i, "_dist100000.rds")
    write_rds(df2, opfn)
    assign(paste0("pbmc.",i), df2)
}
##head(pbmc.caffeine)

pbmc.caffeine <- readRDS("./6.1_linkPeakstoGenes_2.outs/1_seurat.cluster.combined.mitofilt.peaks.linked_caffeine_dist100000.rds")
pbmc.nicotine <- readRDS("./6.1_linkPeakstoGenes_2.outs/1_seurat.cluster.combined.mitofilt.peaks.linked_nicotine_dist100000.rds")
pbmc.zinc <- readRDS("./6.1_linkPeakstoGenes_2.outs/1_seurat.cluster.combined.mitofilt.peaks.linked_zinc_dist100000.rds")
pbmc.vitA <- readRDS("./6.1_linkPeakstoGenes_2.outs/1_seurat.cluster.combined.mitofilt.peaks.linked_vitA_dist100000.rds")
pbmc.vitD <- readRDS("./6.1_linkPeakstoGenes_2.outs/1_seurat.cluster.combined.mitofilt.peaks.linked_vitD_dist100000.rds")
pbmc.vitE <- readRDS("./6.1_linkPeakstoGenes_2.outs/1_seurat.cluster.combined.mitofilt.peaks.linked_vitE_dist100000.rds")
pbmc.water <- readRDS("./6.1_linkPeakstoGenes_2.outs/1_seurat.cluster.combined.mitofilt.peaks.linked_water_dist100000.rds")
pbmc.etOH <- readRDS("./6.1_linkPeakstoGenes_2.outs/1_seurat.cluster.combined.mitofilt.peaks.linked_etOH_dist100000.rds")


pbmc.caffeine
head(pbmc.caffeine)
head(Links(pbmc.caffeine))
length(Links(pbmc.caffeine))
max(Links(pbmc.caffeine)$pvalue)


# link peaks to genes
## pbmc <- LinkPeaks(
##     object = pbmc,
##     peak.assay = "peaks",
##     expression.assay = "SCT",
##     genes.use = genes,
##     pvalue_cutoff = 0.1,
##     ##distance = 100000
## )




#----links data----
dat.caffeine <- Links(pbmc.caffeine) %>% as.data.frame()
dat.nicotine <- Links(pbmc.nicotine) %>% as.data.frame()
dat.zinc <- Links(pbmc.zinc) %>% as.data.frame()
dat.vitA <- Links(pbmc.vitA) %>% as.data.frame()
dat.vitD <- Links(pbmc.vitD) %>% as.data.frame()
dat.vitE <- Links(pbmc.vitE) %>% as.data.frame()
dat.water <- Links(pbmc.water) %>% as.data.frame()
dat.etOH <- Links(pbmc.etOH) %>% as.data.frame()

head(dat.caffeine)

nrow(dat.caffeine)
nrow(dat.nicotine)
nrow(dat.zinc)
nrow(dat.vitA)
nrow(dat.vitD)
nrow(dat.vitE)
nrow(dat.water)
nrow(dat.etOH)

max(dat.caffeine$zscore)
min(dat.caffeine$zscore)
quantile(dat.caffeine$zscore)
length(unique(dat.caffeine$peak))
##nrow(dat.caffeine %>% dplyr::filter(zscore > 20))

opfn <- paste0("./6.1_linkPeakstoGenes_2.outs/links_caffeine_dist100000.rds")
write_rds(dat.caffeine, opfn)
opfn <- paste0("./6.1_linkPeakstoGenes_2.outs/links_nicotine_dist100000.rds")
write_rds(dat.nicotine, opfn)
opfn <- paste0("./6.1_linkPeakstoGenes_2.outs/links_zinc_dist100000.rds")
write_rds(dat.zinc, opfn)
opfn <- paste0("./6.1_linkPeakstoGenes_2.outs/links_vitA_dist100000.rds")
write_rds(dat.vitA, opfn)
opfn <- paste0("./6.1_linkPeakstoGenes_2.outs/links_vitD_dist100000.rds")
write_rds(dat.vitD, opfn)
opfn <- paste0("./6.1_linkPeakstoGenes_2.outs/links_vitE_dist100000.rds")
write_rds(dat.vitE, opfn)
opfn <- paste0("./6.1_linkPeakstoGenes_2.outs/links_water_dist100000.rds")
write_rds(dat.water, opfn)
opfn <- paste0("./6.1_linkPeakstoGenes_2.outs/links_etOH_dist100000.rds")
write_rds(dat.etOH, opfn)




##dat2.caffeine <- dat.caffeine %>% arrange(desc(zscore)) %>% dplyr::filter(zscore > 20)

dat2.caffeine <- dat.caffeine %>% arrange(desc(zscore)) %>% distinct(peak, .keep_all = TRUE)
dat2.nicotine <- dat.nicotine %>% arrange(desc(zscore)) %>% distinct(peak, .keep_all = TRUE)
dat2.zinc <- dat.zinc %>% arrange(desc(zscore)) %>% distinct(peak, .keep_all = TRUE)
dat2.vitA <- dat.vitA %>% arrange(desc(zscore)) %>% distinct(peak, .keep_all = TRUE)
dat2.vitD <- dat.vitD %>% arrange(desc(zscore)) %>% distinct(peak, .keep_all = TRUE)
dat2.vitE <- dat.vitE %>% arrange(desc(zscore)) %>% distinct(peak, .keep_all = TRUE)
dat2.water <- dat.water %>% arrange(desc(zscore)) %>% distinct(peak, .keep_all = TRUE)
dat2.etOH <- dat.etOH %>% arrange(desc(zscore)) %>% distinct(peak, .keep_all = TRUE)

nrow(dat2.caffeine)
nrow(dat2.nicotine)
nrow(dat2.zinc)
nrow(dat2.vitA)
nrow(dat2.vitD)
nrow(dat2.vitE)
nrow(dat2.water)
nrow(dat2.etOH)




##----comparison analysis----##
dat2.caffeine$treats <- "caffeine"
dat2.nicotine$treats <- "nicotine" 
dat2.zinc$treats <- "zinc" 
dat2.vitA$treats <- "vitA" 
dat2.vitD$treats <- "vitD" 
dat2.vitE$treats <- "vitE" 
dat2.water$treats <- "water" 
dat2.etOH$treats <- "etOH" 


head(dat2.caffeine)
nrow(dat2.caffeine)
#unique(dat2.caffeine$gene)


dat3.caffeine <- dat2.caffeine %>% distinct(gene, .keep_all = TRUE) 
dat3.nicotine <- dat2.nicotine %>% distinct(gene, .keep_all = TRUE)
dat3.zinc <- dat2.zinc %>% distinct(gene, .keep_all = TRUE)
dat3.vitA <- dat2.vitA %>% distinct(gene, .keep_all = TRUE)
dat3.vitD <- dat2.vitD %>% distinct(gene, .keep_all = TRUE)
dat3.vitE <- dat2.vitE %>% distinct(gene, .keep_all = TRUE)
dat3.water <- dat2.water %>% distinct(gene, .keep_all = TRUE)
dat3.etOH <- dat2.etOH %>% distinct(gene, .keep_all = TRUE)

head(dat3.caffeine)
head(dat3.nicotine)

nrow(dat3.caffeine)
nrow(dat3.nicotine)
nrow(dat3.zinc)
nrow(dat3.vitA)
nrow(dat3.vitD)
nrow(dat3.vitE)
nrow(dat3.water)
nrow(dat3.etOH)



end






## #----Coverage plot----
## head(pbmc@meta.data)
## colnames(pbmc@meta.data)

## table(pbmc@meta.data$predicted.id.12)


## cell <- c("0"="CD4Naive", "1"="TCM", "2"="NKcell", "3"="TEM",
##           "4"="Bcell", "5"="CD8Naive", "6"="Monocyte", "7"="dnT", "8"="MAIT",
##           "9"="Platelet", "10"="DC")

## pbmc@meta.data$MCls <- cell[as.character(pbmc@meta.data$wsnn_res.0.1)]
## colnames(pbmc@meta.data)

## #meta <- meta%>%
## #       mutate(bti=paste(wsnn_res.0.1, MCls, sep="_"))

## head(x = Idents(object = pbmc))

## Idents(pbmc) <- "MCls"
## idents.plot <- c(unique(pbmc@meta.data$MCls))
## #idents.plot <- c("B naive", "B intermediate", "B memory", "CD14 Mono", "CD16 Mono", "CD8 TEM", "CD8 Naive")
## idents.plot

## head(genes)

## p1 <- CoveragePlot(
##       object = pbmc,
##       region = "MS4A1",
##       features = "MS4A1",
##       expression.assay = "SCT",
##       idents = idents.plot,
##       extend.upstream = 500,
##       extend.downstream = 10000
##     )
## png(paste0("./6.1_linkPeakstoGenes.outs/coverage_plot_MS4A1.png"), width=1500, height=1500, pointsize=16, res=225)
## print(p1)
## dev.off()



## p1 <- CoveragePlot(
##       object = pbmc,
##       region = "MS4A1",
##       features = "MS4A1",
##       expression.assay = "SCT",
##       idents = idents.plot,
##       extend.upstream = 500,
##       extend.downstream = 10000
##     )
## p2 <- CoveragePlot(
##       object = pbmc,
##       region = "CD79A",
##       features = "CD79A",
##       expression.assay = "SCT",
##       idents = idents.plot,
##       extend.upstream = 500,
##       extend.downstream = 10000
##     )
## p3 <- CoveragePlot(
##       object = pbmc,
##       region = "CD14",
##       features = "CD14",
##       expression.assay = "SCT",
##       idents = idents.plot,
##       extend.upstream = 500,
##       extend.downstream = 10000
##     )
## p4 <- CoveragePlot(
##       object = pbmc,
##       region = "MS4A7",
##       features = "MS4A7",
##       expression.assay = "SCT",
##       idents = idents.plot,
##       extend.upstream = 500,
##       extend.downstream = 10000
##     )
## p5 <- CoveragePlot(
##       object = pbmc,
##       region = "GNLY",
##       features = "GNLY",
##       expression.assay = "SCT",
##       idents = idents.plot,
##       extend.upstream = 500,
##       extend.downstream = 10000
##     )
## p6 <- CoveragePlot(
##       object = pbmc,
##       region = "NKG7",
##       features = "NKG7",
##       expression.assay = "SCT",
##       idents = idents.plot,
##       extend.upstream = 500,
##       extend.downstream = 10000
##     )
## p7 <- CoveragePlot(
##       object = pbmc,
##       region = "CD3D",
##       features = "CD3D",
##       expression.assay = "SCT",
##       idents = idents.plot,
##       extend.upstream = 500,
##       extend.downstream = 10000
##     )
## p8 <- CoveragePlot(
##       object = pbmc,
##       region = "CD8A",
##       features = "CD8A",
##       expression.assay = "SCT",
##       idents = idents.plot,
##       extend.upstream = 500,
##       extend.downstream = 10000
##     )



## png(paste0("./6.1_linkPeakstoGenes.outs/coverage_plot_8markers.png"), width=2000, height=2000, pointsize=16, res=125)
## p <- patchwork::wrap_plots(p1, p2, p3, p4, p5, p6, p7, p8, ncol = 4)
## print(p)
## dev.off()




## ##################################################################################
## ## Finding co-accessible networks with Cicero
## ## https://satijalab.org/signac/articles/cicero.html
## ##################################################################################
## library(monocle3)

## if (!requireNamespace("remotes", quietly = TRUE))
##         install.packages("remotes")
## remotes::install_github("cole-trapnell-lab/cicero-release", ref = "monocle3")

## library(cicero)

## pbmc.cds <- as.cell_data_set(x = pbmc)

## head(pbmc.cds)


## pbmc.cicero <- make_cicero_cds(pbmc.cds, reduced_coordinates = reducedDims(pbmc.cds)$UMAP)
## head(pbmc.cicero)



## #----Find Cicero connections----
## # get the chromosome sizes from the Seurat object
## genome <- seqlengths(pbmc[["ATAC"]])
## head(genome)

## # use chromosome 1 to save some time
## # omit this step to run on the whole genome
## #genome <- genome[1]
## #genome


## # convert chromosome sizes to a dataframe
## genome.df <- data.frame("chr" = names(genome), "length" = genome)
## head(genome.df)

## # run cicero
## conns <- run_cicero(pbmc.cicero, genomic_coords = genome.df, sample_num = 100)

## head(conns)
## nrow(conns)

## opfn <- paste0("./6.1_linkPeakstoGenes.outs/2_conns.rds")
## write_rds(conns, opfn)

## head(conns %>% dplyr::filter(coaccess>0))
## nrow(conns %>% dplyr::filter(coaccess>0))


## #----peak anno----
















#----MCls----
MCls <- unique(pbmc@meta.data$MCls)
MCls

for (i in MCls){
    df <- subset(x = pbmc, subset = MCls == i)
    df2 <- LinkPeaks(
        object = df,
        peak.assay = "peaks",
        expression.assay = "SCT",
        genes.use = genes,
        ##pvalue_cutoff = 0.1,
        distance = 100000
        )
    opfn <- paste0("./6.1_linkPeakstoGenes_2.outs/1_seurat.cluster.combined.mitofilt.peaks.linked_", i, "_dist100000.rds")
    write_rds(df2, opfn)
    assign(paste0("pbmc.",i), df2)
    }
##head(pbmc.caffeine)

















inds <- unique(pbmc@meta.data$SNG.BEST.GUESS)
inds

for (i in inds){
    df <- subset(x = pbmc, subset = SNG.BEST.GUESS == i)
    df2 <- LinkPeaks(
        object = df,
        peak.assay = "peaks",
        expression.assay = "SCT",
        genes.use = genes,
        ##pvalue_cutoff = 0.1,
        distance = 100000
        )
    opfn <- paste0("./6.1_linkPeakstoGenes_2.outs/1_seurat.cluster.combined.mitofilt.peaks.linked_", i, "_dist100000.rds")
    write_rds(df2, opfn)
    assign(paste0("pbmc.",i), df2)
}
##head(pbmc.caffeine)






end
