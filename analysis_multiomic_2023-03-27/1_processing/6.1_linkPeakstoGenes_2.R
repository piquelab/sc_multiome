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



#### all combinations ####
treats <- unique(pbmc@meta.data$treats)
treats
MCls <- unique(pbmc@meta.data$MCls)
MCls
inds <- unique(pbmc@meta.data$SNG.BEST.GUESS)
inds

df <- subset(x = pbmc, subset = treats == "caffeine")
df2 <- subset(x = df, subset = MCls == "NKcell")
df3 <- subset(x = df2, subset = SNG.BEST.GUESS == "AL-045")

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













## # link peaks to genes
## pbmc <- LinkPeaks(
##     object = pbmc,
##     peak.assay = "peaks",
##     expression.assay = "SCT",
##     genes.use = genes,
##     pvalue_cutoff = 0.1,
##     ##distance = 100000
## )


## head(Links(pbmc))
## length(Links(pbmc))
## max(Links(pbmc)$pvalue)

## opfn <- paste0("./6.1_linkPeakstoGenes.outs/1_seurat.cluster.combined.mitofilt.peaks.linked_dist100000.rds")
## write_rds(pbmc, opfn)



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
