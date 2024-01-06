###
library(Matrix)
library(tidyverse)
library(data.table)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(Signac)
library(SeuratWrappers)
library(SeuratObject)
library(ArchR)
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

outdir <- "./ArchR_output/"
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)




#################################################
### creat a gene activity matrix using signac ###
#################################################


fn <- "../sc_multiome_data/1_processing/3_Clustering.outs/1_seurat.cluster.combined.mitofilt.rds"
sc <- read_rds(fn)
meta <- sc@meta.data
meta2 <- meta%>%dplyr::select(gex_barcode, EXP, wsnn_res.0.1)%>%
   mutate(NEW_BARCODE=paste0(EXP, "#", gex_barcode))  

## DefaultAssay(sc) <- "ATAC"
## gene.activities <- GeneActivity(sc)
## ann <- Annotation(sc)


###
basefolder <- "/wsu/home/groups/piquelab/scGxE/counts_cellranger_arc/"

expNames <- dir(basefolder, "^scGxE_1.*")

folders <- paste0(basefolder, expNames)

ind <- dir.exists(folders)

###
### framments files used for analysis
folders <- folders[ind]
expNames <- expNames[ind]


###
## addArchRThreads(threads=1)

addArchRGenome("hg38")

###
### 
inputFiles <- paste0(folders, "/outs/atac_fragments.tsv.gz", sep="")
names(inputFiles) <- expNames

###
### creating arrow Files
for (i in 2:length(expNames)){
   ###
   ArrowFiles <- createArrowFiles(inputFiles=inputFiles[i], sampleNames=expNames[i],
      filterTSS=4, filterFrags=1000, addTileMat=T, addGeneScoreMat=T)
}                               


###
### Create ArchRProject
ArrowFiles <- paste(expNames, ".arrow", sep="")
proj <- ArchRProject(ArrowFiles=ArrowFiles, outputDirectory="ArchR_output", copyArrows=TRUE)

### subsetting
idx <- BiocGenerics::which(proj$cellNames%in%meta2$NEW_BARCODE)
cellSel <- proj$cellNames[idx]
proj2 <- proj[cellSel,]

### add cluster in cellColData
Cluster <- as.character(meta2$wsnn_res.0.1)
names(Cluster) <- meta2$NEW_BARCODE


Cluster2 <- Cluster[cellSel]
proj2$cluster <- Cluster2

### Dimensionality reduction and clustering 
proj2 <- addIterativeLSI(ArchRProj=proj2, useMatrix="TileMatrix", name="IterativeLSI")

proj2 <- addClusters(input=proj2, reducedDims="IterativeLSI")

proj2 <- addUMAP(ArchRProj=proj2, reducedDims="IterativeLSI")

### Imputation
proj2 <- addImputeWeights(proj2)


###
### output
proj2 <- saveArchRProject(ArchRProj=proj2)




###
proj3 <- loadArchRProject(path="ArchR_output")
### We have already added  gene score
### proj2 <- addGeneScoreMatrix(proj2, force=T)

x <- getCellColData(proj3)
x <- x%>%as.data.frame()%>%
    rownames_to_column(var="NEW_BARCODE_ArchR")%>%
    dplyr::rename("cluster_signac"="cluster", "cluster_ArchR"="Clusters")

opfn <- "1_meta.data.rds"
write_rds(x, opfn)


###
mat <- getMatrixFromProject(ArchRProj=proj3)
opfn <- "1.2_genescore.rds"
write_rds(mat, file=opfn)
