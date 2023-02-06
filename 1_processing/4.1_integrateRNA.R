#module load R/test_4.0.3

library("rhdf5")
library("corpcor")
library(Matrix)
library(MASS)
library(scales)
library(tidyverse)
library(parallel)
library(data.table)
library(future)
library(purrr)
library(furrr)
##
library(GenomicRanges)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(Signac)
library(EnsDb.Hsapiens.v75)
###
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(viridis)
theme_set(theme_grey())

rm(list=ls())

###########################
### infer gene activity ###
###########################

outdir <- "./4_Integrate.outs/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)

###
atac <- read_rds("./3_Clustering.outs/2_seurat.cluster.rds")
gene.activities <- GeneActivity(atac)

atac[["ACTIVITY"]] <- CreateAssayObject(counts=gene.activities)
atac2 <- NormalizeData(object=atac,
                       assay="ACTIVITY",
                       normalization.method="LogNormalize")
opfn <- "./4_Integrate.outs/1_seurat.activity.rds"
write_rds(atac2, file=opfn)


## summary plots ###
DefaultAssay(atac2) <- "ACTIVITY"

x0 <- c("MS4A1", "CD79A", "MS4A7", "CD14", "GNLY", "NKG7", "CD3D", "CD8A")

fig1 <- FeaturePlot(object=atac2, 
                    features=x0, ncol=4, raster=F)&
        #scale_color_gradient("",low="lightgrey",high="blue")+
        theme_bw()+
           theme(legend.title=element_blank(),
                 legend.key.size=grid::unit(0.5,"lines"),
                 plot.title=element_text(size=12, hjust=0.5))
png("./4_Integrate.outs/Figure1.1_MCls.feature.png", width=950, height=500, res=100)
print(fig1)
dev.off()   
                     

#######################################
#### Integration ATAC with RNA data ###
#######################################

###rna data
sc <- read_rds("/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/5_IdenCelltype_output/1_SCAIP.spliced.rds")
sc <- subset(sc, subset=BATCH=="SCAIP6")

sc2 <- read_rds("/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/5_IdenCelltype_output/4_SCAIP.MCls.Harmony.rds")
sc2 <- subset(sc2, subset=BATCH=="SCAIP6")
x <- sc2@meta.data

sc3 <- read_rds("./3_Clustering.outs/3_seurat.rna.rds")
x3 <- sc3@meta.data%>%dplyr::select("NEW_BARCODE","predicted.celltype.l1")

x <- x%>%left_join(x3,by="NEW_BARCODE")
rownames(x) <- x$NEW_BARCODE

sc <- AddMetaData(sc, x) ##use for following analysis

### normalized datat and dimensional reduction
sc <- sc%>%
      NormalizeData()%>%
      FindVariableFeatures(nfeatures=3000)%>%
      ScaleData()%>%
      RunPCA(npcs=100)%>%
      RunUMAP(dims=1:50)

atac <- read_rds("./4_Integrate.outs/1_seurat.activity.rds")
atac <- atac%>%
        NormalizeData()%>%
        ScaleData(features=rownames(atac))

### Find anchors between RNA and atac 
anchors <- FindTransferAnchors(
        reference=sc, query=atac, features=VariableFeatures(object=sc),
        reference.assay="RNA", query.assay ="ACTIVITY", reduction="cca")

opfn <- "./4_Integrate.outs/anchors.rds"
write_rds(anchors, file=opfn)

### transfer
pred1 <- TransferData(anchorset=anchors, refdata=sc$MCls,
                     weight.reduction=atac[["lsi"]], dims = 2:50)
pred2 <- TransferData(anchorset=anchors, refdata=sc$predicted.celltype.l1,
                     weight.reduction=atac[["lsi"]], dims = 2:50)

pred1 <- pred1%>%mutate(barcode=rownames(pred1))%>%dplyr::select("barcode","predicted.id")%>%dplyr::rename("predicted.id1"="predicted.id")
pred2 <- pred2%>%mutate(barcode=rownames(pred2))%>%dplyr::select("barcode","predicted.id")%>%dplyr::rename("predicted.id2"="predicted.id")

x <- atac@meta.data
xNew <- x%>%left_join(pred1,by="barcode")%>%left_join(pred2, by="barcode")

rownames(xNew) <- xNew$barcode
atac2 <- AddMetaData(atac, metadata=xNew)

opfn <- "./4_Integrate.outs/2_seurat.atac.annot.rds"
write_rds(atac2, file=opfn)



###############
### summary ###
###############
atac <- read_rds("./4_Integrate.outs/2_seurat.atac.annot.rds")

### UMAP annotated by predicted.id1
p1 <- DimPlot(object=atac, reduction="umap", label=T, raster=F)+NoLegend()+
      theme_bw()+
      theme(legend.position="none")

p2 <- DimPlot(atac, reduction="umap", group.by="predicted.id1",
      label=T, label.size=2.5, raster=F, repel=T)+
      NoLegend()+
      theme_bw()+
      theme(legend.position="none", plot.title=element_blank())
figfn <- "./4_Integrate.outs/Figure2.1_annot.UMAP.png"
png(figfn, width=700, height=400, res=120)
print(plot_grid(p1, p2, ncol=2))
dev.off()

###
p1 <- DimPlot(object=atac, reduction="umap", label=T, raster=F)+NoLegend()+
      theme_bw()+
      theme(legend.position="none")

p2 <- DimPlot(atac, reduction="umap", group.by="predicted.id2",
      label=T, label.size=2.5, raster=F, repel=T)+
      NoLegend()+
      theme_bw()+
      theme(legend.position="none", plot.title=element_blank())
figfn <- "./4_Integrate.outs/Figure2.2_annot.UMAP.png"
png(figfn, width=700, height=400, res=120)
print(plot_grid(p1, p2, ncol=2))
dev.off()






