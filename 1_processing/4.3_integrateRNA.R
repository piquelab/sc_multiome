### using cicero

## library("rhdf5")
## library("corpcor")
## library(Matrix)
## library(MASS)
## library(scales)
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
library(cicero, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(monocle3)
###
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(viridis)
theme_set(theme_grey())

outdir <- "./4.3_Integrate.outs/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)
rm(list=ls())



##############################################
### integrate annotated sc-RNA and sc-ATAC ###
##############################################

atac <- read_rds("./4.2_Integrate.outs/1.1_scATAC.cicero.rds")
cluster0 <- as.character(c(0, 1, 3, 6:10, 12:14))
atac2 <- subset(atac, subset=seurat_clusters%in%cluster0)
DefaultAssay(atac2) <- "ACTIVITY"
atac2 <- atac2%>%
   NormalizeData(assay="ACTIVITY",
   normalization.method="LogNormalize",
   scale.factor=median(atac2$nCount_ACTIVITY))


###
### read sc-RNA data
## sc <- read_rds("./3_Clustering.outs/2_seurat.RNA.rds") 
## sc <- sc%>%
##    NormalizeData()%>%
##    FindVariableFeatures(nfeatures=2000)%>%
##    ScaleData()%>%
##    RunPCA(npcs=100)%>%
##    RunUMAP(dims=1:50)
### add annotation into data
## meta <- read_rds("./2_RNA.outs/3.1_scHOLD.RNAseq.annot.rds")@meta.data
## x <- meta%>%dplyr::select(barcode, predicted.celltype.l1.score, predicted.celltype.l1, predicted.celltype.l2.score, predicted.celltype.l2)
## meta <- sc@meta.data
## meta$barcode <- rownames(meta)
## metaNew <- meta%>%left_join(x, by="barcode")
## rownames(metaNew) <- metaNew$barcode
## sc <- AddMetaData(sc, metaNew)


###
### transfer cell type annotation from sc-RNAseq to sc-ATAC data
sc <- read_rds("./3_Clustering.outs/2_seurat.RNA.rds") 
anchors <- FindTransferAnchors(reference=sc, query=atac2,
   features=VariableFeatures(object=sc),
   reference.assay="RNA", query.assay ="ACTIVITY", reduction="cca")
opfn <- "./4.3_Integrate.outs/2.0_anchors.rds"
write_rds(anchors, opfn)

### Using predicted.celltype.l1
pred <- TransferData(anchorset=anchors, refdata=sc$predicted.celltype.l1,
    weight.reduction=atac[["lsi"]], dims = 2:50)
pred <- pred%>%mutate(barcode=rownames(pred))
opfn <- "./4.3_Integrate.outs/2.1_pred.rds"
write_rds(pred, file=opfn)

### Using predicted.celltype.l2
pred <- TransferData(anchorset=anchors, refdata=sc$predicted.celltype.l2,
    weight.reduction=atac[["lsi"]], dims = 2:50)
pred <- pred%>%mutate(barcode=rownames(pred))
opfn <- "./4.3_Integrate.outs/2.2_pred.rds"
write_rds(pred, file=opfn)


###
### Add annotation into atac data
## atac <- read_rds("./4.2_Integrate.outs/1.1_scATAC.cicero.rds")
x <- atac@meta.data
###
pred1 <- read_rds("./4.2_Integrate.outs/2.1_pred.rds")%>%
    dplyr::select(predicted.id, barcode)%>%
    dplyr::rename("predicted.celltype.l1"="predicted.id")
###
pred2 <- read_rds("./4.2_Integrate.outs/2.2_pred.rds")%>%
    dplyr::select(predicted.id, barcode)%>%
    dplyr::rename("predicted.celltype.l2"="predicted.id")

xNew <- x%>%
   left_join(pred1, by=c("NEW_BARCODE"="barcode"))%>%
   left_join(pred2, by=c("NEW_BARCODE"="barcode"))
rownames(xNew) <- xNew$NEW_BARCODE

atac <- AddMetaData(atac, xNew)
opfn <- "./4.2_Integrate.outs/3_scATAC.annot.rds"
write_rds(atac, opfn)

### define cell type based on cluster
atac <- read_rds("./4.2_Integrate.outs/3_scATAC.annot.rds")
cell <- c("0"="Tcell", "1"="Tcell", "2"="Monocyte", "3"="Tcell",
   "4"="Bcell", "5"="NKcell", "6"="Tcell", "7"="Tcell", "8"="Tcell",
   "9"="Tcell", "10"="Tcell", "11"="DC", "12"="Tcell",
   "13"="Tcell", "14"="Tcell")
atac$MCls <- cell[as.character(atac$seurat_clusters)]
DefaultAssay(atac) <- "ATAC"
opfn <- "./4.2_Integrate.outs/3_scATAC.annot.rds"
write_rds(atac, opfn)


##################################
### summary annottaion results ###
##################################

fn <- "./4.2_Integrate.outs/3_scATAC.annot.rds"
atac <- read_rds(fn)

p0 <- DimPlot(atac, reduction="umap.atac", label=T, raster=F)+
   theme_bw()## +
   ## ## ## guides(col=guide_legend(override.aes=list(size=2),ncol=3))+
   ## theme(legend.title=element_blank(),
   ##       legend.key.size=grid::unit(0.8,"lines"))
figfn <- "./4.2_Integrate.outs/Figure1.0_atac.cluster.png"
png(figfn, width=500, height=500, res=120)
print(p0)
dev.off()

###
p1 <- DimPlot(atac, reduction="umap.atac", group.by="predicted.celltype.l1",
   label=T,  raster=F, repel=T)+
   theme_bw()+
   theme(plot.title=element_blank(),
         legend.title=element_blank(),
         ## legend.text=element_text(size=8),
         legend.key.size=grid::unit(1,"lines"))
figfn <- "./4.2_Integrate.outs/Figure1.1_atac.pred1.png"
png(figfn, width=500, height=450, res=120)
print(p1)
dev.off()

###
p2 <- DimPlot(atac, reduction="umap.atac", group.by="predicted.celltype.l2",
   label=T, label.size=2.5, raster=F, repel=T)+
   theme_bw()+
   theme(plot.title=element_blank(),
         legend.title=element_blank(),
         legend.text=element_text(size=8),
         legend.key.size=grid::unit(0.8,"lines"))

figfn <- "./4.2_Integrate.outs/Figure1.2_atac.pred2.png"
png(figfn, width=650, height=400, res=120)
print(p2)
dev.off()


###
### heatmap
fn <- "./4.2_Integrate.outs/3_scATAC.annot.rds"
atac <- read_rds(fn)
meta <- atac@meta.data


df1 <- meta%>%
   group_by(predicted.celltype.l1, seurat_clusters)%>%
   summarize(Freq=n(),.groups="drop")%>%
   group_by(seurat_clusters)%>%
   mutate(Perc=Freq/sum(Freq))

p1 <- ggplot(df1, aes(x=seurat_clusters, y=predicted.celltype.l1, fill=Perc))+
   geom_tile()+
   scale_fill_gradient("Fraction of cells",
      low="#ffffc8", high="#7d0025", na.value=NA)+
   xlab("Cluster")+ylab("predicted.celltype.l1")+
   theme_bw()+theme(axis.text.x=element_text(hjust=0.5, size=10))

###
figfn <- "./4.2_Integrate.outs/Figure2.1_atac.heatmap.png"
png(figfn, width=800, height=600, res=120)
print(p1)
dev.off()


###
df2 <- meta%>%
   group_by(predicted.celltype.l2, seurat_clusters)%>%
   summarize(Freq=n(),.groups="drop")%>%
   group_by(seurat_clusters)%>%
   mutate(Perc=Freq/sum(Freq))

p2 <- ggplot(df2, aes(x=seurat_clusters, y=predicted.celltype.l2, fill=Perc))+
   geom_tile()+
   scale_fill_gradient("Fraction of cells",
      low="#ffffc8", high="#7d0025", na.value=NA)+
   xlab("Cluster")+ylab("predicted.celltype.l2")+
   theme_bw()+theme(axis.text.x=element_text(hjust=0.5))

###
figfn <- "./4.2_Integrate.outs/Figure2.2_atac.heatmap.png"
png(figfn, width=600, height=600, res=120)
print(p2)
dev.off()




###
### umap color by cell types
p0 <- DimPlot(atac, reduction="umap.atac", group.by="MCls", label=T, raster=F)+
   theme_bw()+
   scale_color_manual(values=c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
                               "NKcell"="#aa4b56", "Tcell"="#ffaa00",
                               "DC"="#828282"))+
   theme(legend.position="none",
         plot.title=element_blank())

   ## ## ## guides(col=guide_legend(override.aes=list(size=2),ncol=3))+
   ## theme(legend.title=element_blank(),
   ##       legend.key.size=grid::unit(0.8,"lines"))
figfn <- "./4.2_Integrate.outs/Figure3.0_umap.atac.MCls.png"
png(figfn, width=420, height=500, res=120)
print(p0)
dev.off()
