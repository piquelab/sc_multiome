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
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(viridis)
library(ggrastr)
theme_set(theme_grey())


rm(list=ls())

outdir <- "./3_Clustering.outs/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)

## ##################
## ### clean data ###
## ##################
## ### subset SNG cells
## ## atac <- read_rds("./1_Merge.outs/1_seurat.merge.rds")
## ## demux <- read_rds("/nfs/rprdata/julong/sc-atac/demux.2021-01-23/outs.summary/1_demuxlet.rds")
## ## meta <- atac@meta.data
## ## meta <- meta%>%left_join(demux, by=c("barcode"="NEW_BARCODE"))
## ## meta0 <- meta%>%dplyr::filter(DROPLET.TYPE=="SNG")
## ## rownames(meta0) <- meta0$barcode
## ## ##
## ## atac2 <- subset(atac, cells=meta0$barcode)
## ## atac2@meta.data <- meta0
## ## opfn <- "./3_Clustering.outs/1_seurat.SNG.rds"
## ## write_rds(atac2, opfn)





## ###########################
## ### Clustering analysis ###
## ###########################

## rm(list=ls())

## atac <- read_rds("./2_demux.outs/2_seurat.merge.SNG.rds")

## colnames(atac@meta.data)

## atac2 <- atac%>%
##    RunTFIDF()%>%
##    FindTopFeatures(min.cutoff="q0")%>%
##    RunSVD(n=100)%>%
##    RunUMAP(reduction="lsi", dims=2:50,
##            reduction.name="umap.atac", reduction.key="atacUMAP")
## ##
## atac3 <- atac2%>%
##    FindNeighbors(reduction="lsi", dims=2:50)%>%
##    FindClusters(resolution=0.15, verbose=FALSE, algorithm=3)       

## write_rds(atac3, file="./3_Clustering.outs/1_seurat.cluster.rds")

## fig1 <- DepthCor(atac2)+
##    theme(plot.title=element_text(size=10),
##          plot.subtitle=element_text(size=10))
## figfn <- "./3_Clustering.outs/Figure1.0_depthcor.png"
## png(figfn, width=600, height=400, res=120)
## print(fig1)
## dev.off()


## ###############
## ### summary ###
## ###############

## atac2 <- read_rds("./3_Clustering.outs/1_seurat.cluster.rds")

## #------------------------------------------
## fig1 <- DimPlot(object=atac2, label=TRUE, raster=F)+
##     ## NoLegend()+
##     theme_bw()## +
##     ## theme(legend.position="none")
## figfn <- "./3_Clustering.outs/Figure1.1_cluster.png"
## png(figfn, width=500, height=500, res=120)
## print(fig1)
## dev.off()


## #------------------------------------------
## ### data for umap plot
## umap <- as.data.frame(atac2[["umap.atac"]]@cell.embeddings)
## x <- atac2@meta.data
## df2 <- data.frame(UMAP_1=umap[,1],
##    UMAP_2=umap[,2],
##    seurat_clusters=x$seurat_clusters,
##    treat=gsub(".*-ATAC-|_.*", "", x$NEW_BARCODE),
##    Batch=gsub("-.*", "", x$NEW_BARCODE))
## ## alpha <- c("LPS"="a", "LPS-DEX"="b", "PHA"="c", "PHA-DEX"="d", "CTRL"="e")
## ## atac2$treat.alpha <- alpha[atac2$treat]

## ### 2,
## fig2 <- ggplot(df2, aes(x=UMAP_1, y=UMAP_2, colour=factor(seurat_clusters)))+
##    rasterise(geom_point(size=0.1),dpi=300)+
##    facet_wrap(~factor(Batch),ncol=2)+
##    ## scale_colour_manual(values=col0,
##    ##     guide=guide_legend(override.aes=list(size=2)))+
##    guides(col=guide_legend(override.aes=list(size=2),ncol=1))+
##    theme_bw()+
##    theme(legend.title=element_blank(),
##          legend.background=element_blank(),
##          legend.box.background=element_blank(),
##          legend.key.size=grid::unit(1,"lines"),
##          strip.text=element_text(size=12))
## ###
## png("./3_Clustering.outs/Figure1.2_BATCH.png", width=700, height=450, res=120)
## print(fig2)
## dev.off()

## #------------------------------------------
## ### 3,
## alpha <- c("LPS"="a","LPS-DEX"="b", "PHA"="c", "PHA-DEX"="d", "CTRL"="e")
## df2$treat.alpha <- alpha[df2$treat]
## Labtreat <- c("a"="LPS", "b"="LPS+DEX", "c"="PHA", "d"="PHA+DEX", "e"="CTRL")

## fig3 <- ggplot(df2, aes(x=UMAP_1, y=UMAP_2, colour=factor(seurat_clusters)))+
##    rasterise(geom_point(size=0.1), dpi=300)+
##    facet_wrap(~factor(treat.alpha), ncol=3, labeller=as_labeller(Labtreat))+
##    guides(col=guide_legend(override.aes=list(size=2),ncol=3))+
##    theme_bw()+
##    theme(strip.text=element_text(size=12),
##          legend.title=element_blank(),
##          legend.position=c(0.85,0.25),
##          legend.background=element_blank(),
##          legend.box.background=element_blank(),
##          legend.key.size=grid::unit(1,"lines"))

## figfn <- "./3_Clustering.outs/Figure1.3_treat.png"
## png(figfn, width=700, height=650, res=120)
## print(fig3)
## dev.off()


## ## atac2$Batch <- gsub("-.*", "", colnames(atac2))
## ## fig2 <- DimPlot(object=atac2, raster=F)+
## ##    facet_grid(~Batch)+ 
## ##    theme_bw()+
## ##    theme(legend.position="none",
## ##          plot.title=element_blank())
## ## figfn <- "./3_Clustering.outs/Figure1.2_BATCH.png"
## ## png(figfn, width=650, height=350, res=120)
## ## print(fig2)
## ## dev.off()

## ## col1 <- c("e"="#828282", 
## ##            "a"="#fb9a99", "b"="#e31a1c",
## ##            "c"="#a6cee3", "d"="#1f78b4")
## ## atac2$treat <- gsub(".*-ATAC-|_.*", "", colnames(atac2))
## ## alpha <- c("LPS"="a", "LPS-DEX"="b", "PHA"="c", "PHA-DEX"="d", "CTRL"="e")
## ## atac2$treat.alpha <- alpha[atac2$treat]

## ## fig3 <- DimPlot(object=atac2, group.by="treat.alpha", raster=F)+NoLegend()+
## ##         facet_wrap(~treat.alpha, ncol=2, 
## ##                    labeller=as_labeller(c("a"="LPS", "b"="LPS-DEX", "c"="PHA", "d"="PHA-DEX", "e"="CTRL")))+
## ##         scale_colour_manual(values=col1)+
## ##         theme_bw()+
## ##         theme(legend.position="none",
## ##               plot.title=element_blank())
## ## figfn <- "./3_clustering.outs/Figure2.3_treat.png"
## ## png(figfn, width=650, height=650, res=120)
## ## print(fig3)
## ## dev.off()









## ############
## #### RNA ###
## ############
## rm(list=ls())
## ### annotation cell type using reference data
## ref <- LoadH5Seurat("../pbmc_multimodal.h5seurat")
## sc <- read_rds("/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/5_IdenCelltype_output/1_SCAIP.spliced.rds")
## sc2 <- subset(sc, subset=BATCH=="SCAIP6")
## sc2 <- SCTransform(sc2, verbose=FALSE)


## anchors <- FindTransferAnchors(reference=ref, query=sc2,
##    normalization.method="SCT",
##    reference.reduction="spca",
##    dims=1:50)


## ###
## sc3 <- MapQuery(anchorset=anchors, query=sc2, reference=ref,
##    refdata = list(celltype.l1="celltype.l1",
##    celltype.l2="celltype.l2",
##    predicted_ADT="ADT"),
##    reference.reduction="spca",
##    reduction.model="wnn.umap")
## x <- sc3@meta.data


## ###dimensional reduction 
## sc2 <- subset(sc, subset=BATCH=="SCAIP6")
## sc4.0 <- sc2%>%
##    NormalizeData()%>%
##    FindVariableFeatures(selection.method="vst", nfeatures=2000)%>%
##    ScaleData()%>%
##    RunPCA(npcs=100)%>%
##    RunUMAP(dims=1:50)
## sc4 <- sc4.0%>%
##    FindNeighbors(dims=1:50)%>%
##    FindClusters(resolution=0.15)
## x <- sc4@meta.data
## ###
## metaNew <- cbind(sc3@meta.data,x[,39:40])
## sc4 <- AddMetaData(sc4,metaNew)
## ###
## opfn <- "./3_Clustering.outs/2_seurat.RNA.rds"
## write_rds(sc4, file=opfn)










## ###############
## ### summary ###
## ###############

## fn <- "./3_Clustering.outs/2_seurat.RNA.rds"
## sc <- read_rds(fn)

## #------------------------------------------
## p1 <- DimPlot(object=sc, reduction="umap", label=T, raster=F)+
##    theme_bw()+
##    ## guides(col=guide_legend(override.aes=list(size=2),ncol=3))+
##    theme(legend.title=element_blank(),
##          legend.key.size=grid::unit(1,"lines"))    
## figfn <- "./3_Clustering.outs/Figure2.1_RNA.cluster.png" 
## png(figfn, width=420, height=400, res=120)
## print(p1)
## dev.off()

## #------------------------------------------
## ###
## p2 <- DimPlot(sc, reduction="umap", group.by="predicted.celltype.l1",
##    label=T,  raster=F, repel=T)+
##    theme_bw()+
##    theme(plot.title=element_blank(),
##          legend.title=element_blank(),
##          ## legend.text=element_text(size=8),
##          legend.key.size=grid::unit(1,"lines"))
## figfn <- "./3_Clustering.outs/Figure2.2_RNA.pred1.png" 
## png(figfn, width=450, height=400, res=120)
## print(p2)
## dev.off()

## #------------------------------------------
## ###
## p3 <- DimPlot(sc, reduction="umap", group.by="predicted.celltype.l2",
##    label=T, label.size=2.5, raster=F, repel=T)+
##    theme_bw()+
##    theme(plot.title=element_blank(),
##          legend.title=element_blank(),
##          legend.text=element_text(size=8),
##          legend.key.size=grid::unit(0.8,"lines"))

## figfn <- "./3_Clustering.outs/Figure2.3_RNA.pred2.png"
## png(figfn, width=650, height=400, res=120)
## print(p3)
## dev.off()







## #------------------------------------------
## #------------------------------------------
## ###
## ### heatmap
## fn <- "./3_Clustering.outs/2_seurat.RNA.rds"
## sc <- read_rds(fn)
## meta <- sc@meta.data

## #------------------------------------------
## df1 <- meta%>%
##    group_by(predicted.celltype.l1, seurat_clusters)%>%
##    summarize(Freq=n(),.groups="drop")%>%
##    group_by(seurat_clusters)%>%
##    mutate(Perc=Freq/sum(Freq))

## p1 <- ggplot(df1, aes(x=seurat_clusters, y=predicted.celltype.l1, fill=Perc))+
##    geom_tile()+
##    scale_fill_gradient("Fraction of cells",
##       low="#ffffc8", high="#7d0025", na.value=NA)+
##    xlab("Cluster")+ylab("predicted.celltype.l1")+
##    theme_bw()+theme(axis.text.x=element_text(hjust=0.5, size=10))

## ###
## figfn <- "./3_Clustering.outs/Figure3.1_heatmap.png"
## png(figfn, width=800, height=600, res=120)
## print(p1)
## dev.off()

## #------------------------------------------
## ###
## df2 <- meta%>%
##    group_by(predicted.celltype.l2, seurat_clusters)%>%
##    summarize(Freq=n(),.groups="drop")%>%
##    group_by(seurat_clusters)%>%
##    mutate(Perc=Freq/sum(Freq))

## p2 <- ggplot(df2, aes(x=seurat_clusters, y=predicted.celltype.l2, fill=Perc))+
##    geom_tile()+
##    scale_fill_gradient("Fraction of cells",
##        low="#ffffc8", high="#7d0025", na.value=NA)+
##    xlab("Cluster")+ylab("predicted.celltype.l2")+
##    theme_bw()+theme(axis.text.x=element_text(hjust=0.5))

## ###
## figfn <- "./3_Clustering.outs/Figure3.2_heatmap.png"
## png(figfn, width=600, height=600, res=120)
## print(p2)
## dev.off()








## ########################################################
## ### UMAP of cell type-related marker gene expression ###
## ########################################################

## rm(list=ls())

## fn <- "./3_Clustering.outs/2_seurat.RNA.rds"
## sc <- read_rds(fn)

## x0 <- c("MS4A1", "CD79A", "MS4A7", "CD14", "GNLY", "NKG7", "CD3D", "CD8A")
##   ## "IL7R", "CCR7", "S100A4", "CD8A", "TNFRSF18", "ID3")
## MCls <- rep(c("B cells", "Monocytes", "NK cells", "T cells"), each=2)
## anno <- data.frame(MCls=MCls, symbol=x0)

## #------------------------------------------
## #### plot
## figs_ls <- lapply(1:nrow(anno), function(i){
##    oneMCl <- anno$MCls[i]
##    gene <- anno$symbol[i]
##    fig0 <- FeaturePlot(sc, features=gene)+
##       scale_color_gradient(gene, low="lightgrey", high="blue")+
##       ggtitle(oneMCl)+
##       theme_bw()+
##       theme(axis.text=element_text(size=6),
##             axis.title=element_text(size=8),
##             legend.title=element_text(size=6),
##             legend.key.size=grid::unit(0.4,"lines"),
##             legend.text=element_text(size=6),
##             plot.title=element_text(size=8, hjust=0.5))
##    fig0
## })

## ###
## figfn <- "./3_Clustering.outs/Figure3.3_RNA.feature.png"
## png(figfn, width=900, height=450, res=120)
## plot_grid(figs_ls[[1]], figs_ls[[3]], figs_ls[[5]], figs_ls[[7]],
##           figs_ls[[2]], figs_ls[[4]], figs_ls[[6]], figs_ls[[8]],
##    align="hv", nrow=2, ncol=4, byrow=T)
## dev.off()




## #------------------------------------------
## #### previous UMAP and cell type annotation
## ## sc <- read_rds("/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/5_IdenCelltype_output/4_SCAIP.MCls.Harmony.rds")
## ## sc2 <- subset(sc, subset=BATCH=="SCAIP6")
## ## p1 <- DimPlot(object=sc2, reduction="umap", label=T, raster=F)+NoLegend()+
## ##       theme_bw()+
## ##       theme(legend.position="none")

## ## p2 <- DimPlot(sc2, reduction="umap", group.by="MCls",
## ##               label=T, label.size=2.5, raster=F, repel=T)+
## ##       NoLegend()+
## ##       theme_bw()+
## ##       theme(legend.position="none", plot.title=element_blank())
## ## figfn <- "./3_Clustering.outs/Figure4.1_annot.UMAP.png"
## ## png(figfn, width=700, height=400, res=120)
## ## print(plot_grid(p1, p2, ncol=2))
## ## dev.off()




















##############################################
##  combined ##
##############################################

###########################
### Clustering analysis ###
###########################
rm(list=ls())

combined.2 <- read_rds("./2_demux.outs/2_seurat.merge.SNG.combined.rds")

combined.2
colnames(combined.2@meta.data) 


## # combined after clustering - ideally not needed but incase
## #----------------------------------------------------------
## combined.3 <- read_rds("../1_processing/3_Clustering.outs/1_seurat.cluster.combined.rds") 

## combined.3

## colnames(combined.3@meta.data)

## x <- combined.2@meta.data
## x <- x %>% dplyr::select(nucleosome_signal, nucleosome_percentile, TSS.enrichment, TSS.percentile)
## #head(x)
## combined.3 <- AddMetaData(combined.3, x)
## colnames(combined.3@meta.data)

## #combined.3 <- subset(x = combined.3, subset = percent.mt.RNA < 20)



#---------------------
combined <- combined.2
combined

DefaultAssay(combined) <- "RNA"
combined[["percent.mt.RNA"]] <- PercentageFeatureSet(combined, pattern = "^MT-")
max(combined@meta.data$percent.mt.RNA)

DefaultAssay(combined) <- "ATAC"
combined[["percent.mt.ATAC"]] <- PercentageFeatureSet(combined, pattern = "^MT-")
max(combined@meta.data$percent.mt.ATAC)

combined <- subset(x = combined, subset = percent.mt.RNA < 20)

combined
colnames(combined@meta.data)


#################################################
# Gene expression data processing
# We can normalize the gene expression data using SCTransform, and reduce the dimensionality using PCA.
#################################################
DefaultAssay(combined) <- "RNA"
combined <- SCTransform(combined)
combined <- RunPCA(combined)
## png("./PCA_noclustering.png")
## DimPlot(combined, reduction = "pca", label = TRUE)
## dev.off()


## fig1 <- DepthCor(atac2)+
##    theme(plot.title=element_text(size=10),
##          plot.subtitle=element_text(size=10))
## figfn <- "./3_Clustering.outs/Figure1.0_depthcor.png"
## png(figfn, width=600, height=400, res=120)
## print(fig1)
## dev.off()



#################################################
# DNA accessibility data processing
# Here we process the DNA accessibility assay the same way we would process a scATAC-seq dataset, by performing latent semantic indexing (LSI).
#################################################
DefaultAssay(combined) <- "ATAC"
combined <- RunTFIDF(combined)

combined <- FindTopFeatures(combined, min.cutoff = "q0")

combined <- RunSVD(combined, n=100)

colnames(combined@meta.data)

## png("./LSI_noclustering.png")
## DimPlot(combined, reduction = "lsi", label = TRUE)
## dev.off()

## atac2 <- atac%>%
##    RunTFIDF()%>%
##    FindTopFeatures(min.cutoff="q0")%>%
##    RunSVD(n=100)%>%
##    RunUMAP(reduction="lsi", dims=2:50,
##            reduction.name="umap.atac", reduction.key="atacUMAP")
## ##
## atac3 <- atac2%>%
##    FindNeighbors(reduction="lsi", dims=2:50)%>%
##    FindClusters(resolution=0.15, verbose=FALSE, algorithm=3)       

## write_rds(atac3, file="./3_Clustering.outs/1_seurat.cluster.rds")


## fig1 <- DepthCor(atac2)+
##    theme(plot.title=element_text(size=10),
##          plot.subtitle=element_text(size=10))
## figfn <- "./3_Clustering.outs/Figure1.0_depthcor.png"
## png(figfn, width=600, height=400, res=120)
## print(fig1)
## dev.off()



##################################################
# Annotating cell types
# To annotate cell types in the dataset we can transfer cell labels from an existing PBMC reference dataset using tools in the Seurat package. See the Seurat reference mapping vignette for more information.
# Weâll use an annotated PBMC reference dataset from Hao et al. (2020), available for download here: https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat
# Note that the SeuratDisk package is required to load the reference dataset. Installation instructions for SeuratDisk can be found here.To annotate cell types in the dataset we can transfer cell labels from an existing PBMC reference dataset using tools in the Seurat package. See the Seurat reference mapping vignette for more information.
# Weâll use an annotated PBMC reference dataset from Hao et al. (2020), available for download here: https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat
# Note that the SeuratDisk package is required to load the reference dataset. Installation instructions for SeuratDisk can be found here.
##################################################
library(SeuratDisk)

# load PBMC reference
reference <- LoadH5Seurat("/nfs/rprdata/julong/sc-atac/analyses.2021-02-05/pbmc_multimodal.h5seurat")

reference

colnames(reference@meta.data)

table(reference@meta.data$celltype.l1)
table(reference@meta.data$celltype.l2)
table(reference@meta.data$celltype.l3)

DefaultAssay(combined) <- "SCT"
combined2 <- combined
combined2

# transfer cell type labels from reference to query
transfer_anchors <- FindTransferAnchors(
          reference = reference,
          query = combined2,
          normalization.method = "SCT",
          reference.reduction = "spca",
          recompute.residuals = FALSE,
          dims = 1:50
        )


## sc3 <- MapQuery(anchorset=anchors, query=sc2, reference=ref,
##    refdata = list(celltype.l1="celltype.l1",
##    celltype.l2="celltype.l2",
##    predicted_ADT="ADT"),
##    reference.reduction="spca",
##    reduction.model="wnn.umap")
## x <- sc3@meta.data


predictions <- TransferData(anchorset = transfer_anchors,
                            refdata = reference$celltype.l2,
                            weight.reduction = combined2[['pca']],
                            dims = 1:50
                            )
#colnames(predictions)

predictions2 <- TransferData(anchorset = transfer_anchors,
                            refdata = reference$celltype.l1,
                            weight.reduction = combined2[['pca']],
                            dims = 1:50
                            )



combined2 <- AddMetaData(
          object = combined2,
          metadata = predictions %>% dplyr::select(predicted.id)
        )
#colnames(combined2@meta.data)

combined2$predicted.id.12 <- combined2$predicted.id
combined2$predicted.id <- NULL

combined2 <- AddMetaData(
          object = combined2,
          metadata = predictions2 %>% dplyr::select(predicted.id)
        )

combined2$predicted.id.11 <- combined2$predicted.id
combined2$predicted.id <- NULL

table(combined2$predicted.id.12)
table(combined2$predicted.id.11)

# set the cell identities to the cell type predictions
# by default the cells will be grouped by the following identifier
Idents(combined2) <- "predicted.id"

# set a reasonable order for cell types to be displayed when plotting
levels(combined2) <- c("CD4 Naive", "CD4 TCM", "CD4 CTL", "CD4 TEM", "CD4 Proliferating",
                                                               "CD8 Naive", "dnT",
                                                              "CD8 TEM", "CD8 TCM", "CD8 Proliferating", "MAIT", "NK", "NK_CD56bright",
                                                              "NK Proliferating", "gdT",
                                                              "Treg", "B naive", "B intermediate", "B memory", "Plasmablast",
                                                              "CD14 Mono", "CD16 Mono",
                                                              "cDC1", "cDC2", "pDC", "HSPC", "Eryth", "ASDC", "ILC", "Platelet")





############################################
# Joint UMAP visualization
############################################
# -----------------Clustering and visualization for joint analysis------------------------------------
# build a joint neighbor graph using both assays
combined2 <- FindMultiModalNeighbors(
      object = combined2,
      reduction.list = list("pca", "lsi"),
      dims.list = list(1:50, 2:50),
      modality.weight.name = "RNA.weight",
      verbose = TRUE
    )

head(combined2@meta.data)

# -----------------build a joint UMAP visualization------------------------------------
combined2 <- RunUMAP(
      object = combined2,
      nn.name = "weighted.nn",
      assay = "RNA",
      verbose = TRUE
    )


# png("./UMAP_joint_noclustering.png")
# DimPlot(combined2, reduction = "umap", label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
# dev.off()

for (x in c(0.2)) {
      combined2 <- FindClusters(combined2, graph.name = "wsnn", algorithm = 3, resolution = x, verbose = FALSE)
        png(paste("./3_Clustering.outs/UMAP_joint_clustering_mitofilt_", x, ".png", sep=""))
        p <- DimPlot(combined2, reduction = "umap", group.by = paste("wsnn_res.", x, sep=""), label = TRUE, repel = TRUE) + NoLegend()
        print(p)
        dev.off()
      }

colnames(combined2@meta.data)

# png("./UMAP_joint_clustering.png")
# DimPlot(pbmc.combined, reduction = "umap", group.by = "wsnn_res.0.15", label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
# dev.off()

## png("./3_Clustering.outs/UMAP_joint_labels_0.1.png")
## DimPlot(combined2, reduction = "umap", group.by = "predicted.id.11", label = TRUE, repel = TRUE) + NoLegend()
## dev.off()

## png("./3_Clustering.outs/UMAP_joint_labels_12_0.1.png")
## DimPlot(combined2, reduction = "umap", group.by = "predicted.id.12", label = TRUE, repel = TRUE) + NoLegend()
## dev.off()

png("./3_Clustering.outs/UMAP_joint_labels_mitofilt.png")
DimPlot(combined2, reduction = "umap", group.by = "predicted.id.11", label = TRUE, repel = TRUE) + NoLegend()
dev.off()

png("./3_Clustering.outs/UMAP_joint_labels_12_mitofilt.png")
DimPlot(combined2, reduction = "umap", group.by = "predicted.id.12", label = TRUE, repel = TRUE) + NoLegend()
dev.off()

colnames(combined2@meta.data)

#write_rds(combined2, file="./3_Clustering.outs/1_seurat.cluster.combined.rds")

write_rds(combined2, file="./3_Clustering.outs/1_seurat.cluster.combined.mitofilt.rds")

combined2 <- read_rds("./3_Clustering.outs/1_seurat.cluster.combined.mitofilt.rds")


#------------------------------------------
#------------------------------------------
###
### heatmap
#fn <- "./3_Clustering.outs/2_seurat.RNA.rds"
fn <- "./3_Clustering.outs/1_seurat.cluster.combined.mitofilt.rds"
combined2 <- read_rds(fn)
meta <- combined2@meta.data

colnames(meta)

table(meta$predicted.id.12)

#------------------------------------------
df1 <- meta%>%
   group_by(meta$predicted.id.l2, wsnn_res.0.1)%>%
   summarize(Freq=n(),.groups="drop")%>%
   group_by(wsnn_res.0.1)%>%
   mutate(Perc=Freq/sum(Freq))


p1 <- ggplot(df1, aes(x=wsnn_res.0.1, y="predicted.celltype.l2", fill=Perc))+
   geom_tile()+
   scale_fill_gradient("Fraction of cells",
      low="#ffffc8", high="#7d0025", na.value=NA)+
   xlab("Cluster")+ylab("predicted.celltype.l1")+
   theme_bw()+theme(axis.text.x=element_text(hjust=0.5, size=10))

###
figfn <- "./3_Clustering.outs/Figure3.1_heatmap_mitofilt.png"
png(figfn, width=800, height=600, res=120)
print(p1)
dev.off()

#------------------------------------------
###
df2 <- meta%>%
   group_by(meta$predicted.celltype.l2, wsnn_res.0.1)%>%
   summarize(Freq=n(),.groups="drop")%>%
   group_by(wsnn_res.0.1)%>%
   mutate(Perc=Freq/sum(Freq))

p2 <- ggplot(df2, aes(x=wsnn_res.0.1, y="predicted.celltype.l2", fill=Perc))+
   geom_tile()+
   scale_fill_gradient("Fraction of cells",
       low="#ffffc8", high="#7d0025", na.value=NA)+
   xlab("Cluster")+ylab("predicted.celltype.l2")+
   theme_bw()+theme(axis.text.x=element_text(hjust=0.5))

###
figfn <- "./3_Clustering.outs/Figure3.2_heatmap_mitofilt.png"
png(figfn, width=600, height=600, res=120)
print(p2)
dev.off()





############################################
# UMAP_by_celltypes/treatments
############################################
# # find all markers of cluster 2
# cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
# head(cluster2.markers, n = 5)


# # find all markers distinguishing cluster 5 from clusters 0 and 3
# cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
# head(cluster5.markers, n = 5)


# find markers for every cluster compared to all remaining cells, report only the positive
# ones

pbmc.combined <- combined2

pbmc.markers <- FindAllMarkers(pbmc.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
      group_by(cluster) %>%
      slice_max(n = 2, order_by = avg_log2FC)


write.table(pbmc.markers, file="pbmc_markers.txt", sep="\t")

# vln_plot
############################################
png("./3_Clustering.outs/Bcell_markers_vln.png", width=960, height=480)
VlnPlot(pbmc.combined, features = c("MS4A1", "CD79A")) + NoLegend()
dev.off()

# png("./3_Clustering.outs/Bcell_marker_MS4A1_vln.png")
# VlnPlot(pbmc.combined, features = "MS4A1") + NoLegend()
# dev.off()
#
# png("./3_Clustering.outs/Bcell_marker_CD79A_vln.png")
# VlnPlot(pbmc.combined, features = "CD79A") + NoLegend()
# dev.off()

#------------------------------------------------------------
png("./3_Clustering.outs/Monocyte_markers_vln.png", width=960, height=480)
VlnPlot(pbmc.combined, features = c("CD14", "MS4A7")) + NoLegend()
dev.off()

# png("./3_Clustering.outs/Monocyte_marker_CD14_vln.png")
# VlnPlot(pbmc.combined, features = "CD14") + NoLegend()
# dev.off()
#
# png("./3_Clustering.outs/Monocyte_marker_MS4A7_vln.png")
# VlnPlot(pbmc.combined, features = "MS4A7") + NoLegend()
# dev.off()

#------------------------------------------------------------
png("./3_Clustering.outs/NKcell_markers_vln.png", width=960, height=480)
VlnPlot(pbmc.combined, features = c("GNLY", "NKG7")) + NoLegend()
dev.off()

# png("./3_Clustering.outs/NKcell_marker_GNLY_vln.png")
# VlnPlot(pbmc.combined, features = "GNLY") + NoLegend()
# dev.off()
#
# png("./3_Clustering.outs/NKcell_marker_NKG7_vln.png")
# VlnPlot(pbmc.combined, features = "NKG7") + NoLegend()
# dev.off()

#------------------------------------------------------------
png("./3_Clustering.outs/Tcell_markers_vln.png", width=960, height=480)
VlnPlot(pbmc.combined, features = c("CD3D", "CD8A")) + NoLegend()
dev.off()

# png("./3_Clustering.outs/Tcell_marker_CD3D_vln.png")
# VlnPlot(pbmc.combined, features = "CD3D") + NoLegend()
# dev.off()
#
# png("./3_Clustering.outs/Tcell_marker_CD8A_vln.png")
# VlnPlot(pbmc.combined, features = "CD8A") + NoLegend()
# dev.off()



# feature_plot_markers
############################################
# png("./3_Clustering.outs/feature_plot.png")
# FeaturePlot(pbmc.combined, reduction = "umap", features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
#                                "CD8A"))
# dev.off()

png("./3_Clustering.outs/feature_plot_mitofilt.png", width=20, height=10, units="in", res=1200)
FeaturePlot(pbmc.combined, reduction = "umap", features = c("MS4A1", "CD14", "GNLY", "CD3D", "CD79A", "MS4A7", "NKG7", "CD8A"), ncol = 4)
dev.off()

# png("./3_Clustering.outs/feature_plot.png")
# FeaturePlot(pbmc.combined, reduction = "umap", features = c("MS4A1", "CD14", "GNLY", "CD3D", "CD79A", "MS4A7", "NKG7", "CD8A"))
# dev.off()

# png("./3_Clustering.outs/feature_plot_Bcell.png", width=10, height=5, units="in", res=300) #can adjust res only if width and height specifying in "in"
# FeaturePlot(pbmc.combined, reduction = "umap", features = c("MS4A1", "CD79A")) + NoLegend()
# dev.off()

png("./3_Clustering.outs/feature_plot_Bcell.png", width=960, height=480) #specifying widt and height in pixels
FeaturePlot(pbmc.combined, reduction = "umap", features = c("MS4A1", "CD79A")) + NoLegend()
dev.off()

png("./3_Clustering.outs/feature_plot_Monocyte.png", width=960, height=480)
FeaturePlot(pbmc.combined, reduction = "umap", features = c("CD14", "MS4A7")) + NoLegend()
dev.off()

png("./3_Clustering.outs/feature_plot_NKcell_.png", width=960, height=480)
FeaturePlot(pbmc.combined, reduction = "umap", features = c("GNLY", "NKG7")) + NoLegend()
dev.off()

png("./3_Clustering.outs/feature_plot_Tcell.png", width=960, height=480)
FeaturePlot(pbmc.combined, reduction = "umap", features = c("CD3D", "CD8A")) + NoLegend()
dev.off()



# feature_plot_treatments
############################################
colnames(pbmc.combined@meta.data)
head(pbmc.combined@meta.data$treats)

#treats.list <- unique(pbmc.combined@meta.data$treats)
#treats.list

#cell <- c("caffeine"="1", "nicotine"="2", "zinc"="3", "vitA"="4",
#          "vitE"="5", "vitD"="6", "water"="7", "etOH"="8")
cell <- c("caffeine"=1, "nicotine"=2, "zinc"=3, "vitA"=4,
          "vitE"=5, "vitD"=6, "water"=7, "etOH"=8)
pbmc.combined@meta.data$treats.num <- cell[as.character(meta$treats)]
colnames(pbmc.combined@meta.data)
head(pbmc.combined@meta.data$treats.num)
str(pbmc.combined@meta.data$treats.num)

#pbmc.combined <- transform(pbmc.combined@meta.data, treats.num = as.numeric(treats.num))


png("./3_Clustering.outs/dataset_mitofilt.png", width=24, height=3, units="in", res=600)
FeaturePlot(pbmc.combined, reduction = "umap", features = "treats.num", split.by = "treats", ncol=4, by.col=TRUE) + NoLegend()
dev.off()


png("./3_Clustering.outs/dataset_mitofilt_groupby.png")
DimPlot(pbmc.combined, reduction = "umap", group.by = 'treats', pt.size = 0.1)
dev.off()


############################################
pbmc.markers %>%
      group_by(cluster) %>%
      top_n(n = 10, wt = avg_log2FC) -> top10
png("./top10.png")
DoHeatmap(pbmc.combined, features = top10$gene) + NoLegend()
dev.off()
