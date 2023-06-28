##
library(Matrix)
library(tidyverse)
library(data.table)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(Signac)
library(SeuratWrappers)
library(SeuratObject)
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

outdir <- "./1_processing.outs/"
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)



### scmultiome data
fn <- "./sc_multiome_data/1_processing/3_Clustering.outs/1_seurat.cluster.combined.mitofilt.rds"    
sc <- read_rds(fn)
##
meta <- sc@meta.data
umap <- Embeddings(sc, reduction="umap")

MCls_name <- c("0"="0-CD4Naive", "1"="1-TCM", "2"="2-NKcell", "3"="3-TEM",
    "4"="4-Bcell", "5"="5-CD8Naive", "6"="6-Monocyte", "7"="7-dnT",
    "8"="8-MAIT", "9"="9-Platelets", "10"="10-DC")

col_MCls <- c("0-CD4Naive"="#ffaa00", "1-TCM"="pink",  "2-NKcell"="#aa4b56", 
   "3-TEM"="blue", "4-Bcell"="#4daf4a", "5-CD8Naive"="green", "6-Monocyte"="#984ea3",
   "7-dnT"="black", "8-MAIT"="#ffffb2", "9-Platelets"="#8dd3c7", "10-DC"="#d4b9da")

              
plotDF <- data.frame(umap[,1:2], cluster=as.character(meta$wsnn_res.0.1))%>%
     mutate(MCls=MCls_name[cluster], treats=meta$treats)



### main
p1 <- ggplot(plotDF, aes(x=UMAP_1, y=UMAP_2))+
        geom_point(aes(colour=factor(MCls)), size=0.1)+
        scale_colour_manual(values=col_MCls, guide="none")+
        theme_bw()+
        theme(legend.position="none",
              axis.title=element_text(size=12),
              axis.text=element_text(size=12))

figfn <- paste(outdir, "Figure1.1_umap_cl.png", sep="")
png(figfn, width=420, height=420, res=120)
p1
dev.off()



###
### marker genes expression

## Naive CD4+, IL7R, CCR7
## Memory CD4+, IL7R, S100A4
## CD8+, CD8A,
## NK cell, GNLY, NKG7
## B cell, MS4A1
## CD14+ Mono, CD14, LYZ
## FCGR3A+ Mono, FCGR3A, MS4A7
## DC, FCER1A, CST3
## Platelet, PPBP

gene <- rownames(sc)

x0 <- c("CD3D", "CD3E", "CD3G", "GNLY", "MS4A1", "CD14")
x0 <- c("CD3D", "GNLY", "MS4A1", "CD14")



               

##################
### supp-plots ###
##################


## UMAP facet by treatment
treat_val <- c("caffeine"=1, "nicotine"=2, "vitA"=3, "vitD"=4, "vitE"=5, "zinc"=6, "etOH"=7, "water"=7)
plotDF2 <- plotDF%>%
    mutate(treat_value=as.numeric(treat_val[treats]), treat2=fct_reorder(treats, treat_value))

p2 <- ggplot(plotDF2, aes(x=UMAP_1, y=UMAP_2))+
   geom_point(aes(colour=factor(MCls)), size=0.1)+
   facet_wrap(~factor(treat2), scales="fixed", ncol=4)+
   theme_bw()+
   scale_colour_manual(values=col_MCls, guide=guide_legend(override.aes=list(size=3)))+
   theme(legend.title=element_blank(),
         legend.key.size=unit("0.5", "cm"),
         legend.background=element_blank(),
         legend.key=element_blank(),
         legend.box.background=element_blank(),
         axis.title=element_text(size=12),
         axis.text=element_text(size=12),
         strip.text=element_text(size=14))

figfn <- paste(outdir, "FigS1_1_umap_treat.png", sep="")
png(figfn, width=920, height=520, res=120)
p2
dev.off()


##########################################
### S1_2_UMAP for cell type annotation ###
##########################################

MCls_name <- c("0"="0-CD4Naive", "1"="1-TCM", "2"="2-NKcell", "3"="3-TEM",
    "4"="4-Bcell", "5"="5-CD8Naive", "6"="6-Monocyte", "7"="7-dnT",
    "8"="8-MAIT", "9"="9-Platelets", "10"="10-DC")

col_MCls <- c("0-CD4Naive"="#ffaa00", "1-TCM"="pink",  "2-NKcell"="#aa4b56", 
   "3-TEM"="blue", "4-Bcell"="#4daf4a", "5-CD8Naive"="green", "6-Monocyte"="#984ea3",
   "7-dnT"="black", "8-MAIT"="#feb24c", "9-Platelets"="#8dd3c7", "10-DC"="#d4b9da")

comb <- c("CD4 Naive"="CD4 Naive", "CD4 TCM"="TCM", "CD4 TEM"="TEM",
    "CD8 Naive"="CD8 Naive", "CD8 TCM"="TCM", "CD8 TEM"="TEM",
    "CD4 CTL"="CD4 CTL",  "ILC"="ILC",                             
    "MAIT"="MAIT", "gdT"="gdT", "dnT"="dnT", "Treg"="Treg",    
    "B intermediate"="B cell",
    "B memory"="B cell",
    "B naive"="B cell",
    "CD14 Mono"="Mono",
    "CD16 Mono"="Mono",
    "pDC"="DC", "cDC1"="DC", "cDC2"="DC", "ASDC"="DC", "HSPC"="HSPC",
    "NK"="NK", "NK Proliferating"="NK", "NK_CD56bright"="NK",
    "Plasmablast"="Plasmablast", "Platelet"="Platelet", "Eryth"="Eryth")


### data for plots
x <- data.frame(umap[,1:2], cluster=as.character(meta$wsnn_res.0.1))%>%
    mutate(MCls=MCls_name[cluster], predicted.id=meta$predicted.id.12,
           predicted.id2=comb[as.character(predicted.id)], NEW_BARCODE=meta$NEW_BARCODE)
rownames(x) <- meta$NEW_BARCODE
sc2 <- AddMetaData(sc, metadata=x)

mycol1 <- c("#ffaa00", "pink",  "#d4b9da", "#aa4b56", "blue", "#4daf4a", "green", "#984ea3",
   "black", "#ffffb2", "#8dd3c7")

## mycol1 <- c("0"="#ffaa00", "1"="pink",  "10"="#d4b9da", "2"="#aa4b56", "3"="blue",
##    "4"="#4daf4a", "5"="green", "6"="#984ea3",
##    "7"="black", "8"="#feb24c", "9"="#8dd3c7")

p1 <- DimPlot(sc2, reduction="umap", cols=mycol1, group.by="cluster",
              label=T, repel=T, pt.size=T, label.size=3.5, raster=F)+
    ggtitle("Seurat clusters")+
    theme_bw()+
    theme(legend.position="none",
          plot.title=element_text(hjust=0.5, size=12),
          axis.title=element_text(size=10),
          axis.text=element_text(size=10))


##
mycol2 <- c("#4daf4a", "#045a8d", "#ffaa00", "green",
  "#d4b9da", "black", "#de2d26",
   "grey", "#f0027f", "#74a9cf",
    "#ffffb2", "#984ea3", "#aa4b56", "#01665e", "#8dd3c7", "pink", "blue", "#f768a1")

### annotation
## mycol2 <- c("B cell"="#4daf4a", "CD4 CTL"="#045a8d", "CD4 Naive"="#ffaa00", "CD8 Naive"="green",
##             "DC"="#d4b9da", "dnT"="black", "Eryth"="#54278f",
##    "gdT"="grey", "HSPC"="#f1eef6", "ILC"="#74a9cf",
##     "MAIT"="#feb24c", "Mono"="#984ea3", "NK"="#aa4b56", "Plasmablast"="#006d2c",
##    "Platelet"="#8dd3c7", "TCM"="pink", "TEM"="blue", "Treg"="#f768a1")
 
p2 <- DimPlot(sc2, reduction="umap", cols=mycol2, group.by="predicted.id2",
              label=T, repel=T, pt.size=T, label.size=2.5, raster=F)+
    ggtitle("Automatic annotation")+
    theme_bw()+
    theme(legend.position="none",
          plot.title=element_text(hjust=0.5, size=12),
          axis.title=element_text(size=10),
          axis.text=element_text(size=10))


###
### heatmap

CL_value <- c("0"=1, "1"=2, "3"=3, "5"=4, "7"=5, "8"=6,
              "2"=7, "4"=8, "6"=9, "10"=10, "9"=11)

MCl_value <- c("CD4 Naive"=1, "TCM"=2, "TEM"=3, "CD8 Naive"=4, "dnT"=8,
   "MAIT"=5, "Treg"=6, "gdT"=7,  "CD4 CTL"=9,
   "NK"=10, "B cell"=11, "Mono"=12, "DC"=13, "Platelet"=14, "Eryth"=15)

###
### plot data for heatmap
x2 <- x%>%dplyr::select(NEW_BARCODE, predicted.id2, cluster)

df0 <- x2%>%group_by(predicted.id2)%>%summarize(n_MCls=n())%>%filter(n_MCls>100)

df2 <- x2%>%group_by(predicted.id2, cluster)%>%
    summarize(Freq=n(), .groups="drop")%>%filter(predicted.id2%in%df0$predicted.id2)%>%
    group_by(cluster)%>%
    mutate(Perc=Freq/sum(Freq))%>%ungroup()

###
### order x-axis and y-axis
df2 <- df2%>%
    mutate(CL_val=as.numeric(CL_value[as.character(cluster)]),
           cluster2=fct_reorder(as.character(cluster), CL_val),
           predicted.id2_val=as.numeric(MCl_value[as.character(predicted.id2)]),
           predicted.id3=fct_reorder(predicted.id2, predicted.id2_val))


p3 <- ggplot(df2, aes(x=cluster2, y=predicted.id3, fill=Perc))+
    geom_tile()+
    scale_fill_gradient(low="#ffffcc", high="#e31a1c",
      na.value=NA, limits=c(0,1), n.breaks=5,
      guide=guide_legend(keywidth=grid::unit(0.5,"lines"), keyheight=grid::unit(1,"lines") ))+
    xlab("Cluster id")+ ##ylab("Predicted.celltype")+
    ggtitle("Fraction of cells")+
    theme_bw()+
    theme(legend.title=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_text(size=10),
          axis.text.y=element_text(size=10),
          plot.title=element_text(hjust=0.5, size=12))


###
### output

figfn <- paste(outdir, "FigS1_2_comb_ref2.png", sep="")
png(figfn, width=1250, height=460, res=120)
plot_grid(p1, p2, p3, ncol=3, labels="AUTO", label_fontface="plain", label_x=0.1,
    align="h", axis="tb", rel_widths=c(1,1,1.3))
dev.off()
