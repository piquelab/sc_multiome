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
library(SummarizedExperiment)
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

## gene <- rownames(sc)

## x0 <- c("CD3D", "CD3E", "CD3G", "GNLY", "MS4A1", "CD14")
## x0 <- c("CD3D", "GNLY", "MS4A1", "CD14")






#############################################################
### Heatmap and genomic track for cell types marker genes ###
#############################################################

### single cell RNA+ATAC
 
rm(list=ls())

outdir <- "./1_processing.outs/"

fn <- "./sc_multiome_data/1_processing/3_Clustering.outs/1_seurat.cluster.combined.mitofilt.rds"    
sc <- read_rds(fn)
DefaultAssay(sc) <- "RNA"
Idents(sc) <- sc$wsnn_res.0.1


###
### Find differentially expressed marker genes

## markers <- FindAllMarkers(sc, only.pos=T, min.pct=0.1, logfc.threshold=0.1)

## opfn <- paste(outdir, "1.0_wsnn_res0.1_allmarkers.rds", sep="")
## markers <- markers%>%group_by(cluster)%>%arrange(desc(avg_log2FC), .by_group=T)%>%as.data.frame()
## write_rds(markers, opfn)

## ###
## fn <- "./1_processing.outs/1.0_wsnn_res0.1_allmarkers.rds"
## x <- read_rds(fn)
## x2 <- x%>%group_by(cluster)%>%slice_max(order_by=abs(avg_log2FC), n=50)%>%ungroup()
## x2 <- x2%>%mutate(MCls=MCls_name[as.character(cluster)])
## opfn <- paste(outdir, "1.2_top50.xlsx", sep="")
## write.xlsx(x2, file=opfn)




x <- sc@assays$RNA@counts
rn <- rownames(x)

###
### cell-type value
MCls_name <- c("0"="0-CD4Naive", "1"="1-TCM", "2"="2-NKcell", "3"="3-TEM",
    "4"="4-Bcell", "5"="5-CD8Naive", "6"="6-Monocyte", "7"="7-dnT",
    "8"="8-MAIT", "9"="9-Platelets", "10"="10-DC")

col_MCls <- c("0-CD4Naive"="#ffaa00", "1-TCM"="pink",  "2-NKcell"="#aa4b56", 
   "3-TEM"="blue", "4-Bcell"="#4daf4a", "5-CD8Naive"="green", "6-Monocyte"="#984ea3",
   "7-dnT"="black", "8-MAIT"="#ffffb2", "9-Platelets"="#8dd3c7", "10-DC"="#d4b9da")

MCls_val <- c("0-CD4Naive"=1, "5-CD8Naive"=2, "1-TCM"=3, "3-TEM"=4, "7-dnT"=5, "8-MAIT"=6,
              "2-NKcell"=7, "4-Bcell"=8, "6-Monocyte"=9, "10-DC"=10, "9-Platelets"=11)


###
###
## CD4 naive, "IL7R", "CCR7", "CD62L", "CD45RA"  
## "S100A4", "CCR10", "TNFRSF18", ""
## CD8 Naive, "CD8A", "CD3D", "ID3"
## NK cell, "GNLY", "NKG7"
## B cell, "MS4A1", "CD79A",
## Mono, "CD14", "MS4A7", "CD16", "LYZ", "FCGR3A", "S100A8", "S100A9"
## DC, "FCER1A", "CST3", 
## Platelet, "PPBP", 

###
### gene value
geneList <- c("TSHZ2", "FHIT", "PDE3B", "THEMIS", "CAMK4", "BACH2", "SLC4A10",
     "IL7R", "CCR7", "CD3D",
     "LINGO2", "NKG7", "GNLY", "EBF1", "MS4A1", "CD79A","VCAN", "CD14", "LYZ", "MS4A7", "S100A8", "S100A9",
     "CCL22", "FCER1A", "CST3", "IGHA1")
nlen <- length(geneList)

gene_val <- 1:nlen
names(gene_val) <- geneList



###
### plot data 
data <- FetchData(sc, vars=geneList, slot="data")

meta <- sc@meta.data
meta <- meta%>%mutate(MCls=MCls_name[as.character(wsnn_res.0.1)])




### data for plots 
MCls <- sort(unique(meta$MCls))
DF_plot <- map_dfr(MCls, function(ii){
   ##
   cellSel <- meta%>%filter(MCls==ii)%>%pull(NEW_BARCODE)
   xi <- data[cellSel,]
   DF <- data.frame(MCls=ii, gene=colnames(data),
      yave=apply(xi, 2, mean), percent=apply(xi>0, 2, mean))
   ##DF <- DF%>%mutate(yave2=ifelse(percent<1e-6, yave, yave/percent))
   DF 
})
rownames(DF_plot) <- NULL

 
DF_plot <- DF_plot%>%group_by(gene)%>%
    mutate(yave2=yave/max(yave))%>%as.data.frame()
###
DF_plot2 <- DF_plot%>%
   mutate(MCls_value=as.numeric(MCls_val[MCls]),
          gene_value=as.numeric(gene_val[gene]),
          MCl2=fct_reorder(MCls, MCls_value),
          gene2=fct_reorder(gene, gene_value),
          yave2=ifelse(yave<1e-6, NA, yave2),
          gr2=case_when(percent<=0.15~"1",
                        percent>0.15&percent<=0.5~"2",
                        percent>0.5&percent<=0.75~"3",
                        percent>0.75~"4"),
          gr2=ifelse(percent<0.05, NA, gr2))

          ## size2=ifelse(percent2<0.6, NA, percent2),
          ## yave2=ifelse(yave2<0.6, NA, yave2),
          ## gr2=ifelse(is.na(yave2), NA, "gr2"))

###
### setting colors 
mycol <- viridis(20)
 
p <- ggplot(DF_plot2, aes(x=MCl2, y=gene2))+
     geom_point(aes(color=yave2, size=factor(gr2)),  shape=19)+
     scale_color_gradientn(name="Relative\nexpression",
          colours=mycol, na.value=NA,
          breaks=c(0, 1),  limits=c(0,1), labels=c("Min.", "Max."),
          guide=guide_colourbar(barwidth=grid::unit(0.8,"lines"),
              barheight=grid::unit(5,"lines"), order=1)
          )+
    scale_size_manual(name="Percents(%)", values=c("1"=2.5, "2"=3.5, "3"=5.5, "4"=6.5),
        labels=c("1"="<15", "2"="<50", "3"="<75", "4"="<100"), na.value=NA,
        guide=guide_legend(order=2))+
##    scale_size_binned(range=c(0.1, 7), limits=c(0.6,1), guide="none")+
    ##scale_color_manual(values=c("gr2"="black"), na.value=NA,  guide="none")+    
     ##scale_size_manual(values=c("a1"=0.1, "a2"=4), guide="none")+ 
    theme(axis.title=element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, size=10),
          axis.text.y=element_text(size=10),
          legend.background=element_blank(),
          legend.key=element_blank(), 
          panel.background=element_blank(),
          panel.border=element_rect(color="black", fill=NA, size=1.5),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank())
          #legend.text=element_text(size=8))
 
##
figfn <- paste(outdir, "Figure1.3_dotplot.png", sep="")
ggsave(figfn, p, width=700, height=850, units="px", dpi=120)



#####################################
### Create a gene activity matrix ###
#####################################

## fn <- "./sc_multiome_data/1_processing/5.1_reCallPeak.outs/3_scATAC.annot.mitofilt.ndup.rds"
## sc <- read_rds(fn)
meta2 <- meta%>%dplyr::select(NEW_BARCODE, gex_barcode, EXP, MCls)%>%
    mutate(NEW_BARCODE_ArchR=paste0(EXP, "#", gex_barcode))

meta_ArchR <- read_rds("./1.2_ArchR_process/1_meta.data.rds")%>%
    dplyr::select(NEW_BARCODE_ArchR, cluster_signac, cluster_ArchR)        

meta_ArchR <- meta_ArchR%>%inner_join(meta2, by="NEW_BARCODE_ArchR")
id1 <- meta_ArchR$NEW_BARCODE_ArchR


###
### gene score
mat <- read_rds("./1.2_ArchR_process/1.2_genescore.rds")

id2 <- colnames(mat)
rn <- rowData(mat)$name

data <- assays(mat)$GeneScoreMatrix
rownames(data) <- rn

geneList2 <- geneList[geneList%in%rn] 
data2 <- data[geneList2,]
data2 <- t(data2)


### data for plots 
MCls <- sort(unique(meta$MCls))
DF_plot <- map_dfr(MCls, function(ii){
   ##
   cellSel <- meta_ArchR%>%filter(MCls==ii)%>%pull(NEW_BARCODE_ArchR)
   xi <- data2[cellSel,]
   DF <- data.frame(MCls=ii, gene=colnames(data2),
      yave=apply(xi, 2, mean), percent=apply(xi>0, 2, mean))
   ##DF <- DF%>%mutate(yave2=ifelse(percent<1e-6, yave, yave/percent))
   DF 
})
rownames(DF_plot) <- NULL

   
DF_plot <- DF_plot%>%group_by(gene)%>%
    mutate(yave2=yave/max(yave))%>%as.data.frame()
###
DF_plot2 <- DF_plot%>%
   mutate(MCls_value=as.numeric(MCls_val[MCls]),
          gene_value=as.numeric(gene_val[gene]),
          MCl2=fct_reorder(MCls, MCls_value),
          gene2=fct_reorder(gene, gene_value),
          yave2=ifelse(yave<0.1, NA, yave2),
          gr2=case_when(percent<=0.25~"1",
                        percent>0.25&percent<=0.5~"2",
                        percent>0.5&percent<=0.75~"3",
                        percent>0.75~"4"),
          gr2=ifelse(percent<0.1, NA, gr2))

          ## size2=ifelse(percent2<0.6, NA, percent2),
          ## yave2=ifelse(yave2<0.6, NA, yave2),
          ## gr2=ifelse(is.na(yave2), NA, "gr2"))

###
### setting colors 
mycol <- viridis(20)
 
p2 <- ggplot(DF_plot2, aes(x=MCl2, y=gene2))+
     geom_point(aes(color=yave2, size=factor(gr2)),  shape=19)+
     scale_color_gradientn(name="Relative\nexpression",
          colours=mycol, na.value=NA,
          breaks=c(0, 1),  limits=c(0,1), labels=c("Min.", "Max."),
          guide=guide_colourbar(barwidth=grid::unit(0.8,"lines"),
              barheight=grid::unit(5,"lines"), order=1)
          )+
    scale_size_manual(name="Percents(%)", values=c("1"=2.5, "2"=3.5, "3"=5, "4"=6.5),
        labels=c("1"="<25", "2"="<50", "3"="<75", "4"="<100"), na.value=NA,
        guide=guide_legend(order=2))+
##    scale_size_binned(range=c(0.1, 7), limits=c(0.6,1), guide="none")+
    ##scale_color_manual(values=c("gr2"="black"), na.value=NA,  guide="none")+    
     ##scale_size_manual(values=c("a1"=0.1, "a2"=4), guide="none")+ 
    theme(axis.title=element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, size=10),
          axis.text.y=element_text(size=10),
          legend.background=element_blank(),
          legend.key=element_blank(), 
          panel.background=element_blank(),
          panel.border=element_rect(color="black", fill=NA, size=1.5),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank())
          #legend.text=element_text(size=8))
 
##
figfn <- paste(outdir, "Figure1.3_dotplot_genescore.png", sep="")
ggsave(figfn, p2, width=700, height=850, units="px", dpi=120)



               

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

figfn <- paste(outdir, "FigS1_2_umap_treat.png", sep="")
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

figfn <- paste(outdir, "FigS1_3_comb_ref2.png", sep="")
png(figfn, width=1250, height=460, res=120)
plot_grid(p1, p2, p3, ncol=3, labels="AUTO", label_fontface="plain", label_x=0.1,
    align="h", axis="tb", rel_widths=c(1,1,1.3))
dev.off()





#############################################
### summary of RNA and ATAC data quanlity ###
#############################################


rm(list=ls())

outdir <- "./1_processing.outs/"

fn <- "./sc_multiome_data/1_processing/3_Clustering.outs/1_seurat.cluster.combined.mitofilt.rds"    
sc <- read_rds(fn)
DefaultAssay(sc) <- "RNA"


### RNA counts data 
x <- sc@assays$RNA@counts
rnz <- rowSums(x)
x <- x[rnz>0,]

meta <- sc@meta.data

meta2 <- sc@meta.data%>%select(NEW_BARCODE, sampleID=SNG.BEST.GUESS, treats)%>%
   mutate(nCount_RNA=colSums(x), nFeature_RNA=colSums(x>0))

###
### plot data
dd2 <- meta2%>%group_by(treats, sampleID)%>%
    summarise(ny=n(), nreads=mean(nCount_RNA), ngenes=mean(nFeature_RNA), .groups="drop")


### output data for summary statistics 
summ2 <- dd2%>%group_by(treats)%>%
    summarise(nindi=n(), ncells=median(ny), reads=median(nreads), genes=median(ngenes), .groups="drop")

opfn <- paste(outdir, "TableS1_1_summ_data.xlsx", sep="")
write.xlsx(summ2, file=opfn, overwrite=T)





###
### supp plots 

col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
       "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3", "etOH"="grey50", "water"="grey50")
treat_val <- c("caffeine"=1, "nicotine"=2, "vitA"=3, "vitD"=4, "vitE"=5, "zinc"=6, "etOH"=7, "water"=7)

dd2 <- dd2%>%mutate(treat_value=as.numeric(treat_val[as.character(treats)]),
                   treat2=fct_reorder(treats, treat_value))
 

p1 <- ggplot(dd2, aes(x=treat2, y=ny, fill=treat2))+
   geom_violin()+xlab("")+ylab("")+ylim(0,1500)+
   geom_boxplot(width=0.2, color="grey", outlier.shape=NA)+ 
   ggtitle("#Cells per individual")+
   scale_fill_manual(values=col2)+
   ##scale_x_discrete(labels=lab1)+ylim(0,3000)+
   theme_bw()+
   theme(legend.position="none",
   plot.title=element_text(hjust=0.5, size=14),
   axis.text.x=element_text(angle=45, hjust=1, size=12),
   axis.text.y=element_text(size=12))

###
p2 <- ggplot(dd2,aes(x=treat2, y=nreads, fill=treat2))+
   geom_violin()+xlab("")+ylab("")+ylim(1000, 6000)+
   geom_boxplot(width=0.2, color="grey", outlier.shape=NA)+     
   ggtitle("#Reads per cell")+
   scale_fill_manual(values=col2)+
   ##scale_x_discrete(labels=lab1)+
   theme_bw()+
   theme(legend.position="none",
         plot.title=element_text(hjust=0.5, size=14),
         axis.text.x=element_text(angle=45, hjust=1, size=12),
         ## axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5, size=12),
         axis.text.y=element_text(size=12))

###
p3 <- ggplot(dd2,aes(x=treat2, y=ngenes, fill=treat2))+
    geom_violin()+xlab("")+ylab("")+ylim(500,2500)+
    geom_boxplot(width=0.2, color="grey", outlier.shape=NA)+     
    ggtitle("#Genes per cell")+
    scale_fill_manual(values=col2)+
    ##scale_x_discrete(labels=lab1)+
    theme_bw()+
    theme(legend.position="none",
          plot.title=element_text(hjust=0.5, size=14),
          axis.text.x=element_text(angle=45, hjust=1, size=12), 
          ## axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5, size=12),
          axis.text.y=element_text(size=12))


## png("./2_kb2_output/Figure2.5_violin.png", width=800, height=500, res=120) 
figfn <- paste(outdir, "FigS1_1_violin.pdf", sep="")
pdf(figfn, width=10, height=4.5)
print(plot_grid(p1, p2, p3, labels="AUTO", label_fontface="plain", label_x=0.1,  ncol=3))
dev.off()



################################
#### summarize the peak data ###
################################

fn <- "./sc_multiome_data/1_processing/5.1_reCallPeak.outs/3_scATAC.annot.mitofilt.ndup.rds"
sc <- read_rds(fn)

## DefaultAssay(sc) <- "RNA"

### ATAC counts data 
x <- sc@assays$ATAC@counts
rnz <- rowSums(x)
x <- x[rnz>0,]

meta <- sc@meta.data

meta2 <- sc@meta.data%>%select(NEW_BARCODE, sampleID=SNG.BEST.GUESS, treats)%>%
   mutate(nCount_ATAC=colSums(x), nFeature_ATAC=colSums(x>0))

###
### plot data
dd2 <- meta2%>%group_by(treats, sampleID)%>%
    summarise(ny=n(), nreads=mean(nCount_ATAC), ngenes=mean(nFeature_ATAC), .groups="drop")


### output data for summary statistics 
summ2 <- dd2%>%group_by(treats)%>%
    summarise(nindi=n(), ncells=median(ny), reads=median(nreads), genes=median(ngenes), .groups="drop")

opfn <- paste(outdir, "TableS1_2_ATAC_summ_data.xlsx", sep="")
write.xlsx(summ2, file=opfn, overwrite=T)


###
### supp plots for ATAC data summary quality 
col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
       "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3", "etOH"="grey50", "water"="grey50")
treat_val <- c("caffeine"=1, "nicotine"=2, "vitA"=3, "vitD"=4, "vitE"=5, "zinc"=6, "etOH"=7, "water"=7)

dd2 <- dd2%>%mutate(treat_value=as.numeric(treat_val[as.character(treats)]),
                   treat2=fct_reorder(treats, treat_value))
 
## p1 <- ggplot(dd2, aes(x=treat2, y=ny, fill=treat2))+
##    geom_violin()+xlab("")+ylab("")+ylim(0,1500)+
##    geom_boxplot(width=0.2, color="grey", outlier.shape=NA)+ 
##    ggtitle("#Cells per individual")+
##    scale_fill_manual(values=col2)+
##    ##scale_x_discrete(labels=lab1)+ylim(0,3000)+
##    theme_bw()+
##    theme(legend.position="none",
##    plot.title=element_text(hjust=0.5, size=14),
##    axis.text.x=element_text(angle=45, hjust=1, size=12),
##    axis.text.y=element_text(size=12))

###
p2 <- ggplot(dd2,aes(x=treat2, y=nreads, fill=treat2))+
   geom_violin()+xlab("")+ylab("")+ylim(10000, 20000)+
   geom_boxplot(width=0.15, color="grey", outlier.shape=NA)+     
   ggtitle("#Reads per cell")+
   scale_fill_manual(values=col2)+
   ##scale_x_discrete(labels=lab1)+
   theme_bw()+
   theme(legend.position="none",
         plot.title=element_text(hjust=0.5, size=14),
         axis.text.x=element_text(angle=45, hjust=1, size=12),
         ## axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5, size=12),
         axis.text.y=element_text(size=12))

###
p3 <- ggplot(dd2,aes(x=treat2, y=ngenes, fill=treat2))+
    geom_violin()+xlab("")+ylab("")+ylim(7000,15000)+
    geom_boxplot(width=0.15, color="grey", outlier.shape=NA)+     
    ggtitle("#Peaks per cell")+
    scale_fill_manual(values=col2)+
    ##scale_x_discrete(labels=lab1)+
    theme_bw()+
    theme(legend.position="none",
          plot.title=element_text(hjust=0.5, size=14),
          axis.text.x=element_text(angle=45, hjust=1, size=12), 
          ## axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5, size=12),
          axis.text.y=element_text(size=12))


## png("./2_kb2_output/Figure2.5_violin.png", width=800, height=500, res=120) 
figfn <- paste(outdir, "FigS1_4_ATAC_violin.pdf", sep="")
pdf(figfn, width=8, height=4.8)
print(plot_grid(p2, p3, labels="AUTO", label_fontface="plain", label_x=0.1,  ncol=2))
dev.off()







### summary number of cells for each cell-types
df2 <- meta%>%select(NEW_BARCODE, clusters=wsnn_res.0.1)%>%
    group_by(clusters)%>%summarise(ny=n(),.groups="drop")
df2 <- df2%>%mutate(percent=round(ny/sum(ny), 3))%>%arrange(desc(ny))

opfn <- paste(outdir, "TableS1_2_cells_MCls.xlsx", sep="")
write.xlsx(df2, file=opfn, overwrite=T)





###
###
a0 <- 2
y1 <- a0+rnorm(12)*1.5
y2 <- 2*a0+rnorm(12)*1.6

df2 <- data.frame(x=rep(c("gr1", "gr2"), each=12), y=c(y1, y2))
df2 <- df2%>%mutate(gr=x)

p0 <- ggplot(df2, aes(x, y, fill=gr))+
   geom_violin(aes(color=gr))+ 
   geom_boxplot(outlier.shape=NA, width=0.3, color="grey30", lwd=0.25)+
   ##geom_jitter(width=0.2, size=0.5)+
   scale_fill_manual(values=c("gr1"="grey70", "gr2"="maroon3"))+
   scale_color_manual(values=c("gr1"="grey70", "gr2"="maroon3"))+    
   theme_bw()+
   theme(legend.position="none",
         axis.title=element_blank(),
         axis.text=element_blank(),
         axis.ticks=element_blank())
 
 figfn <- paste(outdir, "Figure1_demo_peak.png", sep="")
 ggsave(figfn, p0, width=320, height=320, units="px", dpi=120)




### gene expression 
a0 <- 5
y1 <- a0+rnorm(12)*1.5
y2 <- 2*a0+rnorm(12)*1.6

df2 <- data.frame(x=rep(c("gr1", "gr2"), each=12), y=c(y1, y2))
df2 <- df2%>%mutate(gr=x)

p1 <- ggplot(df2, aes(x, y, fill=gr))+
   geom_violin(aes(color=gr))+ 
   geom_boxplot(outlier.shape=NA, width=0.3, color="grey30", lwd=0.25)+
   ##geom_jitter(width=0.2, size=0.5)+
   scale_fill_manual(values=c("gr1"="grey70", "gr2"="maroon3"))+
   scale_color_manual(values=c("gr1"="grey70", "gr2"="maroon3"))+    
   theme_bw()+
   theme(legend.position="none",
         axis.title=element_blank(),
         axis.text=element_blank(),
         axis.ticks=element_blank())
 
 figfn <- paste(outdir, "Figure1_demo_gene.png", sep="")
 ggsave(figfn, p1, width=320, height=320, units="px", dpi=120)


