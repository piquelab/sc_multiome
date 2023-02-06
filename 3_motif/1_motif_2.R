###
library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(SeuratObject)
library(Signac)
library(SeuratWrappers)
library(SeuratData)
library(qqman)
library(JASPAR2020)
library(JASPAR2022)
library(TFBSTools)
## library(JASPAR2022, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/") 
## library(TFBSTools, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
library(BSgenome.Hsapiens.1000genomes.hs37d5, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(ChIPseeker, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
###
library(ggplot2)
library(cowplot) ##lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(grid)
library(gridExtra)
library(ggExtra)
library(RColorBrewer)
library(ggsci)
library(viridis)
library(ComplexHeatmap)
library(circlize)
theme_set(theme_grey())

outdir <- "./1_motif.outs/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


###########################
### obtain motif object ###
###########################
##
opfn <- "./1_motif.outs/1_scATAC.motif.rds"
atac <- read_rds(opfn)

##rm(list=ls())
atac
atac[["ATAC"]]
head(atac@meta.data)
table(atac@meta.data$MCls)


motif <- Motifs(atac)
motif







## ########################################
## ### MOTIF EXAMPLE                    ###
## ########################################
## #----MOTIF FILE SAVED IN THIS SCRIPT----
## motif <- Motifs(atac)
## motif
## #head(motif)

## motif.col <- colnames(motif)
## motif.col
## length(motif.col)

## #----count----
## motif.matrix <- GetMotifData(object = atac)
## motif.matrix
## #motif_data <- SetMotifData(object = atac, assay = 'ATAC', slot = 'data', new.data = motif.matrix)
## #motif_data

## motif.matrix.df <- as.data.frame(motif.matrix)
## head(motif.matrix.df)

## motif.matrix.peak <- rownames(motif.matrix) 
## motif.matrix.peak

## cvt <- str_split(rownames(motif.matrix), "-", simplify=T)
## table(cvt[,1])

## #----example----
## example <- c("MA0137.3", "MA0517.1", "MA0144.2",
##              "MA0105.4", "MA0778.1", "MA0099.3",
##              "MA1563.2", "MA1134.1", "MA1141.1")


## example %in% motif.col

## #motif.matrix <- as.matrix(motif)

## fig <- MotifPlot(object=atac,motifs=example, facet="wrap",ncol=3,nrow=3)

## figfn <- "./1_motif.outs/Figure1.2_motif.example.png"
## png(figfn, width=800, height=600,res=120)
## print(fig)
## dev.off()







############################################################################
############################################################################
### enrichment analysis for motif ###
############################################################################
############################################################################
##----enriched motif file----
enriched.motif <- read_rds("./1_motif.outs/2_motif.enrich_macs2_0.1_cn_control.rds")
res <- enriched.motif %>% as.data.frame()

head(res)

dim(res)
colnames(res)


#max(res$qvalue.fisher)
#min(res$qvalue.fisher)
max(res$fold.enrichment)
min(res$fold.enrichment)
#nrow(res)
#nrow(res%>%dplyr::filter(qvalue.fisher<=0.1))



res <- res %>% drop_na(pvalue) %>% mutate(MCls=gsub("_", "-", MCls)) %>% mutate(comb=paste(MCls, contrast, sep="_"))
res <- res %>% mutate(logFC=log2(fold.enrichment))

res.not0 <- res %>% dplyr::filter(fold.enrichment>0)
head(res.not0)
maxlogFC <- max(res.not0$logFC)
maxlogFC
minlogFC <- min(res.not0$logFC)
minlogFC


res$logFC <- ifelse(res$logFC==-Inf, minlogFC, res$logFC) 
head(res)
nrow(res)
max(res$logFC)
min(res$logFC)
table(res$MCls, res$contrast)

res_MCls <- unique(res$MCls)
#res_MCls
res_contrast <- unique(res$contrast)
#res_contrast
res_comb <- unique(sort(res$comb))
#res_comb


res_MCls.comb <- rep(unique(res$MCls), each=length(unique(res$contrast)))
#res_MCls.comb
res_contrast.comb <- rep(unique(res$contrast), times=length(unique(res$MCls)))
#res_contrast.comb
res_total.comb <- length(unique(res$contrast))*length(unique(res$MCls))
#res_total.comb
res_dataset_res <- data.frame(MCls=res_MCls.comb, contrast=res_contrast.comb)
#res_dataset_res
res_dataset_res[1,1]



######################################################
# significant motifs
######################################################
#----sig_res_annotation-----
padj_cutoff <- 0.1
logFC_cutoff <- 0.5

sig_res <- res %>% dplyr::filter(qvalue.hyper < padj_cutoff, abs(logFC) > logFC_cutoff)
head(sig_res)
nrow(sig_res)

sig_res_maxlogFC <- max(sig_res$logFC)
sig_res_maxlogFC
sig_res_minlogFC <- min(sig_res$logFC)
sig_res_minlogFC
table(sig_res$MCls, sig_res$contrast)

sig_res_MCls <- unique(sig_res$MCls)
#sig_res_MCls
sig_res_contrast <- unique(sig_res$contrast)
#sig_res_contrast
sig_res_comb <- unique(sort(sig_res$comb))
#sig_res_comb


sig_res_MCls.comb <- rep(unique(sig_res$MCls), each=length(unique(sig_res$contrast)))
#sig_res_MCls.comb
sig_res_contrast.comb <- rep(unique(sig_res$contrast), times=length(unique(sig_res$MCls)))
#sig_res_contrast.comb
sig_res_total.comb <- length(unique(sig_res$contrast))*length(unique(sig_res$MCls))
#sig_res_total.comb
sig_res_dataset <- data.frame(MCls=sig_res_MCls.comb, contrast=sig_res_contrast.comb)
#sig_res_dataset
#sig_res_dataset[1,1]

#write.csv(sig_res, paste0("./1_DiffRNA_2.outs/sig_genes.csv"), quote = FALSE, row.names = FALSE)
#sig_res <- read.csv("./1_DiffRNA_2.outs/sig_genes.csv")






#################################################
#----dataframe & heatmap for only significant----
#################################################
## motifs <- res$motif
## names(motifs) <- res$motif.name
## motifs <- motifs[!duplicated(motifs)]
## head(motifs)
## length(motifs)


sig_motifs <- sig_res$motif
names(sig_motifs) <- sig_res$motif.name
sig_motifs <- sig_motifs[!duplicated(sig_motifs)]
head(sig_motifs)
length(sig_motifs)

res_comb
sig_res_comb

mat <- map_dfc(res_comb, function(ii){
    enrich2 <- res%>%dplyr::filter(comb==ii)
    z <- enrich2$logFC
    #res2 <- res %>% dplyr::filter(comb==ii)
    #b <- res2$estimate
    #names(b) <- res2$peak
    #b[sigdps.qval]
    names(z) <- enrich2$motif.name
    z[names(sig_motifs)]
})
head(mat)
print("Number of significant in mat")
nrow(mat)


mat <- as.matrix(mat)
rownames(mat) <- sig_motifs #DP
colnames(mat) <- res_comb
head(mat)
max(mat)
min(mat)

##mat.mono <- (mat[,"6-Monocyte_vitD"]) %>% as.data.frame() %>%sort()
##mat.mono

#this is the reight way!!
#mat[!is.na(mat)] = 1
#head(mat)
#mat[is.na(mat)] = 1 #since this is a ratio and 1 indicates no change
mat[is.na(mat)] = 0  #if using for logFC
head(mat)
max(mat)
min(mat)


mat.ct <- mat
head(mat.ct)

#----annotation as cell type----
## colnames(mat)
## cvt0 <- str_split(colnames(mat), "_", simplify=T)
## cvt0
## cvt0[,1]

comb.celltype <- colnames(mat.ct)
comb.celltype

cvt.celltype <- as.data.frame(str_split(comb.celltype, "_", simplify=T)) %>% mutate(comb=paste0(V1, "_", V2))
cvt.celltype
#cvt.celltype[,1]

column_ha_ct <- HeatmapAnnotation(
    celltype=cvt.celltype[,1],
    treatment=cvt.celltype[,2],
#    celltype=rep(MCls, each=length(unique(sig_res$contrast))),
#    treatment=rep(contrast, times=length(unique(sig_res$MCls))),
    col=list(celltype=c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
                    "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00", "1-TCM"="pink",
                    "3-TEM"="blue", "5-CD8Naive"="green", "7-dnT"="black"),
             treatment=c("caffeine"="red", "nicotine"="tan",
                         "vitA"="tan4", "vitD"="seagreen4",
                         "vitE"="salmon3", "water"="grey", "zinc"="maroon3")))


#----annotation and data as treatment----
cvt.treat <- cvt.celltype[order(cvt.celltype$V2),]
cvt.treat

mat.treat <- mat.ct[, cvt.treat$comb]
head(mat.treat)

#annotation as treatment
column_ha_treat <- HeatmapAnnotation(
    celltype=cvt.treat[,1],
    treatment=cvt.treat[,2],
    col=list(celltype=c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
                    "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00", "1-TCM"="pink",
                    "3-TEM"="blue", "5-CD8Naive"="green", "7-dnT"="black"),
             treatment=c("caffeine"="red", "nicotine"="tan",
                         "vitA"="tan4", "vitD"="seagreen4",
                         "vitE"="salmon3", "water"="grey", "zinc"="maroon3")))


#----breaks and color----
#breaks <- c(seq(0,1,length.out=50), quantile(b2, probs=seq(0, 1, length.out=49)), 38.93) 
#breaks
#col_fun <-  colorRamp2(breaks,
#   colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100))


#breaks <- (seq(-6, 6, length.out=13)) 
#breaks


#col_fun <-  colorRamp2(breaks,
#   colorRampPalette(rev(brewer.pal(n=11, name="RdBu")))(13))


#col_fun <-  colorRampPalette(brewer.pal(n=9, name="Reds"))(100)

y <- as.numeric(mat.treat)
max(y)
abs(min(y))
#length(y)

scale <- quantile(abs(y),probs=0.99)
scale


mybreaks <- seq(-scale, scale, length.out=100)
mybreaks
names(mybreaks) <- NULL

col_fun <- colorRampPalette(rev(brewer.pal(n=11, name="RdBu")))(100)


#----figure as celltype----
fig <- Heatmap(mat.ct,
    col_fun,
#   Rowv = NA,
#   Colv = NA, 
   cluster_rows=T, cluster_columns=F,
   top_annotation=column_ha_ct,
   heatmap_legend_param=list(title="fold.enrichment",
      title_gp=gpar(fontsize=10),
      labels_gp=gpar(fontsize=10)),
   show_row_names=T, show_column_names=T,
   row_names_gp=gpar(fontsize=2),
   column_names_gp=gpar(fontsize=10),
   raster_device="png")


figfn <- "./1_motif.outs/Figure4.2_significant_heatmap_macs2_0.1_cn_celltype_allcols_control.png"
png(figfn, height=5000, width=3500, res=300)
#set.seed(0)
fig <- draw(fig)
dev.off()


#----figure as treat----
fig <- Heatmap(mat.treat,
    col_fun,
#   Rowv = NA,
#   Colv = NA, 
   cluster_rows=T, cluster_columns=F,
   top_annotation=column_ha_treat,
   heatmap_legend_param=list(title="fold.enrichment",
      title_gp=gpar(fontsize=10),
      labels_gp=gpar(fontsize=10)),
   show_row_names=T, show_column_names=T,
   row_names_gp=gpar(fontsize=2),
   column_names_gp=gpar(fontsize=10),
   raster_device="png")


figfn <- "./1_motif.outs/Figure4.2_significant_heatmap_macs2_0.1_cn_treat_allcols_control.png"
png(figfn, height=5000, width=3500, res=300)
#set.seed(0)
fig <- draw(fig)
dev.off()


#######################################################
#----dataframe & heatmap for only significant top10----
#######################################################
#sig_res_comb2 <- "6-Monocyte_vitD"
sig_res_filt <- map_dfr(sig_res_comb, function(ii){
    enrich2 <- sig_res %>% dplyr::filter(comb==ii)
    enrich2 <- enrich2 %>% dplyr::arrange(qvalue.fisher)
    #enrich2 <- enrich2 %>% dplyr::arrange(pvalue)
    enrich2 <- enrich2[1:10,]
})
head(sig_res_filt)
nrow(sig_res_filt)

#colnames(sig_res_filt)
#sig_res_filt %>% dplyr::select("motif.name", "qvalue.fisher", "observed", "logFC")

sig_motifs <- sig_res_filt$motif
names(sig_motifs) <- sig_res_filt$motif.name
sig_motifs <- sig_motifs[!duplicated(sig_motifs)]
head(sig_motifs)
length(sig_motifs)

mat <- map_dfc(res_comb, function(ii){
    enrich2 <- res %>% dplyr::filter(comb==ii)
    z <- enrich2$logFC
    names(z) <- enrich2$motif.name
    z[names(sig_motifs)]
})
head(mat)
length(mat[[1]])
#min(mat)
#max(mat)


#****
mat <- do.call(cbind, mat)
colnames(mat) <- res_comb
rownames(mat) <- names(sig_motifs)
head(mat)
length(mat)
nrow(mat)
max(mat)
min(mat)

mat2 <- mat

#this is the reight way!!
#mat[!is.na(mat)] = 1
#head(mat)
#mat[is.na(mat)] = 1 #since this is a ratio and 1 indicates no change
mat2[is.na(mat2)] = 0  #if using for logFC
head(mat2)
nrow(mat2)
max(mat2)
min(mat2)


#mat.df <- mat2 %>% as.data.frame()
#head(mat.df)
#mat.mono <- mat.df %>% dplyr::select("6-Monocyte_vitD") %>% mutate(logFC = mat.mono$"6-Monocyte_vitD")
#head(mat.mono)
#mat.mono2 <- mat.mono[order(mat.mono$logFC, decreasing=TRUE),]
#head(mat.mono2)





mat.ct <- mat2
head(mat.ct)

#----annotation as cell type----
## colnames(mat)
## cvt0 <- str_split(colnames(mat), "_", simplify=T)
## cvt0
## cvt0[,1]

comb.celltype <- colnames(mat.ct)
comb.celltype

cvt.celltype <- as.data.frame(str_split(comb.celltype, "_", simplify=T)) %>% mutate(comb=paste0(V1, "_", V2))
cvt.celltype
#cvt.celltype[,1]

column_ha_ct <- HeatmapAnnotation(
    celltype=cvt.celltype[,1],
    treatment=cvt.celltype[,2],
#    celltype=rep(MCls, each=length(unique(sig_res$contrast))),
#    treatment=rep(contrast, times=length(unique(sig_res$MCls))),
    col=list(celltype=c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
                    "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00", "1-TCM"="pink",
                    "3-TEM"="blue", "5-CD8Naive"="green", "7-dnT"="black"),
             treatment=c("caffeine"="red", "nicotine"="tan",
                         "vitA"="tan4", "vitD"="seagreen4",
                         "vitE"="salmon3", "water"="grey", "zinc"="maroon3")))


#----annotation and dataas treatment----
cvt.treat <- cvt.celltype[order(cvt.celltype$V2),]
cvt.treat

mat.treat <- mat.ct[, cvt.treat$comb]
head(mat.treat)

#annotation as treatment
column_ha_treat <- HeatmapAnnotation(
    celltype=cvt.treat[,1],
    treatment=cvt.treat[,2],
    col=list(celltype=c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
                    "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00", "1-TCM"="pink",
                    "3-TEM"="blue", "5-CD8Naive"="green", "7-dnT"="black"),
             treatment=c("caffeine"="red", "nicotine"="tan",
                         "vitA"="tan4", "vitD"="seagreen4",
                         "vitE"="salmon3", "water"="grey", "zinc"="maroon3")))


#----breaks and color----
#breaks <- c(seq(0,1,length.out=50), quantile(b2, probs=seq(0, 1, length.out=49)), 38.93) 
#breaks
#col_fun <-  colorRamp2(breaks,
#   colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100))


#breaks <- (seq(-6, 6, length.out=13)) 
#breaks


#col_fun <-  colorRamp2(breaks,
#   colorRampPalette(rev(brewer.pal(n=11, name="RdBu")))(13))

##col_fun <-  colorRampPalette(brewer.pal(n=9, name="Reds"))(100)


y <- as.numeric(mat.ct)
length(y)
max(y)
min(y)


scale <- quantile(abs(y),probs=0.99)
scale
length(y[abs(y)>scale])

mybreaks <- seq(-scale, scale, length.out=100)
mybreaks
names(mybreaks) <- NULL

col_fun <- colorRampPalette(rev(brewer.pal(n=11, name="RdBu")))(100)



#----figure as celltype----
fig <- Heatmap(mat.ct,
    col_fun,
#   Rowv = NA,
#   Colv = NA, 
   cluster_rows=T, cluster_columns=F,
   top_annotation=column_ha_ct,
   heatmap_legend_param=list(title="fold.enrichment",
      title_gp=gpar(fontsize=14),
      labels_gp=gpar(fontsize=14)),
   show_row_names=T, show_column_names=T,
   row_names_gp=gpar(fontsize=7),
   column_names_gp=gpar(fontsize=14),
   raster_device="png")


figfn <- "./1_motif.outs/Figure5.1_top10_heatmap_macs2_0.1_cn_celltype_allcols_control.png"
png(figfn, height=5000, width=3500, res=300)
#set.seed(0)
fig <- draw(fig)
dev.off()


#----figure as treat----
fig <- Heatmap(mat.treat,
    col_fun,
#   Rowv = NA,
#   Colv = NA, 
   cluster_rows=T, cluster_columns=F,
   top_annotation=column_ha_treat,
   heatmap_legend_param=list(title="fold.enrichment",
      title_gp=gpar(fontsize=14),
      labels_gp=gpar(fontsize=14)),
   show_row_names=T, show_column_names=T,
   row_names_gp=gpar(fontsize=7),
   column_names_gp=gpar(fontsize=14),
   raster_device="png")


figfn <- "./1_motif.outs/Figure5.1_top10_heatmap_macs2_0.1_cn_treat_allcols_control.png"
png(figfn, height=5000, width=3500, res=350)
#set.seed(0)
fig <- draw(fig)
dev.off()

   

## enrich%>%mutate(LFC=log2(fold.enrichment))%>%
##    dplyr::filter(qvalue.hyper<0.1, LFC>0.5)%>%
##    group_by(condition)%>%summarise(ny=n(),.groups="drop")
## head(enrich)





################
### qq plots ###
################
dfNew <- map_dfr(res_comb, function(ii){
  res2 <- res%>%dplyr::filter(comb==ii)
  ngene <- nrow(res2)
  res2 <- res2%>%
      arrange(pvalue)%>%
      mutate(observed=-log10(pvalue), expected=-log10(ppoints(ngene)))
  res2
})

head(dfNew)
table(dfNew$contrast)

lab1 <- c("caffeine"="caffeine", "nicotine"="nicotine", "zinc"="zinc",
          "vitA"="vitA", "vitD"="vitD", "vitE"="vitE", "water"="water")

## lab2 <- c("Bcell"="B cell", "Monocyte"= "Monocyte", "NKcell"="NK cell",
##    "Tcell"="T cell")

lab2 <- c("4_Bcell"="4_Bcell", "6_Monocyte"="6_Monocyte",
          "2_NKcell"="2_NKcell", "0_CD4Naive"="0_CD4Naive", "1_TCM"="1_TCM",
          "3_TEM"="3_TEM", "5_CD8Naive"="5_CD8Naive", "7_dnT"="7_dnT")

p <- ggplot(dfNew, aes(x=expected,y=observed))+
   ggrastr::rasterise(geom_point(size=0.3, color="grey40"),dpi=300)+
   geom_abline(colour="red")+
   facet_grid(MCls~contrast, scales="free",
      labeller=labeller(contrast=lab1, MCls=lab2))+
   xlab(bquote("Expected"~ -log[10]~"("~italic(plain(P))~")"))+
   ylab(bquote("Observed"~-log[10]~"("~italic(plain(P))~")"))+
   theme_bw()+
   theme(strip.text=element_text(size=12))

figfn <- "./1_motif.outs/Figure2.2_qq_macs2_0.1_cn_control.png"
png(figfn, width=1600, height=1600, res=175)
print(p)
dev.off()






#------------------plot_loop--------------------------------
#nrow(res) #if qq plot has been called before this step, this already is with dropped NA values
res.1 <- dfNew %>% drop_na(pvalue)
head(res.1)

## res.1.not0 <- res %>% dplyr::filter(fold.enrichment>0)
## head(res.not0)
## maxlogFC <- max(res.not0$logFC)
## maxlogFC
## minlogFC <- min(res.not0$logFC)
## minlogFC
## ## res %>% dplyr::filter(logFC < -3.93, logFC != -Inf)

## res$logFC <- ifelse(res$logFC==-Inf, minlogFC, res$logFC)
## head(res)
## nrow(res)
## max(res$logFC)
## min(res$logFC)

## res.1 <- res.1 %>% dplyr::filter ifelse(res$logFC==-Inf, minlogFC, res$logFC) 

## max(res.1$pvalue)
## max(res.1 %>% dropres.1$observed)

## max.y <- floor(max(-log10(res.1$pvalue)))/2
## max.y

png(paste0("./1_motif.outs/Figure2.3_qq_macs2_0.1_cn_clusters_together_control.png"), width=2000, height=4000, pointsize=22, res=200)
par(mar=c(4,4,2,2),mgp=c(2,1,0))
x <- matrix(1:8, 4, 2, byrow=T)
layout(x)
MCls <- unique(dfNew$MCls)

for (oneMCl in MCls){
        ##1
        res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="caffeine")
        plot(res.2$expected, res.2$observed,
             main=oneMCl,
             pch = 19,
             #ylim = c(0, max.y),
             ylab="observed -log10(p)",
             xlab="expected -log10(p)",
             col="red")
        ##2
        res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="nicotine")
        points(res.2$expected, res.2$observed,  col="tan", pch = 19)
        ##3
        res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="zinc")
        points(res.2$expected, res.2$observed,  col="maroon3", pch = 19)
        ##4
        res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="vitA")
        points(res.2$expected, res.2$observed,  col="tan4", pch = 19)
        ##5
        res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="vitD")
        points(res.2$expected, res.2$observed,  col="seagreen4", pch = 19)
        ##6
        ##res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="water")
        ##points(res.2$expected, res.2$observed,  col="grey", pch = 19)
        ##7
        res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="vitE")
        points(res.2$expected, res.2$observed,  col="salmon3", pch = 19)
        ##control
        abline(0, 1, col = "black")
        ##legend
        legend("topleft",
               ##inset=c(-0.2,0),
               cex = 1.2,
               pch = 19,
               ##c("caff","nic", "zinc", "vitA", "vitD", "vitE", "water"),
               ##fill=c("red","tan", "maroon3", "tan4", "seagreen4", "salmon3", "grey"))
               c("caff","nic", "zinc", "vitA", "vitD", "vitE"),
               fill=c("red","tan", "maroon3", "tan4", "seagreen4", "salmon3"))
}
dev.off()










############################################################################
############################################################################
### motif enrichment analysis directionally, up and down ###
############################################################################
############################################################################
#----MOTIF FILE SAVED IN THIS SCRIPT----
#fn <- "./1_motif.outs/1_scATAC.motif.rds" 
#atac <- read_rds(fn)
## motif <- Motifs(atac)
## pfm <- GetMotifData(object=motif, slot="pwm")
## motif <- SetMotifData(object=motif, slot="pwm", new.data=pfm)

#----RESDP FILE----
## differential peaks

fn <- "../2_Differential/1.2_DiffPeak.outs/2.0_DESeq.results_0.1_cn_sampleID+new_treat.rds"
resDP <- read_rds(fn) %>% as.data.frame() %>% mutate(direction=ifelse(estimate>0, 1, 0))

resDP <- resDP %>% as.data.frame() %>% mutate(direction=ifelse(estimate>0, 1, 0))
head(resDP)

### fisher test
cal.fisher <- function(df){
  ###  
  resfisher <- map_dfr(1:nrow(df),function(i){  
     dmat <- matrix(as.numeric(df[i,]),2,2)
     colnames(dmat) <- c("interest", "not.interest")
     rownames(dmat) <- c("in.motif", "not.motif")
     res <- fisher.test(dmat, alternative="greater")
     res2 <- data.frame(odds=res$estimate, pval.fisher=res$p.value)
     res2
  })
  resfisher
}

    
#MCls <- rep(c("Bcell", "Monocyte", "NKcell", "Tcell"), each=4)
#contrast <- rep(c("LPS", "LPS-DEX", "PHA", "PHA-DEX"), times=4)

MCls <- rep(unique(resDP$MCls), each=length(unique(resDP$contrast)))
#MCls
contrast <- rep(unique(resDP$contrast), times=length(unique(resDP$MCls)))
#contrast
total.comb <- length(unique(resDP$contrast))*length(unique(resDP$MCls))
#total.comb
dataset <- data.frame(MCls=MCls, contrast=contrast)
head(dataset)
dataset[1,1]




###enrichment motif analysis
enriched.motif <- lapply(1:total.comb, function(i){
###
    cell0 <- dataset[i,1]
    contrast0 <- dataset[i,2]
    cat(i, cell0, contrast0, "\n") 
###
    enrich <- lapply(c(0,1),function(ii){ 
        res2 <- resDP%>%
            dplyr::filter(MCls==cell0, contrast==contrast0)
        top.DP <- res2%>%
            dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5, direction==ii)%>%
            dplyr::pull(gene)%>%as.character()
        bg.DP <- res2%>%dplyr::pull(gene)%>%as.character() 
        n.interest <- length(top.DP)
        n.not <- length(bg.DP)-n.interest
### 
        if(n.interest>5){
            enrich2 <- FindMotifs(
                object=atac,
                features=top.DP,
                background=bg.DP)
            df <- data.frame("interest.in.motif"=enrich2$observed,
                             "interest.not.motif"=n.interest-enrich2$observed,
                             "not.interest.in.motif"=enrich2$background-enrich2$observed)
            df <- df%>%mutate("not.interest.not.motif"=n.not-not.interest.in.motif)
            fisher <- cal.fisher(df)
            enrich2 <- cbind(enrich2, fisher)
            enrich2$qvalue.hyper <- p.adjust(enrich2$pvalue,"BH")
            enrich2$qvalue.fisher <- p.adjust(enrich2$pval.fisher,"BH")
            enrich2$MCls <- cell0
            enrich2$contrast <- contrast0
            enrich2$direction <- ii
        }else{
            enrich2 <- NA   
        }
        enrich2
    })    
   ##return results 
    if ( sum(is.na(enrich))==2){
        enrich <- NA
    }else{   
        enrich <-enrich[!is.na(enrich)]
        enrich <- do.call(rbind,enrich)
    }   
    enrich 
})

head(enriched.motif)

enriched.motif <- enriched.motif[!is.na(enriched.motif)]
enriched.motif <- do.call(rbind,enriched.motif)
head(enriched.motif)

opfn <- "./1_motif.outs/3_motif.enrich.direction_macs2_0.1_cn_control.rds"
write_rds(enriched.motif, opfn)

enriched.motif <- read_rds("./1_motif.outs/3_motif.enrich.direction_macs2_0.1_cn.rds")
head(enriched.motif)



######################
###    heatmap     ###
######################
rm(list=ls())

direction2 <- c("0"="Down","1"="Up")

enrich <- read_rds("./1_motif.outs/3_motif.enrich.direction_macs2_0.1_cn.rds")
#enrich <- enriched.motif

enrich <- enrich%>%
   drop_na(fold.enrichment)%>%
   mutate(dir2=direction2[as.character(direction)],
          condition=paste(MCls, contrast, dir2, sep="_"))

head(enrich)

condition <- unique(enrich$condition)
condition

#----build matrix for heatmap----
motif <- enrich$motif
names(motif) <- enrich$motif.name
motif <- motif[!duplicated(motif)]


###
mat <- lapply(condition, function(ii){
   enrich2 <- enrich%>%dplyr::filter(condition==ii)
   z <- enrich2$fold.enrichment
#   names(z) <- enrich2$motif
#   z[motif]
   names(z) <- enrich2$motif.name
   z[names(motif)]
})

mat <- do.call(cbind, mat)
colnames(mat) <- condition
#head(mat)

#----not using this step----
#rownames(mat) <- motif
#----not using this step----

###
b <- as.vector(mat)

b2 <- b[b>1&b<2.62]

breaks <- c(seq(0, 1, length.out=50),
            quantile(b2, probs=seq(0, 1, length.out=49)),7.71) 

col_fun <-  colorRamp2(breaks,
   colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100))

## column_ha <- HeatmapAnnotation(
##     celltype=gsub("_.*", "", condition),
##     treatment=gsub(".*cell_|.*cyte_|_(D|U).*", "", condition),
##     col=list(celltype=c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
##                     "NKcell"="#aa4b56", "Tcell"="#ffaa00"),
##              treatment=c("LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
##                         "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")))

MCls <- unique(enrich$MCls)
MCls

treats <- unique(enrich$contrast)
treats

#done above
## condition <- unique(enrich$condition)
## condition

celltype <- str_split(condition, "_", simplify=T)
celltype

celltype.2 <- paste(celltype[,1], "_", celltype[,2], sep="")
celltype.2

#treatment=gsub(".*cell_|.*cyte_|_(D|U).*", "", condition)
#treatment
treatment.2 <- celltype[,3]
treatment.2

column_ha <- HeatmapAnnotation(
    celltype=celltype.2,
    treatment=treatment.2,
    col=list(celltype=c("4_Bcell"="#4daf4a", "6_Monocyte"="#984ea3",
                    "2_NKcell"="#aa4b56", "0_CD4Naive"="#ffaa00", "1_TCM"="pink",
                    "3_TEM"="blue", "5_CD8Naive"="green", "7_dnT"="black"),
             treatment=c("caffeine"="red", "nicotine"="tan",
                         "vitA"="tan4", "vitD"="seagreen4",
                         "vitE"="salmon3", "water"="grey", "zinc"="maroon3"))
)


## fig <- Heatmap(mat, col=col_fun,
##    cluster_rows=T, cluster_columns=F,
##    top_annotation=column_ha,
##    heatmap_legend_param=list(title="fold.enrichment",
##       title_gp=gpar(fontsize=10),
##       labels_gp=gpar(fontsize=10)),
##    show_row_names=F, show_column_names=T,
##    column_names_gp=gpar(fontsize=5),
##    raster_device="png")

## figfn <- "./1_motif.outs/Figure3.1_heatmap_macs2_0.1_cn.png"
## png(figfn, height=800, width=700, res=120)
## set.seed(0)
## fig <- draw(fig)
## dev.off()

fig <- Heatmap(mat,
   cluster_rows=T, cluster_columns=F,
   top_annotation=column_ha,
   heatmap_legend_param=list(title="fold.enrichment",
      title_gp=gpar(fontsize=10),
      labels_gp=gpar(fontsize=10)),
   show_row_names=F, show_column_names=T,
   column_names_gp=gpar(fontsize=5),
   raster_device="png")

figfn <- "./1_motif.outs/Figure3.1_heatmap_macs2_0.1_cn_2.png"
png(figfn, height=3200, width=2800, res=225)
#set.seed(0)
fig <- draw(fig)
dev.off()



enrich%>%mutate(LFC=log2(fold.enrichment))%>%
   dplyr::filter(qvalue.hyper<0.1, LFC>0.5)%>%
   group_by(condition)%>%summarise(ny=n(),.groups="drop")




##############################
### barplot enriched motif ###
##############################
direction2 <- c("0"="Down","1"="Up")

enrich <- read_rds("./1_motif.outs/3_motif.enrich.direction_macs2_0.1_cn_control.rds")
head(enrich)
nrow(enrich)
##table(enrich$MCls)
##table(enrich$contrast)

enrich <- enrich%>%
   drop_na(fold.enrichment)%>%
   mutate(LFC=log2(fold.enrichment))

###
res <- enrich%>%dplyr::filter(qvalue.hyper<0.1, fold.enrichment>1.41)

sigs <- res%>%group_by(MCls, contrast, direction) %>% summarise(ny=n(), .groups="drop")
head(sigs)

###
sig4 <- sigs%>%mutate(ny2=ifelse(direction==0, -ny, ny))
sig4

breaks_value <- pretty(c(-600, 600), 5)
breaks_value

## p <- ggplot(sig4, aes(x=MCls, y=ny2))+
##    geom_bar(aes(fill=factor(MCls), alpha=factor(direction)), stat="identity")+
##    scale_fill_manual(values=c("Bcell"="#4daf4a",
##       "Monocyte"="#984ea3", "NKcell"="#aa4b56", "Tcell"="#ffaa00"))+
##    scale_alpha_manual(values=c("0"=0.5, "1"=1))+
##    geom_hline(yintercept=0, color="grey60")+
##    geom_text(aes(x=MCls, y=ny2, label=abs(ny2),
##       vjust=ifelse(direction==1, -0.2, 1.2)), size=3)+
##    scale_y_continuous("", breaks=breaks_value, limits=c(-600,650),
##                       labels=abs(breaks_value))+
##    facet_grid(~contrast,
##       labeller=labeller(contrast=c("LPS"="LPS","LPS-DEX"="LPS+DEX",
##                                    "PHA"="PHA", "PHA-DEX"="PHA+DEX")))+
##    theme_bw()+
##    theme(legend.position="none",
##          axis.title.x=element_blank(),
##          axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5))


p <- ggplot(sig4, aes(x=MCls, y=ny2))+
    geom_bar(aes(fill=factor(MCls), alpha=factor(direction)), stat="identity")+
    scale_fill_manual(values=c("4_Bcell"="#4daf4a", "6_Monocyte"="#984ea3",
                               "2_NKcell"="#aa4b56", "0_CD4Naive"="#ffaa00",
                               "1_TCM"="#ffaa00", "3_TEM"="#ffaa00",
                               "5_CD8Naive"="#ffaa00", "7_dnT"="#ffaa00"))+
    scale_alpha_manual(values=c("0"=0.5, "1"=1))+
    geom_hline(yintercept=0, color="grey60")+
    geom_text(aes(x=MCls, y=ny2, label=abs(ny2),
                  vjust=ifelse(direction==1, -0.2, 1.2)), size=3)+
    scale_y_continuous("", breaks=breaks_value, limits=c(-600,650),
                       labels=abs(breaks_value))+
    facet_grid(~contrast,
               ## labeller=labeller(contrast=c("4_Bcell"="4_Bcell",
               ##                              "6_Monocyte"="6_Monocyte",
               ##                              "2_NKcell"="2_NKcell",
               ##                              "0_CD4Naive"="0_CD4Naive",
               ##                              "1_TCM"="1_TCM",
               ##                              "3_TEM"="3_TEM",
               ##                              "5_CD8Naive"="5_CD8Naive",
               ##                              "7_dnT"="7_dnT")))+
               labeller=labeller(contrast=c("caffeine"="caffeine",
                                            "nicotine"="nicotine",
                                            "zinc"="zinc",
                                            "vitA"="vitA",
                                            "vitD"="vitD",
                                            "vitE"="vitE")))+
    theme_bw()+
    theme(legend.position="none",
          axis.title.x=element_blank(),
          axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5))


###
figfn <- "./1_motif.outs/Figure3.2_barplot_macs2_0.1_cn_control.png"
png(filename=figfn, width=3000, height=500, pointsize=24, res=225)
print(p)
dev.off()
                      







#######################################
### dot plots show motif enrichment ###
#######################################
#----PREPARE DATA FOR DOT PLOT----
#----ENRICHED.MOTIF----
enrich <- read_rds("./1_motif.outs/3_motif.enrich.direction_macs2_0.1_cn.rds")
head(enrich)

#### label
## clst <- paste(rep(c("LPS", "PHA", "LPS-DEX", "PHA-DEX"), each=4),
##    rep(c("Bcell", "Monocyte", "NKcell", "Tcell"), times=4), sep=".")


MCls <- unique(enrich$MCls)
MCls
length(MCls)

treats <- unique(enrich$contrast)
treats
length(treats)

## enrich <- enrich%>%
##    drop_na(fold.enrichment)%>%
##    mutate(dir2=direction2[as.character(direction)],
##           condition=paste(MCls, contrast, dir2, sep="_"))

#condition <- unique(enrich$condition)
#condition


clst <- paste(rep(unique(enrich$contrast), each=length(unique(enrich$MCls))),
   rep(unique(enrich$MCls), times=length(unique(enrich$contrast))), sep=".")
clst
length(clst)

clst <- paste(rep(clst,times=2), rep(c("Up","Down"), each=length(clst)), sep=".")
clst
length(clst)


###
ii <- paste(rep(c("A","B","C","D", "E", "F", "G"),each=length(unique(enrich$MCls))),rep(1:length(unique(enrich$MCls)),times=length(unique(enrich$contrast))), sep="")
ii

clst2 <- paste(rep(c("X","Y"),each=length(clst)/2), rep(ii,times=2), sep=".")
clst2

names(clst2) <- clst
clst2

lab2 <- gsub("-", "+", clst)
lab2

names(lab2) <- clst2 
lab2


###
drt2 <- c("0"="Down", "1"="Up")

enrich <- enrich%>%mutate(direction2=drt2[as.character(direction)],
   cluster=paste(contrast, MCls, direction2, sep="."),
   newCluster=clst2[cluster])
head(enrich)


###
## topmotif <- enrich%>%
##    dplyr::filter(fold.enrichment>1.41, qvalue.fisher<0.1)%>%dplyr::pull(motif)%>%unique()
##    ## group_by(cluster)%>%
##    ## top_n(n=5, wt=fold.enrichment)%>%ungroup()%>%as.data.frame()
## head(topmotif)
## length(topmotif)

## head(enrich2)
## topmotif <- enrich2%>%dplyr::pull(motif)
## head(topmotif)

## enrich3 <- enrich%>%
##    dplyr::filter(fold.enrichment>1.41, qvalue.fisher<0.1)
## head(enrich3)
## nrow(enrich3)


topmotif <- enrich%>%
       ## dplyr::filter(fold.enrichment>1.41, qvalue.fisher<0.1)%>%
       ## dplyr::pull(motif)%>%unique()
       group_by(cluster)%>%
       top_n(n=6, wt=fold.enrichment)%>%ungroup()%>%dplyr::pull(motif)%>%unique()
head(topmotif)

enrich3 <- enrich%>%
   dplyr::filter(motif%in%topmotif, fold.enrichment>1.41, qvalue.fisher<0.1)
head(enrich3)
nrow(enrich3)




#----MOTIF ACTIVITIES/CHROMVAR----
res <- read_rds("./2_motif.activities.outs/3_motif.diff.results_macs2_0.1_cn.rds")
head(res)


#motif <- read_rds("./2_motif.activities.outs/2_motif.ave_macs2_0.1_cn.rds")
#head(motif)

topmotif2 <- res%>%dplyr::filter(qval<0.1, abs(beta)>1.41)%>%dplyr::pull(motif)%>%unique()
head(topmotif2)
length(topmotif2)





#----PLOT----
p <- ggplot(enrich3, aes(x=newCluster, y=motif.name))+
   geom_point(aes(size=fold.enrichment, colour=qvalue.fisher))+
   scale_colour_gradient(name="FDR",
      low="blue", high="red", na.value=NA, trans="reverse", n.breaks=5,
      guide=guide_colourbar(order=1))+
   scale_size_binned("fold enrichment",
      guide=guide_bins(show.limits=T, axis=T,
          axis.show=arrow(length=unit(1.5,"mm"), ends="both"), order=2),
      n.breaks=4)+
   scale_x_discrete(labels=lab2)+ 
   theme_bw()+
   theme(axis.title=element_blank(),
         axis.text.x=element_text(angle=60, hjust=1, size=10),
         axis.text.y=element_text(size=8))

###
figfn <- "./1_motif.outs/Figure3.3_dotplot_macs2_0.1_cn.png"
png(figfn, width=4000, height=4000, res=250)
print(p)
dev.off()




################
### qq plots ###
################
drt2 <- c("0"="Down", "1"="Up")

res <- read_rds("./1_motif.outs/3_motif.enrich.direction_macs2_0.1_cn.rds")%>%
   as.data.frame()%>%
   drop_na(pvalue)%>%
   mutate(direction2=drt2[as.character(direction)],
          comb=paste(MCls, contrast, direction2, sep="_"))
head(res)

comb <- sort(unique(res$comb))
comb

dfNew <- map_dfr(comb, function(ii){
  res2 <- res%>%dplyr::filter(comb==ii)
  ngene <- nrow(res2)
  res2 <- res2%>%
     arrange(pvalue)%>%
      mutate(observed=-log10(pvalue), expected=-log10(ppoints(ngene)))
  res2
})

head(dfNew)

unique(dfNew$MCls)

unique(dfNew$contrast)

unique(dfNew$comb)

vars(comb)

p <- ggplot(dfNew, aes(x=expected,y=observed))+
   ggrastr::rasterise(geom_point(size=0.3, color="grey40"),dpi=300)+
   geom_abline(colour="red")+
   facet_wrap(vars(comb), scales="free", nrow=112, ncol=7)+
   xlab(bquote("Expected"~ -log[10]~"("~italic(plain(P))~")"))+
   ylab(bquote("Observed"~-log[10]~"("~italic(plain(P))~")"))+
   theme_bw()+
   theme(axis.text=element_text(size=7),
         strip.text=element_text(size=8))

figfn <- "./1_motif.outs/Figure3.4_qq_macs2_0.1_cn.png"
png(figfn, width=3200, height=3600, res=225)
print(p)
dev.off()




