##
library(tidyverse)
library(Seurat)
library(parallel)
library(SeuratDisk)
library(SeuratData)
library(SeuratObject)
library(Signac)
library(SeuratWrappers)
library(SeuratData)
library(chromVAR)
library(JASPAR2020)
library(JASPAR2022)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38)
##library(BSgenome.Hsapiens.1000genomes.hs37d5, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
##library(BSgenome.Hsapiens.1000genomes.hs37d5)
##library(ChIPseeker, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(ChIPseeker)
###
library(ggplot2)
##library(cowplot, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(cowplot)
library(grid)
library(gridExtra)
library(ggExtra)
library(RColorBrewer)
library(ggsci)
library(viridis)
library(ComplexHeatmap)
library(circlize)
theme_set(theme_grey())
library(reshape2)
library("ggpubr")


###
outdir <- "./2_motif.activities_2.outs/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


##################################################
#---- chromvar object----
##################################################
atac <- read_rds("./2_motif.activities.outs/1_scATAC.motifActivities.rds")
atac




## #####################################################
## #----differential motif activities----
## #####################################################
## #YtX_sel
## mat <- read_rds("./2_motif.activities.outs/2_motif.ave_macs2_0.1_cn_control.rds")
## head(mat)
## ncol(mat)

## meta <- str_split(colnames(mat), "_", simplify=T)%>%as.data.frame()
## head(meta)

## names(meta) <- c("cluster", "MCls", "treat", "sampleID")
## head(meta)
## unique(meta$treat)


## #plot
## ## b <- as.vector(mat)
## ## breaks <- quantile(b, probs=seq(0, 1, length.out=100))
## ## col_fun <-  colorRamp2(breaks,
## ##                           colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100))


## y <- as.numeric(mat)
## max(y)
## min(y)
## scale <- quantile(abs(y),probs=0.99)
## scale
## mybreaks <- seq(-scale, scale, length.out=100)
## names(mybreaks) <- NULL
## mybreaks

## mycol <- colorRampPalette(rev(brewer.pal(n=11, name="RdBu")))(100)


## ##
## column_ha <- HeatmapAnnotation(
##     celltype=paste0(meta$cluster, "_", meta$MCls),
##     treatment=meta$treat,
##     col=list(celltype=c("4_Bcell"="#4daf4a", "6_Monocyte"="#984ea3",
##                         "2_NKcell"="#aa4b56", "0_CD4Naive"="#ffaa00",
##                         "1_TCM"="pink","3_TEM"="blue",
##                         "5_CD8Naive"="green", "7_dnT"="black", "8_MAIT"="grey",
##                         "9_Platelet"="purple4", "10_DC"="slategray4"),
##              treatment=c("caffeine"="red", "nicotine"="tan",
##                          "vitA"="tan4", "vitD"="seagreen4",
##                          "vitE"="salmon3", "zinc"="maroon3",
##                          "control"="grey")))
##                          ##"water"="grey", "etOH"="darkgrey")))

## fig <- Heatmap(mat,
##                ##col=col_fun,
##                cluster_rows=T, cluster_columns=F,
##                show_row_dend=T, show_column_dend=F,
##                top_annotation=column_ha,
##                heatmap_legend_param=list(title="motif activities",
##                                          title_gp=gpar(fontsize=10),
##                                          labels_gp=gpar(fontsize=10)),
##                show_row_names=F, show_column_names=F,
##                use_raster=F, raster_device="png")

## figfn <- "./2_motif.activities.outs/Figure1.1_heatmap.motif.activities_macs2_0.1_cn_2_nocluster_control.png"
## png(figfn, height=800, width=900, res=120)
## #set.seed(0)
## fig <- draw(fig)
## dev.off()








#######################################
### Differential activities for ind ###
#######################################
#res <- read_rds("./2_motif.activities.outs/3_motif.diff.results_macs2_0.1_cn.rds")
#res <- read_rds("./2_motif.activities.outs/3_motif.diff.results_macs2_0.1_cn_all_cols.rds") #include all values like coeff treatment, coeff control, etc
##res <- read_rds("./2_motif.activities.outs/3_motif.diff.results_macs2_0.1_cn_indfilt_all_cols.rds")


res <- read_rds("./2_motif.activities.outs/3_motif.diff.results_macs2_0.1_cn_indfilt_all_cols.rds")


res <- res %>% mutate(condition=paste0(gsub("_", "-", MCls), "_", contrast))
#res <- res%>%mutate(condition=paste(MCls, contrast, sep="_"))
res <- res %>% mutate(zscore=beta/stderr)

res.filt <- res%>%dplyr::filter(qval<0.1, abs(beta)>0)
topmotif <- sort(unique(res.filt$motif))
print("Unique significant motifs are:")
length(topmotif)

sig_motifs_table <- table(res.filt$MCls, res.filt$contrast)
sig.matrix <-as.matrix(sig_motifs_table) # create a numeric matrix object
sig.melt <- melt(sig.matrix)
sig.p <- ggplot(sig.melt, aes(x = Var2, y = Var1)) +
      geom_tile(aes(fill=value), color="black")+
      scale_fill_gradient(low = "white", high = "white")+
      geom_text(aes(label = value), color = "black", size = 5) +
      ggtitle(paste0("Significant_motifs_coef_0")) +
      theme(axis.text.x = element_text(angle = 90),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    legend.position="none")
figure <- ggarrange(sig.p + font("x.text", size = 14))
png(paste0("./2_motif.activities_2.outs/sig_motifs_table_indfilt_coef_0.png"), width=500, height=500, pointsize=16, res=125)
print(figure)
dev.off()

condition <- unique(sort(res$condition))
condition


##--------------------------------
##----mat data for significant----
##--------------------------------
head(res)

condition <- unique(res$condition)
condition

mat.sigmotifs.qval <- map_dfr(condition, function(ii){
      res2 <- res.filt %>% dplyr::filter(condition==ii)
      })
sigmotif.qval <- unique(mat.sigmotifs.qval$motif)
mat.sigmotifs.qval.2 <- map_dfc(condition, function(ii){
      res2 <- res%>%dplyr::filter(condition==ii)
        b <- res2$zscore
        names(b) <- res2$motif
        b[sigmotif.qval]
      })

mat <- mat.sigmotifs.qval.2
mat <- as.matrix(mat)
mat[is.na(mat)] = 0
colnames(mat) <- condition
rownames(mat) <- sigmotif.qval
head(mat)
print("unique sig motifs")
nrow(mat.sigmotifs.qval.2)

mat2 <- mat

x <- as.data.frame(str_split(condition, "_", simplify=T))
x <- x %>% dplyr::arrange(x$V2) %>% mutate(cond=paste0(V1, "_", V2))
condition.2 <- x$cond
mat3 <- mat2[, condition.2]


#----plot----
#breaks and color scale
y <- as.numeric(mat2)
quantile(abs(y),probs=0.99)
y0 <- y[abs(y)< quantile(abs(y),probs=0.99)] #99% percent quantile(abs(y),probs=0.99)
mybreaks <- c(min(y),quantile(y0,probs=seq(0,1,length.out=98)),max(y))
names(mybreaks) <- NULL
mycol <- colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100)
col1 <- c("caffeine"="red", "nicotine"="tan",
                    "vitA"="tan4", "vitD"="seagreen4",
                    "vitE"="salmon3", "water"="grey", "zinc"="maroon3")
col2 <- c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
                    "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00", "1-TCM"="pink",
                    "3-TEM"="blue", "5-CD8Naive"="green", "7-dnT"="black")
tmp_colors <- list(treatment=col1, celltype=col2) #brewer.pal(4,"Set1")

##
x <- str_split(condition, "_", simplify=T)
tmp_column <- data.frame(treatment=x[,2], celltype=x[,1])
rownames(tmp_column) <- condition
fig1 <- pheatmap(mat2,  col=mycol, breaks=mybreaks, border_color="NA",
                 cluster_rows=T, cluster_cols=F, treeheight_row = 0,
                 annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend=T,
                 show_colnames=T,
                 show_rownames=T,
                 na_col="white",
                 fontsize_col=7,
                 fontsize_row=2)
figfn <- "./2_motif.activities_2.outs/Figure1.2_heatmap.motif.activities_sig_cell_indfilt_zscore.png"
png(figfn, width=1800, height=2100,res=225)
print(fig1)
dev.off()

##
x <- str_split(condition.2, "_", simplify=T)
tmp_column <- data.frame(treatment=x[,2], celltype=x[,1])
rownames(tmp_column) <- condition.2
fig1 <- pheatmap(mat3,  col=mycol, breaks=mybreaks, border_color="NA",
                 cluster_rows=T, cluster_cols=F,treeheight_row = 0,
                 annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend=T,
                 show_colnames=T,
                 show_rownames=T,
                 na_col="white",
                 fontsize_col=7,
                 fontsize_row=2)
figfn <- "./2_motif.activities_2.outs/Figure1.2_heatmap.motif.activities_sig_treat_indfilt_zscore.png"
png(figfn, width=1800, height=2100,res=225)
print(fig1)
dev.off()


#----sig motifs plotting as done in motif file----
#breaks and color scale
#breaks <- (seq(-8, 8, length.out=17))
#col_fun <- colorRamp2(breaks, colorRampPalette(rev(brewer.pal(n=11, name="RdBu")))(17))
#annotation as cell type
x <- str_split(condition, "_", simplify=T)
x
column_ha_ct <- HeatmapAnnotation(
          treatment=x[,2],
          celltype=x[,1],
          col=list(treatment=c("caffeine"="red", "nicotine"="tan",
                               "vitA"="tan4", "vitD"="seagreen4",
                               "vitE"="salmon3", "water"="grey",
                               "zinc"="maroon3"),
                   celltype=c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
                              "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00",
                              "1-TCM"="pink", "3-TEM"="blue",
                              "5-CD8Naive"="green", "7-dnT"="black")))

fig <- Heatmap(mat2,
               ##col=col_fun,
               cluster_rows=T, cluster_columns=F, show_row_dend = FALSE,
               top_annotation=column_ha_ct,
               heatmap_legend_param=list(title="fold.enrichment",
                                         title_gp=gpar(fontsize=10),
                                         labels_gp=gpar(fontsize=10)),
               show_row_names=F,
               show_column_names=T,
               row_names_gp=gpar(fontsize=2),
               column_names_gp=gpar(fontsize=7),
               raster_device="png")
figfn <- "./2_motif.activities_2.outs/Figure1.3_heatmap.motif.activities_sig_cell_indfilt_zscore.png"
png(figfn, width=1800, height=2100,res=225)
print(fig)
dev.off()

#annotation as treatment
x <- str_split(comb.2, "_", simplify=T)
x
column_ha_treat <- HeatmapAnnotation(
          treatment=x[,2],
          celltype=x[,1],
          col=list(treatment=c("caffeine"="red", "nicotine"="tan",
                               "vitA"="tan4", "vitD"="seagreen4",
                               "vitE"="salmon3", "water"="grey",
                               "zinc"="maroon3"),
                   celltype=c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
                              "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00",
                              "1-TCM"="pink", "3-TEM"="blue",
                              "5-CD8Naive"="green", "7-dnT"="black")))
fig <- Heatmap(mat3,
               #col=col_fun,
               cluster_rows=T, cluster_columns=F, show_row_dend = FALSE,
               top_annotation=column_ha_treat,
               heatmap_legend_param=list(title="fold.enrichment",
                                         title_gp=gpar(fontsize=10),
                                         labels_gp=gpar(fontsize=10)),
               show_row_names=F,
               show_column_names=T,
               row_names_gp=gpar(fontsize=2),
               column_names_gp=gpar(fontsize=7),
               raster_device="png")
figfn <- "./2_motif.activities_2.outs/Figure1.3_heatmap.motif.activities_sig_treat_indfilt_zscore.png"
png(figfn, width=1800, height=2100,res=225)
print(fig)
dev.off()

end





#----------------------------
#----mat.10topmotifs.qval----
#----------------------------
mat.10topmotifs.qval <- map_dfr(condition, function(ii){
       res2 <- res.filt %>% dplyr::filter(condition==ii)
          res2 <- res2 %>% arrange(qval)
          res2 <- res2 [1:10, ]
       })
topmotif.qval <- unique(mat.10topmotifs.qval$motif)
mat.10topmotifs.qval.2 <- map_dfc(condition, function(ii){
      res2 <- res%>%dplyr::filter(condition==ii)
        b <- res2$zscore
        names(b) <- res2$motif
        b[topmotif.qval]
          })
print("unique top 10 motifs")
nrow(mat.10topmotifs.qval.2)

mat <- mat.10topmotifs.qval.2
mat <- as.matrix(mat)
mat[is.na(mat)] = 0
colnames(mat) <- condition
rownames(mat) <- topmotif.qval
head(mat)


mat2 <- mat
x <- as.data.frame(str_split(condition, "_", simplify=T))
x <- x %>% dplyr::arrange(x$V2) %>% mutate(cond=paste0(V1, "_", V2))
condition.2 <- x$cond
mat3 <- mat[, condition.2]

##----plot----
##breaks and color scale
y <- as.numeric(mat2)
quantile(abs(y),probs=0.99)
y0 <- y[abs(y)< quantile(abs(y),probs=0.99)] #99% percent quantile(abs(y),probs=0.99)
mybreaks <- c(min(y),quantile(y0,probs=seq(0,1,length.out=98)),max(y))
names(mybreaks) <- NULL
mycol <- colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100)
col1 <- c("caffeine"="red", "nicotine"="tan",
                    "vitA"="tan4", "vitD"="seagreen4",
                    "vitE"="salmon3", "water"="grey", "zinc"="maroon3")
col2 <- c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
                    "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00", "1-TCM"="pink",
                    "3-TEM"="blue", "5-CD8Naive"="green", "7-dnT"="black")
tmp_colors <- list(treatment=col1, celltype=col2) #brewer.pal(4,"Set1")

##
x <- str_split(condition, "_", simplify=T)
tmp_column <- data.frame(treatment=x[,2], celltype=x[,1])
rownames(tmp_column) <- condition
fig1 <- pheatmap(mat2,  col=mycol, breaks=mybreaks, border_color="NA",
                 cluster_rows=T, cluster_cols=F, treeheight_row = 0,
                 annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend=T,
                 show_colnames=T,
                 show_rownames=T,
                 na_col="white",
                 fontsize_col=7,
                 fontsize_row=4)
figfn <- "./2_motif.activities_2.outs/Figure1.2_heatmap.motif.activities_top10_cell_indfilt_zscore.png"
png(figfn, width=1800, height=2100,res=225)
print(fig1)
dev.off()

##
x <- str_split(condition.2, "_", simplify=T)
tmp_column <- data.frame(treatment=x[,2], celltype=x[,1])
rownames(tmp_column) <- condition.2
fig1 <- pheatmap(mat3,  col=mycol, breaks=mybreaks, border_color="NA",
                 cluster_rows=T, cluster_cols=F,treeheight_row = 0,
                 annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend=T,
                 show_colnames=T,
                 show_rownames=T,
                 na_col="white",
                 fontsize_col=7,
                 fontsize_row=4)
figfn <- "./2_motif.activities_2.outs/Figure1.2_heatmap.motif.activities_top10_treat_indfilt_zscore.png"
png(figfn, width=1800, height=2100,res=225)
print(fig1)
dev.off()


##----top10 motifs plotting as done in motif file----
#breaks and color scale
#breaks <- (seq(-8, 8, length.out=17))
#col_fun <- colorRamp2(breaks, colorRampPalette(rev(brewer.pal(n=11, name="RdBu")))(17))
#annotation as cell type
x <- str_split(condition, "_", simplify=T)
x
column_ha_ct <- HeatmapAnnotation(
      treatment=x[,2],
      celltype=x[,1],
      col=list(treatment=c("caffeine"="red", "nicotine"="tan",
                           "vitA"="tan4", "vitD"="seagreen4",
                           "vitE"="salmon3", "water"="grey", "zinc"="maroon3"),
               celltype=c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
                          "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00",
                          "1-TCM"="pink", "3-TEM"="blue",
                          "5-CD8Naive"="green", "7-dnT"="black")))
fig <- Heatmap(mat2,
               #col=col_fun,
               cluster_rows=T, cluster_columns=F, show_row_dend = FALSE,
               top_annotation=column_ha_ct,
               heatmap_legend_param=list(title="fold.enrichment",
                                         title_gp=gpar(fontsize=10),
                                         labels_gp=gpar(fontsize=10)),
               show_row_names=T,
               show_column_names=T,
               row_names_gp=gpar(fontsize=4),
               column_names_gp=gpar(fontsize=7),
               raster_device="png")
figfn <- "./2_motif.activities_2.outs/Figure1.3_heatmap.motif.activities_top10_cell_indfilt_zscore.png"
png(figfn, width=1800, height=2100,res=225)
print(fig)
dev.off()

#annotation as treatment
x <- str_split(condition.2, "_", simplify=T)
x
column_ha_treat <- HeatmapAnnotation(
      treatment=x[,2],
      celltype=x[,1],
      col=list(treatment=c("caffeine"="red", "nicotine"="tan",
                           "vitA"="tan4", "vitD"="seagreen4",
                           "vitE"="salmon3", "water"="grey", "zinc"="maroon3"),
               celltype=c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
                          "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00",
                          "1-TCM"="pink", "3-TEM"="blue",
                          "5-CD8Naive"="green", "7-dnT"="black")))
fig <- Heatmap(mat3,
               #col=col_fun,
               cluster_rows=T, cluster_columns=F, show_row_dend = FALSE,
               top_annotation=column_ha_treat,
               heatmap_legend_param=list(title="fold.enrichment",
                                         title_gp=gpar(fontsize=10),
                                         labels_gp=gpar(fontsize=10)),
               show_row_names=T,
               show_column_names=T,
               row_names_gp=gpar(fontsize=4),
               column_names_gp=gpar(fontsize=7),
               raster_device="png")
figfn <- "./2_motif.activities_2.outs/Figure1.3_heatmap.motif.activities_top10_treat_indfilt_zscore.png"
png(figfn, width=1800, height=2100,res=225)
print(fig)
dev.off()

end


















##################################################
### Differential activities for ind + sampleID ###
##################################################
res <- read_rds("./2_motif.activities.outs/3_motif.diff.results_macs2_0.1_cn_indfilt_incsampleID_all_cols_control.rds")

res <- res %>% mutate(condition=paste0(gsub("_", "-", MCls), "_", contrast))
#res <- res%>%mutate(condition=paste(MCls, contrast, sep="_"))

res <- res %>% mutate(zscore=beta/stderr)

table(res$MCls, res$contrast)

head(res)

res.filt <- res%>%dplyr::filter(qval<0.1, abs(beta)>0)
topmotif <- sort(unique(res.filt$motif))
print("Unique significant motifs are:")
length(topmotif)


sig_motifs_table <- table(res.filt$MCls, res.filt$contrast)
sig_motifs_table

sig.matrix <-as.matrix(sig_motifs_table) # create a numeric matrix object
sig.melt <- melt(sig.matrix)
sig.p <- ggplot(sig.melt, aes(x = Var2, y = Var1)) +
      geom_tile(aes(fill=value), color="black")+
      scale_fill_gradient(low = "white", high = "white")+
      geom_text(aes(label = value), color = "black", size = 5) +
      ggtitle(paste0("Significant_motifs_coef_0")) +
      theme(axis.text.x = element_text(angle = 90),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    legend.position="none")
figure <- ggarrange(sig.p + font("x.text", size = 14))
png(paste0("./2_motif.activities_2.outs/sig_motifs_table_indfilt_incsampleID_coef_0_control.png"), width=500, height=500, pointsize=16, res=125)
print(figure)
dev.off()




##--------------------------------
##----mat data for significant----
##--------------------------------
head(res)
table(res$MCls)
table(res$contrast)

condition <- unique(sort(res$condition))
condition

mat.sigmotifs.qval <- map_dfr(condition, function(ii){
    res2 <- res.filt %>% dplyr::filter(condition==ii)
})

##
sigmotif.qval <- unique(mat.sigmotifs.qval$motif)

##
mat.sigmotifs.qval.2 <- map_dfc(condition, function(ii){
    res2 <- res%>%dplyr::filter(condition==ii)
    b <- res2$zscore
    names(b) <- res2$motif
    b[sigmotif.qval]
})

##
mat <- mat.sigmotifs.qval.2
mat <- as.matrix(mat)
mat[is.na(mat)] = 0
colnames(mat) <- condition
rownames(mat) <- sigmotif.qval
head(mat)
print("unique sig motifs")
nrow(mat.sigmotifs.qval.2)

mat2 <- mat

x <- as.data.frame(str_split(condition, "_", simplify=T))
x <- x %>% dplyr::arrange(x$V2) %>% mutate(cond=paste0(V1, "_", V2))
condition.2 <- x$cond
mat3 <- mat2[, condition.2]

##----plot----
#breaks and color scale
## y <- as.numeric(mat2)
## quantile(abs(y),probs=0.99)
## y0 <- y[abs(y)< quantile(abs(y),probs=0.99)] #99% percent quantile(abs(y),probs=0.99)
## mybreaks <- c(min(y),quantile(y0,probs=seq(0,1,length.out=98)),max(y))
## names(mybreaks) <- NULL

y <- as.numeric(mat2)
max(y)
min(y)
scale <- quantile(abs(y),probs=0.99)
scale
mybreaks <- seq(-scale, scale, length.out=100)
names(mybreaks) <- NULL
mybreaks



##
mycol <- colorRampPalette(rev(brewer.pal(n=11, name="RdBu")))(100)
col1 <- c("caffeine"="red", "nicotine"="tan",
          "vitA"="tan4", "vitD"="seagreen4",
          "vitE"="salmon3", "water"="grey", "zinc"="maroon3")
col2 <- c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
          "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00", "1-TCM"="pink",
          "3-TEM"="blue", "5-CD8Naive"="green", "7-dnT"="black")
tmp_colors <- list(treatment=col1, celltype=col2) #brewer.pal(4,"Set1")

##
x <- str_split(condition, "_", simplify=T)
tmp_column <- data.frame(treatment=x[,2], celltype=x[,1])
rownames(tmp_column) <- condition
fig1 <- pheatmap(mat2,  col=mycol, breaks=mybreaks, border_color="NA",
                 cluster_rows=T, cluster_cols=F, treeheight_row = 0,
                 annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend=T,
                 show_colnames=T,
                 show_rownames=T,
                 na_col="white",
                 fontsize_col=7,
                 fontsize_row=2)

##
figfn <- "./2_motif.activities_2.outs/Figure1.2_heatmap.motif.activities_sig_cell_indfilt_incsampleID_zscore_control.png"
png(figfn, width=1800, height=2100,res=225)
print(fig1)
dev.off()

##
x <- str_split(condition.2, "_", simplify=T)
tmp_column <- data.frame(treatment=x[,2], celltype=x[,1])
rownames(tmp_column) <- condition.2
fig1 <- pheatmap(mat3,  col=mycol, breaks=mybreaks, border_color="NA",
                 cluster_rows=T, cluster_cols=F,treeheight_row = 0,
                 annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend=T,
                 show_colnames=T,
                 show_rownames=T,
                 na_col="white",
                 fontsize_col=7,
                 fontsize_row=2)

##
figfn <- "./2_motif.activities_2.outs/Figure1.2_heatmap.motif.activities_sig_treat_indfilt_incsampleID_zscore_control.png"
png(figfn, width=1800, height=2100,res=225)
print(fig1)
dev.off()


## #----sig motifs plotting as done in motif file----
## #breaks and color scale
## #breaks <- (seq(-8, 8, length.out=17))
## #col_fun <- colorRamp2(breaks, colorRampPalette(rev(brewer.pal(n=11, name="RdBu")))(17))
## #annotation as cell type
## x <- str_split(condition, "_", simplify=T)
## x
## column_ha_ct <- HeatmapAnnotation(
##           treatment=x[,2],
##           celltype=x[,1],
##           col=list(treatment=c("caffeine"="red", "nicotine"="tan",
##                                "vitA"="tan4", "vitD"="seagreen4",
##                                "vitE"="salmon3", "water"="grey",
##                                "zinc"="maroon3"),
##                    celltype=c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
##                               "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00",
##                               "1-TCM"="pink", "3-TEM"="blue",
##                               "5-CD8Naive"="green", "7-dnT"="black")))

## fig <- Heatmap(mat2,
##                ##col=col_fun,
##                cluster_rows=T, cluster_columns=F, show_row_dend = FALSE,
##                top_annotation=column_ha_ct,
##                heatmap_legend_param=list(title="fold.enrichment",
##                                          title_gp=gpar(fontsize=10),
##                                          labels_gp=gpar(fontsize=10)),
##                show_row_names=F,
##                show_column_names=T,
##                row_names_gp=gpar(fontsize=2),
##                column_names_gp=gpar(fontsize=7),
##                raster_device="png")
## figfn <- "./2_motif.activities_2.outs/Figure1.3_heatmap.motif.activities_sig_cell_indfilt_incsampleID.png"
## png(figfn, width=1800, height=2100,res=225)
## print(fig)
## dev.off()

## #annotation as treatment
## x <- str_split(comb.2, "_", simplify=T)
## x
## column_ha_treat <- HeatmapAnnotation(
##           treatment=x[,2],
##           celltype=x[,1],
##           col=list(treatment=c("caffeine"="red", "nicotine"="tan",
##                                "vitA"="tan4", "vitD"="seagreen4",
##                                "vitE"="salmon3", "water"="grey",
##                                "zinc"="maroon3"),
##                    celltype=c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
##                               "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00",
##                               "1-TCM"="pink", "3-TEM"="blue",
##                               "5-CD8Naive"="green", "7-dnT"="black")))

## fig <- Heatmap(mat3,
##                #col=col_fun,
##                cluster_rows=T, cluster_columns=F, show_row_dend = FALSE,
##                top_annotation=column_ha_treat,
##                heatmap_legend_param=list(title="fold.enrichment",
##                                          title_gp=gpar(fontsize=10),
##                                          labels_gp=gpar(fontsize=10)),
##                show_row_names=F,
##                show_column_names=T,
##                row_names_gp=gpar(fontsize=2),
##                column_names_gp=gpar(fontsize=7),
##                raster_device="png")
## figfn <- "./2_motif.activities_2.outs/Figure1.3_heatmap.motif.activities_sig_treat_indfilt_incsampleID.png"
## png(figfn, width=1800, height=2100,res=225)
## print(fig)
## dev.off()

## end
## #--------













##----------------------------
##----mat.10topmotifs.qval----
##----------------------------
head(res.filt)

mat.10topmotifs.qval <- map_dfr(condition, function(ii){
       res2 <- res.filt %>% dplyr::filter(condition==ii)
          res2 <- res2 %>% arrange(qval)
          res2 <- res2 [1:10, ]
       })

topmotif.qval <- unique(mat.10topmotifs.qval$motif)

mat.10topmotifs.qval.2 <- map_dfc(condition, function(ii){
    res2 <- res%>%dplyr::filter(condition==ii)
    b <- res2$beta
    names(b) <- res2$motif
    b[topmotif.qval]
})
print("unique top 10 motifs")
nrow(mat.10topmotifs.qval.2)

mat <- mat.10topmotifs.qval.2
mat <- as.matrix(mat)
mat[is.na(mat)] = 0
colnames(mat) <- condition
rownames(mat) <- topmotif.qval
head(mat)


mat2 <- mat
x <- as.data.frame(str_split(condition, "_", simplify=T))
x <- x %>% dplyr::arrange(x$V2) %>% mutate(cond=paste0(V1, "_", V2))
condition.2 <- x$cond
mat3 <- mat[, condition.2]

##----plot----
##breaks and color scale
## y <- as.numeric(mat2)
## quantile(abs(y),probs=0.99)
## y0 <- y[abs(y)< quantile(abs(y),probs=0.99)] #99% percent quantile(abs(y),probs=0.99)
## mybreaks <- c(min(y),quantile(y0,probs=seq(0,1,length.out=98)),max(y))
## names(mybreaks) <- NULL

y <- as.numeric(mat2)
max(y)
min(y)
scale <- quantile(abs(y),probs=0.99)
scale
mybreaks <- seq(-scale, scale, length.out=100)
names(mybreaks) <- NULL
mybreaks


##
mycol <- colorRampPalette(rev(brewer.pal(n=11, name="RdBu")))(100)
col1 <- c("caffeine"="red", "nicotine"="tan",
          "vitA"="tan4", "vitD"="seagreen4",
          "vitE"="salmon3", "water"="grey", "zinc"="maroon3")
col2 <- c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
          "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00", "1-TCM"="pink",
          "3-TEM"="blue", "5-CD8Naive"="green", "7-dnT"="black")
tmp_colors <- list(treatment=col1, celltype=col2) #brewer.pal(4,"Set1")

##
x <- str_split(condition, "_", simplify=T)
tmp_column <- data.frame(treatment=x[,2], celltype=x[,1])
rownames(tmp_column) <- condition
fig1 <- pheatmap(mat2,  col=mycol, breaks=mybreaks, border_color="NA",
                 cluster_rows=T, cluster_cols=F, treeheight_row = 0,
                 annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend=T,
                 show_colnames=T,
                 show_rownames=T,
                 na_col="white",
                 fontsize_col=7,
                 fontsize_row=4)

figfn <- "./2_motif.activities_2.outs/Figure1.2_heatmap.motif.activities_top10_cell_indfilt_incsampleID_beta_control.png"
png(figfn, width=1800, height=2100,res=225)
print(fig1)
dev.off()

##
x <- str_split(condition.2, "_", simplify=T)
tmp_column <- data.frame(treatment=x[,2], celltype=x[,1])
rownames(tmp_column) <- condition.2
fig1 <- pheatmap(mat3,  col=mycol, breaks=mybreaks, border_color="NA",
                 cluster_rows=T, cluster_cols=F,treeheight_row = 0,
                 annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend=T,
                 show_colnames=T,
                 show_rownames=T,
                 na_col="white",
                 fontsize_col=7,
                 fontsize_row=4)

figfn <- "./2_motif.activities_2.outs/Figure1.2_heatmap.motif.activities_top10_treat_indfilt_incsampleID_beta_control.png"
png(figfn, width=1800, height=2100,res=225)
print(fig1)
dev.off()



## ##----top10 motifs plotting as done in motif file----
## #breaks and color scale
## #breaks <- (seq(-8, 8, length.out=17))
## #col_fun <- colorRamp2(breaks, colorRampPalette(rev(brewer.pal(n=11, name="RdBu")))(17))
## #annotation as cell type
## x <- str_split(condition, "_", simplify=T)
## x
## column_ha_ct <- HeatmapAnnotation(
##       treatment=x[,2],
##       celltype=x[,1],
##       col=list(treatment=c("caffeine"="red", "nicotine"="tan",
##                            "vitA"="tan4", "vitD"="seagreen4",
##                            "vitE"="salmon3", "water"="grey", "zinc"="maroon3"),
##                celltype=c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
##                           "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00",
##                           "1-TCM"="pink", "3-TEM"="blue",
##                           "5-CD8Naive"="green", "7-dnT"="black")))
## fig <- Heatmap(mat2,
##                #col=col_fun,
##                cluster_rows=T, cluster_columns=F, show_row_dend = FALSE,
##                top_annotation=column_ha_ct,
##                heatmap_legend_param=list(title="fold.enrichment",
##                                          title_gp=gpar(fontsize=10),
##                                          labels_gp=gpar(fontsize=10)),
##                show_row_names=T,
##                show_column_names=T,
##                row_names_gp=gpar(fontsize=4),
##                column_names_gp=gpar(fontsize=7),
##                raster_device="png")
## figfn <- "./2_motif.activities_2.outs/Figure1.3_heatmap.motif.activities_top10_cell_indfilt_incsampleID_zscore.png"
## png(figfn, width=1800, height=2100,res=225)
## print(fig)
## dev.off()

## #annotation as treatment
## x <- str_split(condition.2, "_", simplify=T)
## x
## column_ha_treat <- HeatmapAnnotation(
##       treatment=x[,2],
##       celltype=x[,1],
##       col=list(treatment=c("caffeine"="red", "nicotine"="tan",
##                            "vitA"="tan4", "vitD"="seagreen4",
##                            "vitE"="salmon3", "water"="grey", "zinc"="maroon3"),
##                celltype=c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
##                           "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00",
##                           "1-TCM"="pink", "3-TEM"="blue",
##                           "5-CD8Naive"="green", "7-dnT"="black")))
## fig <- Heatmap(mat3,
##                #col=col_fun,
##                cluster_rows=T, cluster_columns=F, show_row_dend = FALSE,
##                top_annotation=column_ha_treat,
##                heatmap_legend_param=list(title="fold.enrichment",
##                                          title_gp=gpar(fontsize=10),
##                                          labels_gp=gpar(fontsize=10)),
##                show_row_names=T,
##                show_column_names=T,
##                row_names_gp=gpar(fontsize=4),
##                column_names_gp=gpar(fontsize=7),
##                raster_device="png")
## figfn <- "./2_motif.activities_2.outs/Figure1.3_heatmap.motif.activities_top10_treat_indfilt_incsampleID_zscore.png"
## png(figfn, width=1800, height=2100,res=225)
## print(fig)
## dev.off()

## end
## #--------

## head(res)







