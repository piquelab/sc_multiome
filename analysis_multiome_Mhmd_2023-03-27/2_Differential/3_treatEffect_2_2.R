##
library(Matrix)
library(tidyverse)
library(DESeq2)
library(annotables)
##
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(ggExtra)
library(gtable)
library(ggsignif)
library(pheatmap)
library(ComplexHeatmap)
library(corrplot)
library(gtable)
library(RColorBrewer)
library(viridis)
library(ggrastr)
library(ggsci)
library(circlize)
library(reshape2)
library("ggpubr")
# library(ggsci)
# library(circlize)

rm(list=ls())

##
outdir <- "./3_treatEffect_2_2.outs/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=F)


#####################################################################
###            heatmap for resDP                                  ###
#####################################################################
#----resDP----
##resDP <- read_rds("../2_Differential/1.2_DiffPeak.outs/2.0_DESeq.results_0.1_cn.rds") %>% as.data.frame()

resDP <- read_rds("../2_Differential/1.2_DiffPeak.outs/2.0_DESeq.results_0.1_cn_sampleID+new_treat.rds") %>% as.data.frame()
nrow(resDP)
nrow(resDP %>% dplyr::filter(is.na(p.adjusted)))

#resDP <- resDP %>% mutate(MCls=gsub("_", "-", MCls)) %>% mutate(comb=paste(MCls, contrast, sep="_"))%>% dplyr::rename("peak"="gene")
#%>% drop_na(p.value)
resDP <- resDP %>% mutate(MCls=gsub("_", "-", MCls)) %>% mutate(comb=paste(contrast, MCls, sep="_"))%>% dplyr::rename("peak"="gene")
#%>% drop_na(p.value)
head(resDP)
nrow(resDP)
max(resDP$estimate)
min(resDP$estimate)


resDP[!is.na(resDP$p.adjusted) & resDP$p.adjusted<0.1 & resDP$estimate > 0.5, "Significance"] <- "Up"
resDP[!is.na(resDP$p.adjusted) & resDP$p.adjusted<0.1 & resDP$estimate < -0.5, "Significance"] <- "Down"
resDP[!is.na(resDP$p.adjusted) & resDP$p.adjusted>=0.1, "Significance"] <- "Not Significant"
resDP[is.na(resDP$p.adjusted),"Significance"] <- "NA"
resDP <- resDP %>% mutate(comb_2=paste(comb, Significance, sep="_"))
head(resDP)
#table(resDP$Significance)
#table(resDP$MCls, resDP$contrast)


resDP[!is.na(resDP$p.adjusted) & resDP$p.adjusted<0.1 & resDP$estimate > 0.25, "Significance.25"] <- "Up"
resDP[!is.na(resDP$p.adjusted) & resDP$p.adjusted<0.1 & resDP$estimate < -0.25, "Significance.25"] <- "Down"
resDP[!is.na(resDP$p.adjusted) & resDP$p.adjusted>=0.1, "Significance.25"] <- "Not Significant"
resDP[is.na(resDP$p.adjusted),"Significance.25"] <- "NA"
#head(resDP)
table(resDP$Significance.25)


resDP[!is.na(resDP$p.adjusted) & resDP$p.adjusted<0.1 & resDP$estimate > 0, "Significance.0"] <- "Up"
resDP[!is.na(resDP$p.adjusted) & resDP$p.adjusted<0.1 & resDP$estimate < 0, "Significance.0"] <- "Down"
resDP[!is.na(resDP$p.adjusted) & resDP$p.adjusted>=0.1, "Significance.0"] <- "Not Significant"
resDP[is.na(resDP$p.adjusted),"Significance.0"] <- "NA"
#head(resDP)
table(resDP$Significance.0)

opfn <- paste0("./3_treatEffect_2_2.outs/resDP_control.rds")
write_rds(resDP, opfn)

resDP <- read_rds(opfn)


#----allDP----
all.DP  <- resDP %>% dplyr::pull(peak)
all.DP <- as.character(unique(all.DP))
head(all.DP)
length(all.DP)


#----sig_resDP----
sig_resDP.5 <- resDP %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)
sig_resDP.25 <- resDP %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<0.1, abs(estimate)>0.25)
sig_resDP.0 <- resDP %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<0.1, abs(estimate)>0)


sig <- table(sig_resDP.5$MCls, sig_resDP.5$contrast)
#sig
sig.matrix <- as.matrix(sig)
sig.melt <- melt(sig.matrix)
sig.p <- ggplot(sig.melt, aes(x = Var2, y = Var1)) +
      geom_tile(aes(fill=value), color="black")+
      scale_fill_gradient(low = "white", high = "white")+
      geom_text(aes(label = value), color = "black", size = 5) +
      ggtitle("Significant_0.1_0.1-0.5") +
      theme(axis.text.x = element_text(angle = 90),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            legend.position="none")
figure <- ggarrange(sig.p +
                    font("x.text", size = 14))
# ncol = 2, nrow = 2)
png(paste0("./3_treatEffect_2_2.outs/plot_sig_treatandind_filt_DP_all_0.1_0.1_0.5_cn_mitofilt_control.png"), width=500, height=500, pointsize=16, res=125)
print(figure)
dev.off()


sig <- table(sig_resDP.25$MCls, sig_resDP.25$contrast)
sig.matrix <- as.matrix(sig)
sig.melt <- melt(sig.matrix)
sig.p <- ggplot(sig.melt, aes(x = Var2, y = Var1)) +
      geom_tile(aes(fill=value), color="black")+
      scale_fill_gradient(low = "white", high = "white")+
      geom_text(aes(label = value), color = "black", size = 5) +
      ggtitle("Significant_0.1_0.1-0.25") +
      theme(axis.text.x = element_text(angle = 90),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    legend.position="none")
figure <- ggarrange(sig.p +
                    font("x.text", size = 14))
# ncol = 2, nrow = 2)
png(paste0("./3_treatEffect_2_2.outs/plot_sig_treatandind_filt_DP_all_0.1_0.1_0.25_cn_mitofilt.png"), width=500, height=500, pointsize=16, res=125)
print(figure)
dev.off()




sig <- table(sig_resDP.0$MCls, sig_resDP.0$contrast)
sig.matrix <- as.matrix(sig)
sig.melt <- melt(sig.matrix)
sig.p <- ggplot(sig.melt, aes(x = Var2, y = Var1)) +
      geom_tile(aes(fill=value), color="black")+
      scale_fill_gradient(low = "white", high = "white")+
      geom_text(aes(label = value), color = "black", size = 5) +
      ggtitle("Significant_0.1_0.1-0") +
      theme(axis.text.x = element_text(angle = 90),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    legend.position="none")
figure <- ggarrange(sig.p +
                                          font("x.text", size = 14))
# ncol = 2, nrow = 2)
png(paste0("./3_treatEffect_2_2.outs/plot_sig_treatandind_filt_DP_all_0.1_0.1_0_cn_mitofilt.png"), width=500, height=500, pointsize=16, res=125)
print(figure)
dev.off()


sig_resDP <- sig_resDP.5


#----to check----
nrow(sig_resDP %>% dplyr::filter(abs(estimate)>5))


#----DP-----
DP  <- sig_resDP %>% dplyr::pull(peak)
DP <- as.character(unique(DP))
head(DP)
length(DP)




#######################################
#----mat data for significant----
#######################################
comb <- unique(resDP$comb)
sigdps.qval <- unique(sig_resDP$peak)

##----mat_estimate----
mat.sigdps.qval.2 <- map_dfc(comb, function(ii){
      res2 <- resDP %>% dplyr::filter(comb==ii)
        ##b <- res2$estimate
        b <- res2$statistic
        names(b) <- res2$peak
        b[sigdps.qval]
      })
print("Number of significant in mat")
nrow(mat.sigdps.qval.2)
mat <- mat.sigdps.qval.2


mat <- as.matrix(mat)
rownames(mat) <- sigdps.qval #DP
colnames(mat) <- comb
head(mat)
max(mat)
min(mat)

##
mat[is.na(mat)] = 0
mat2 <- mat
comb.2 <- sort(comb)
mat3 <- mat2[, comb.2]


#----DP sig plot Fig1.1----
#breaks and color scale
y <- as.numeric(mat2)
max(y)
min(y)

scale <- quantile(abs(y),probs=0.99)
scale

##y0 <- y[abs(y)< quantile(abs(y),probs=0.99)] #99% percent quantile(abs(y),probs=0.99)
##mybreaks <- c(min(y),quantile(y0,probs=seq(0,1,length.out=98)),max(y))

mybreaks <- seq(-scale, scale, length.out=100)
mybreaks

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
x <- str_split(comb, "_", simplify=T)
tmp_column <- data.frame(treatment=x[,1], celltype=x[,2])
rownames(tmp_column) <- comb
fig1 <- pheatmap(mat2,  col=mycol, breaks=mybreaks, border_color="NA",
                 cluster_rows=T, cluster_cols=F,
                 annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend=T,
                 show_colnames=T,
                 show_rownames=F,
                 treeheight_row = 0,
                 na_col="white",
                 fontsize_col=7,
                 fontsize_row=3)
figfn <- "./3_treatEffect_2_2.outs/Figure1.1_heatmap.beta.DP.sig.cell_control_quantile_0.99.png"
png(figfn, width=1800, height=2100,res=225)
print(fig1)
dev.off()
##
x <- str_split(comb.2, "_", simplify=T)
tmp_column <- data.frame(treatment=x[,1], celltype=x[,2])
rownames(tmp_column) <- comb.2
fig1 <- pheatmap(mat3,  col=mycol, breaks=mybreaks, border_color="NA",
                 cluster_rows=T, cluster_cols=F,
                 annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend=T,
                 show_colnames=T,
                 show_rownames=F,
                 treeheight_row = 0,
                 na_col="white",
                 fontsize_col=7,
                 fontsize_row=3)
figfn <- "./3_treatEffect_2_2.outs/Figure1.1_heatmap.beta.DP.sig.treat_control_quantile_0.99.png"
png(figfn, width=1800, height=2100,res=225)
print(fig1)
dev.off()



###########################
### correlation heatmap ###
###########################
#Neworder <- c("Monocyte_LPS", "Monocyte_PHA", "Bcell_LPS", "Bcell_PHA",
   "NKcell_LPS", "NKcell_PHA", "Tcell_LPS", "Tcell_PHA",
   "Monocyte_LPS-DEX", "Monocyte_PHA-DEX", "Bcell_LPS-DEX", "Bcell_PHA-DEX",
   "NKcell_LPS-DEX", "NKcell_PHA-DEX", "Tcell_LPS-DEX", "Tcell_PHA-DEX")

head(comb.2)
head(mat3) #mat3 treat is ordered per treatment

Neworder <- c(comb.2)
Neworder

#corr <- cor(mat3)[Neworder, Neworder]
corr <- cor(mat3, method = "spearman")[Neworder, Neworder]
head(corr)
length(corr)

#as.numeric(cor.test(res1$logFC, res1$estimate, method = "spearman")$p.value)

##breaks and color scale
y <- as.numeric(corr)
##min(y)
##max(y)

##quantile(abs(y),probs=0.99)
##quantile
##y0 <- y[abs(y)< quantile(abs(y),probs=0.99)] #99% percent quantile(abs(y),probs=0.99)
##y0
##seq(0,1,length.out=98)
##mybreaks <- c(min(y),quantile(y0,probs=seq(0,1,length.out=98)),max(y))
##mybreaks <- c(-1,quantile(y0,probs=seq(0,1,length.out=98)),max(y))

mybreaks <- seq(-1, 1, length.out = 100)
mybreaks

names(mybreaks) <- NULL

##
## mycol <- colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100)
col1 <- c("caffeine"="red", "nicotine"="tan",
                     "vitA"="tan4", "vitD"="seagreen4",
                     "vitE"="salmon3", "water"="grey", "zinc"="maroon3")
col2 <- c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
                     "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00", "1-TCM"="pink",
                     "3-TEM"="blue", "5-CD8Naive"="green", "7-dnT"="black")
## tmp_colors <- list(treatment=col1, celltype=col2) #brewer.pal(4,"Set1")

mycol <- colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100)


##
x <- str_split(colnames(corr), "_", simplify=T) %>% as.data.frame()
#head(x)
#colnames(x)
#x <- x[order(x$V2),]
x

tmp_column <- data.frame(celltype=x[,2], treatment=x[,1])
rownames(tmp_column) <- colnames(corr)

tmp_colors <- list(celltype=col2, treatment=col1)

fig2 <- pheatmap(corr, col=mycol,
                 breaks=mybreaks,
                 scale="none", border_color="NA",
                 cluster_rows=F, cluster_cols=F,
                 annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend =T,
                 show_colnames=T, show_rownames=F,
                 fontsize_row=5,
                 fontsize_col=5,
                 na_col="white")

###
figfn <- "./3_treatEffect_2_2.outs/Figure1.2_heatmap.corr.DP.sig.treat_control_2_spearman.png"
png(figfn, width=600, height=600,res=120)
print(fig2)
dev.off()


##----as celltype----
comb.cell <- str_split(comb.2, "_", simplify=T) %>% as.data.frame()
comb.cell <- comb.cell[order(comb.cell$V2),] %>% mutate(V3 = paste0(V1, "_", V2))
head(comb.cell)

head(mat2)

Neworder <- c(comb.cell$V3)
Neworder

corr <- cor(mat2, method = "spearman")[Neworder, Neworder]
head(corr)
length(corr)

##breaks and color scale
y <- as.numeric(corr)
##min(y)
##max(y)

##quantile(abs(y),probs=0.99)
##quantile
##y0 <- y[abs(y)< quantile(abs(y),probs=0.99)] #99% percent quantile(abs(y),probs=0.99)
##y0
##seq(0,1,length.out=98)
##mybreaks <- c(min(y),quantile(y0,probs=seq(0,1,length.out=98)),max(y))
##mybreaks <- c(-1,quantile(y0,probs=seq(0,1,length.out=98)),max(y))

mybreaks <- seq(-1, 1, length.out = 100)
mybreaks

names(mybreaks) <- NULL

##
## mycol <- colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100)
## col1 <- c("caffeine"="red", "nicotine"="tan",
##                     "vitA"="tan4", "vitD"="seagreen4",
##                     "vitE"="salmon3", "water"="grey", "zinc"="maroon3")
## col2 <- c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
##                     "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00", "1-TCM"="pink",
##                     "3-TEM"="blue", "5-CD8Naive"="green", "7-dnT"="black")
## tmp_colors <- list(treatment=col1, celltype=col2) #brewer.pal(4,"Set1")
mycol <- colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100)


##
x <- str_split(colnames(corr), "_", simplify=T) %>% as.data.frame()
#head(x)
#colnames(x)
#x <- x[order(x$V2),]
x

tmp_column <- data.frame(celltype=x[,2], treatment=x[,1])
rownames(tmp_column) <- colnames(corr)

tmp_colors <- list(celltype=col2, treatment=col1)

fig2 <- pheatmap(corr, col=mycol,
                 breaks=mybreaks,
                 scale="none", border_color="NA",
                 cluster_rows=F, cluster_cols=F,
                 annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend =T,
                 show_colnames=T, show_rownames=F,
                 fontsize_row=5,
                 fontsize_col=5,
                 na_col="white")

###
figfn <- "./3_treatEffect_2_2.outs/Figure1.2_heatmap.corr.DP.sig.cell_control_2_spearman.png"
png(figfn, width=600, height=600,res=120)
print(fig2)
dev.off()





















##----mat_zscore----
mat.sigdps.qval.2 <- map_dfc(comb, function(ii){
      res2 <- resDP %>% dplyr::filter(comb==ii)
        b <- res2$statistic
        names(b) <- res2$peak
        b[sigdps.qval]
      })
print("Number of significant in mat")
nrow(mat.sigdps.qval.2)
mat <- mat.sigdps.qval.2


mat <- as.matrix(mat)
rownames(mat) <- sigdps.qval #DP
colnames(mat) <- comb
head(mat)
max(mat)
min(mat)

##
mat[is.na(mat)] = 0
mat2 <- mat
comb.2 <- sort(comb)
mat3 <- mat2[, comb.2]



#----DP sig plot Fig1.1----
#breaks and color scale
y <- as.numeric(mat2)
max(y)
min(y)

scale <- quantile(abs(y),probs=0.99)
scale

##y0 <- y[abs(y)< quantile(abs(y),probs=0.99)] #99% percent quantile(abs(y),probs=0.99)
##mybreaks <- c(min(y),quantile(y0,probs=seq(0,1,length.out=98)),max(y))

mybreaks <- seq(-scale, scale, length.out=100)
mybreaks
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
x <- str_split(comb, "_", simplify=T)
tmp_column <- data.frame(treatment=x[,1], celltype=x[,2])
rownames(tmp_column) <- comb
fig1 <- pheatmap(mat2,  col=mycol, breaks=mybreaks, border_color="NA",
                 cluster_rows=T, cluster_cols=F,
                 annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend=T,
                 show_colnames=T,
                 show_rownames=F,
                 treeheight_row = 0,
                 na_col="white",
                 fontsize_col=7,
                 fontsize_row=3)
figfn <- "./3_treatEffect_2_2.outs/Figure1.1_heatmap.beta.DP.sig.cell_zscore_control_quantile_0.99.png"
png(figfn, width=1800, height=2100,res=225)
print(fig1)
dev.off()
##
x <- str_split(comb.2, "_", simplify=T)
tmp_column <- data.frame(treatment=x[,1], celltype=x[,2])
rownames(tmp_column) <- comb.2
fig1 <- pheatmap(mat3,  col=mycol, breaks=mybreaks, border_color="NA",
                 cluster_rows=T, cluster_cols=F,
                 annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend=T,
                 show_colnames=T,
                 show_rownames=F,
                 treeheight_row = 0,
                 na_col="white",
                 fontsize_col=7,
                 fontsize_row=3)
figfn <- "./3_treatEffect_2_2.outs/Figure1.1_heatmap.beta.DP.sig.treat_zscore_control_quantile_0.99.png"
png(figfn, width=1800, height=2100,res=225)
print(fig1)
dev.off()


## #----sig DP plotting as done in motif file----
## #breaks and color scale
## #breaks <- (seq(-8, 8, length.out=17))
## #col_fun <- colorRamp2(breaks, colorRampPalette(rev(brewer.pal(n=11, name="RdBu")))(17))
## #annotation as cell type
## x <- str_split(comb, "_", simplify=T)
## x
## column_ha_ct <- HeatmapAnnotation(
##       treatment=x[,1],
##       celltype=x[,2],
##       col=list(treatment=c("caffeine"="red", "nicotine"="tan",
##                            "vitA"="tan4", "vitD"="seagreen4",
##                            "vitE"="salmon3", "water"="grey", "zinc"="maroon3"),
##                celltype=c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
##                           "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00",
##                           "1-TCM"="pink", "3-TEM"="blue",
##                           "5-CD8Naive"="green", "7-dnT"="black")))
## fig <- Heatmap(mat2,
## #               col=col_fun,
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
## figfn <- "./3_treatEffect_2_2.outs/Figure1.1.2_heatmap.beta.DP.sig.cell_zscore.png"
## png(figfn, width=1800, height=2100,res=225)
## print(fig)
## dev.off()
## #annotation as treatment
## x <- str_split(comb.2, "_", simplify=T)
## x
## column_ha_treat <- HeatmapAnnotation(
##       treatment=x[,1],
##       celltype=x[,2],
##       col=list(treatment=c("caffeine"="red", "nicotine"="tan",
##                            "vitA"="tan4", "vitD"="seagreen4",
##                            "vitE"="salmon3", "water"="grey", "zinc"="maroon3"),
##                celltype=c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
##                           "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00",
##                           "1-TCM"="pink", "3-TEM"="blue",
##                           "5-CD8Naive"="green", "7-dnT"="black")))
## fig <- Heatmap(mat3,
## #               col=col_fun,
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
## figfn <- "./3_treatEffect_2_2.outs/Figure1.1.2_heatmap.beta.DP.sig.treat_zscore.png"
## png(figfn, width=1800, height=2100,res=225)
## print(fig)
## dev.off()

## end
## #--------



#####################################
### correlation heatmap - zscore  ###
#####################################
#Neworder <- c("Monocyte_LPS", "Monocyte_PHA", "Bcell_LPS", "Bcell_PHA",
   "NKcell_LPS", "NKcell_PHA", "Tcell_LPS", "Tcell_PHA",
   "Monocyte_LPS-DEX", "Monocyte_PHA-DEX", "Bcell_LPS-DEX", "Bcell_PHA-DEX",
   "NKcell_LPS-DEX", "NKcell_PHA-DEX", "Tcell_LPS-DEX", "Tcell_PHA-DEX")

head(comb.2)
head(mat3) #mat3 treat is ordered per treatment

Neworder <- c(comb.2)
Neworder

corr <- cor(mat3)[Neworder, Neworder]
head(corr)
length(corr)

##breaks and color scale
y <- as.numeric(corr)
##min(y)
##max(y)

##quantile(abs(y),probs=0.99)
##quantile
##y0 <- y[abs(y)< quantile(abs(y),probs=0.99)] #99% percent quantile(abs(y),probs=0.99)
##y0
##seq(0,1,length.out=98)
##mybreaks <- c(min(y),quantile(y0,probs=seq(0,1,length.out=98)),max(y))
##mybreaks <- c(-1,quantile(y0,probs=seq(0,1,length.out=98)),max(y))

mybreaks <- seq(-1, 1, length.out = 100)
mybreaks

names(mybreaks) <- NULL

##
## mycol <- colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100)
## col1 <- c("caffeine"="red", "nicotine"="tan",
##                     "vitA"="tan4", "vitD"="seagreen4",
##                     "vitE"="salmon3", "water"="grey", "zinc"="maroon3")
## col2 <- c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
##                     "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00", "1-TCM"="pink",
##                     "3-TEM"="blue", "5-CD8Naive"="green", "7-dnT"="black")
## tmp_colors <- list(treatment=col1, celltype=col2) #brewer.pal(4,"Set1")
mycol <- colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100)


##
x <- str_split(colnames(corr), "_", simplify=T) %>% as.data.frame()
#head(x)
#colnames(x)
#x <- x[order(x$V2),]
x

tmp_column <- data.frame(celltype=x[,2], treatment=x[,1])
rownames(tmp_column) <- colnames(corr)

tmp_colors <- list(celltype=col2, treatment=col1)

fig2 <- pheatmap(corr, col=mycol,
                 breaks=mybreaks,
                 scale="none", border_color="NA",
                 cluster_rows=F, cluster_cols=F,
                 annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend =T,
                 show_colnames=T, show_rownames=F,
                 fontsize_row=5,
                 fontsize_col=5,
                 na_col="white")

###
figfn <- "./3_treatEffect_2_2.outs/Figure1.2_heatmap.corr.DP.sig.treat_zscore_control_2.png"
png(figfn, width=600, height=600,res=120)
print(fig2)
dev.off()


##----as celltype----
comb.cell <- str_split(comb.2, "_", simplify=T) %>% as.data.frame()
comb.cell <- comb.cell[order(comb.cell$V2),] %>% mutate(V3 = paste0(V1, "_", V2))
head(comb.cell)

head(mat2)

Neworder <- c(comb.cell$V3)
Neworder

corr <- cor(mat2)[Neworder, Neworder]
head(corr)
length(corr)

##breaks and color scale
y <- as.numeric(corr)
##min(y)
##max(y)

##quantile(abs(y),probs=0.99)
##quantile
##y0 <- y[abs(y)< quantile(abs(y),probs=0.99)] #99% percent quantile(abs(y),probs=0.99)
##y0
##seq(0,1,length.out=98)
##mybreaks <- c(min(y),quantile(y0,probs=seq(0,1,length.out=98)),max(y))
##mybreaks <- c(-1,quantile(y0,probs=seq(0,1,length.out=98)),max(y))

mybreaks <- seq(-1, 1, length.out = 100)
mybreaks

names(mybreaks) <- NULL

##
## mycol <- colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100)
## col1 <- c("caffeine"="red", "nicotine"="tan",
##                     "vitA"="tan4", "vitD"="seagreen4",
##                     "vitE"="salmon3", "water"="grey", "zinc"="maroon3")
## col2 <- c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
##                     "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00", "1-TCM"="pink",
##                     "3-TEM"="blue", "5-CD8Naive"="green", "7-dnT"="black")
## tmp_colors <- list(treatment=col1, celltype=col2) #brewer.pal(4,"Set1")
mycol <- colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100)


##
x <- str_split(colnames(corr), "_", simplify=T) %>% as.data.frame()
#head(x)
#colnames(x)
#x <- x[order(x$V2),]
x

tmp_column <- data.frame(celltype=x[,2], treatment=x[,1])
rownames(tmp_column) <- colnames(corr)

tmp_colors <- list(celltype=col2, treatment=col1)

fig2 <- pheatmap(corr, col=mycol,
                 breaks=mybreaks,
                 scale="none", border_color="NA",
                 cluster_rows=F, cluster_cols=F,
                 annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend =T,
                 show_colnames=T, show_rownames=F,
                 fontsize_row=5,
                 fontsize_col=5,
                 na_col="white")

###
figfn <- "./3_treatEffect_2_2.outs/Figure1.2_heatmap.corr.DP.sig.cell_zscore_control_2.png"
png(figfn, width=600, height=600,res=120)
print(fig2)
dev.off()



























#################################################
#----mat data for for top 10 significant only----
#################################################
mat.10topdps.qval <- map_dfr(comb, function(ii){
    res2 <- sig_resDP %>% dplyr::filter(comb==ii)
    res2 <- res2 %>% arrange(p.adjusted)
    res2 <- res2 [1:10, ]
})
print("Number of top 10 significant")
nrow(mat.10topdps.qval)

topdps.qval <- unique(mat.10topdps.qval$peak)
mat.10topdps.qval.2 <- map_dfc(comb, function(ii){
    res2 <- resDP %>% dplyr::filter(comb==ii)
    b <- res2$estimate
    names(b) <- res2$peak
    b[topdps.qval]
})
print("Number of top 10 significant in mat")
nrow(mat.10topdps.qval.2)
mat <- mat.10topdps.qval.2

mat <- as.matrix(mat)
rownames(mat) <- topdps.qval
colnames(mat) <- comb
head(mat)
max(mat)
min(mat)

##
mat[is.na(mat)] = 0
mat2 <- mat

comb.2 <- sort(comb)
mat3 <- mat2[, comb.2]


#----top10 DP plot Fig1.1----
#breaks and color scale
y <- as.numeric(mat2)
max(y)
abs(min(y))

scale <- quantile(abs(y),probs=0.97)
scale

length(y)
length(y[abs(y)>scale])


##y0 <- y[abs(y)< quantile(abs(y),probs=0.99)] #99% percent quantile(abs(y),probs=0.99)
##mybreaks <- c(min(y),quantile(y0,probs=seq(0,1,length.out=98)),max(y))

mybreaks <- seq(-scale, scale, length.out=100)
mybreaks
names(mybreaks) <- NULL

mycol <- colorRampPalette(rev(brewer.pal(n=11, name="RdBu")))(100)
col1 <- c("caffeine"="red", "nicotine"="tan",
                    "vitA"="tan4", "vitD"="seagreen4",
                    "vitE"="salmon3", "water"="grey", "zinc"="maroon3")
col2 <- c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
                    "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00", "1-TCM"="pink",
                    "3-TEM"="blue", "5-CD8Naive"="green", "7-dnT"="black")
tmp_colors <- list(treatment=col1, celltype=col2) #brewer.pal(4,"Set1")

##
x <- str_split(comb, "_", simplify=T)
tmp_column <- data.frame(treatment=x[,1], celltype=x[,2])
rownames(tmp_column) <- comb
fig1 <- pheatmap(mat2,  col=mycol, breaks=mybreaks, border_color="NA",
                 cluster_rows=T, cluster_cols=F, treeheight_row = 0,
                 annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend=T,
                 show_colnames=T,
                 show_rownames=T,
                 na_col="white",
                 fontsize_col=7,
                 fontsize_row=2)
figfn <- "./3_treatEffect_2_2.outs/Figure1.1_heatmap.beta.DP.top10.cell_control_quantile_0.97.png"
png(figfn, width=1800, height=2100,res=225)
print(fig1)
dev.off()
##
x <- str_split(comb.2, "_", simplify=T)
tmp_column <- data.frame(treatment=x[,1], celltype=x[,2])
rownames(tmp_column) <- comb.2
fig1 <- pheatmap(mat3,  col=mycol, breaks=mybreaks, border_color="NA",
                 cluster_rows=T, cluster_cols=F,treeheight_row = 0,
                 annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend=T,
                 show_colnames=T,
                 show_rownames=T,
                 na_col="white",
                 fontsize_col=7,
                 fontsize_row=2)
figfn <- "./3_treatEffect_2_2.outs/Figure1.1_heatmap.beta.DP.top10.treat_control_quantile_0.97.png"
png(figfn, width=1800, height=2100,res=225)
print(fig1)
dev.off()






##----mat_zscore----
mat.10topdps.qval <- map_dfr(comb, function(ii){
    res2 <- sig_resDP %>% dplyr::filter(comb==ii)
    res2 <- res2 %>% arrange(p.adjusted)
    res2 <- res2 [1:10, ]
})
print("Number of top 10 significant")
nrow(mat.10topdps.qval)
topdps.qval <- unique(mat.10topdps.qval$peak)
mat.10topdps.qval.2 <- map_dfc(comb, function(ii){
    res2 <- resDP %>% dplyr::filter(comb==ii)
    b <- res2$statistic
    names(b) <- res2$peak
    b[topdps.qval]
})
print("Number of top 10 significant in mat")
nrow(mat.10topdps.qval.2)
mat <- mat.10topdps.qval.2

mat <- as.matrix(mat)
rownames(mat) <- topdps.qval
colnames(mat) <- comb
head(mat)
max(mat)
min(mat)

#OR
mat[is.na(mat)] = 0
mat2 <- mat

comb.2 <- sort(comb)
mat3 <- mat2[, comb.2]


#----top10 DP plot Fig1.1----
#breaks and color scale
y <- as.numeric(mat2)
max(y)
min(y)

scale <- quantile(abs(y),probs=0.97)
scale

##y0 <- y[abs(y)< quantile(abs(y),probs=0.99)] #99% percent quantile(abs(y),probs=0.99)
##mybreaks <- c(min(y),quantile(y0,probs=seq(0,1,length.out=98)),max(y))

mybreaks <- seq(-scale, scale, length.out=100)
mybreaks
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
x <- str_split(comb, "_", simplify=T)
tmp_column <- data.frame(treatment=x[,1], celltype=x[,2])
rownames(tmp_column) <- comb
fig1 <- pheatmap(mat2,  col=mycol, breaks=mybreaks, border_color="NA",
                 cluster_rows=T, cluster_cols=F, treeheight_row = 0,
                 annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend=T,
                 show_colnames=T,
                 show_rownames=T,
                 na_col="white",
                 fontsize_col=7,
                 fontsize_row=2)
figfn <- "./3_treatEffect_2_2.outs/Figure1.1_heatmap.beta.DP.top10.cell_zscore_control_quantile_0.97.png"
png(figfn, width=1800, height=2100,res=225)
print(fig1)
dev.off()
##
x <- str_split(comb.2, "_", simplify=T)
tmp_column <- data.frame(treatment=x[,1], celltype=x[,2])
rownames(tmp_column) <- comb.2
fig1 <- pheatmap(mat3,  col=mycol, breaks=mybreaks, border_color="NA",
                 cluster_rows=T, cluster_cols=F,treeheight_row = 0,
                 annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend=T,
                 show_colnames=T,
                 show_rownames=T,
                 na_col="white",
                 fontsize_col=7,
                 fontsize_row=2)
figfn <- "./3_treatEffect_2_2.outs/Figure1.1_heatmap.beta.DP.top10.treat_zscore_control_quantile_0.97.png"
png(figfn, width=1800, height=2100,res=225)
print(fig1)
dev.off()



end



## #----plotting as done in motif file----
## #breaks and color scale
## #breaks <- (seq(-8, 8, length.out=17))
## #col_fun <- colorRamp2(breaks, colorRampPalette(rev(brewer.pal(n=11, name="RdBu")))(17))
## #annotation as cell type
## x <- str_split(comb, "_", simplify=T)
## x
## column_ha_ct <- HeatmapAnnotation(
##       treatment=x[,1],
##       celltype=x[,2],
##       col=list(treatment=c("caffeine"="red", "nicotine"="tan",
##                            "vitA"="tan4", "vitD"="seagreen4",
##                            "vitE"="salmon3", "water"="grey", "zinc"="maroon3"),
##                celltype=c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
##                           "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00",
##                           "1-TCM"="pink", "3-TEM"="blue",
##                           "5-CD8Naive"="green", "7-dnT"="black")))
## fig <- Heatmap(mat2,
## #               col=col_fun,
##                cluster_rows=T, cluster_columns=F,show_row_dend = FALSE,
##                top_annotation=column_ha_ct,
##                heatmap_legend_param=list(title="fold.enrichment",
##                                          title_gp=gpar(fontsize=10),
##                                          labels_gp=gpar(fontsize=10)),
##                show_row_names=T,
##                show_column_names=T,
##                row_names_gp=gpar(fontsize=2),
##                column_names_gp=gpar(fontsize=7),
##                raster_device="png")
## figfn <- "./3_treatEffect_2_2.outs/Figure1.1.2_heatmap.beta.DP.top10.cell_zscore.png"
## png(figfn, width=1800, height=2100,res=225)
## print(fig)
## dev.off()
## #annotation as treatment
## x <- str_split(comb.2, "_", simplify=T)
## x
## column_ha_treat <- HeatmapAnnotation(
##       treatment=x[,1],
##       celltype=x[,2],
##       col=list(treatment=c("caffeine"="red", "nicotine"="tan",
##                            "vitA"="tan4", "vitD"="seagreen4",
##                            "vitE"="salmon3", "water"="grey", "zinc"="maroon3"),
##                celltype=c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
##                           "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00",
##                           "1-TCM"="pink", "3-TEM"="blue",
##                           "5-CD8Naive"="green", "7-dnT"="black")))
## fig <- Heatmap(mat3,
## #               col=col_fun,
##                cluster_rows=T, cluster_columns=F, show_row_dend = FALSE,
##                top_annotation=column_ha_treat,
##                heatmap_legend_param=list(title="fold.enrichment",
##                                          title_gp=gpar(fontsize=10),
##                                          labels_gp=gpar(fontsize=10)),
##                show_row_names=T,
##                show_column_names=T,
##                row_names_gp=gpar(fontsize=2),
##                column_names_gp=gpar(fontsize=7),
##                raster_device="png")
## figfn <- "./3_treatEffect_2_2.outs/Figure1.1.2_heatmap.beta.DP.top10.treat_zscore.png"
## png(figfn, width=1800, height=2100,res=225)
## print(fig)
## dev.off()

## end
## #--------

























###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################
#### resDE heatmap   ###
###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################
#----resDE----
resDE <- read_rds("../2_Differential/1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn.rds") %>% as.data.frame()


resDE <- read_rds("../2_Differential/1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn_sampleID+treat_allvsetOH.rds") %>% as.data.frame()


resDE <- read_rds("../2_Differential/1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn_sampleID+new_treat.rds") %>% as.data.frame()


#resDE <- resDE %>% mutate(MCls=gsub("_", "-", MCls)) %>% mutate(comb=paste(MCls, contrast, sep="_"))
#%>% drop_na(p.value)

resDE <- resDE %>% mutate(MCls=gsub("_", "-", MCls)) %>% mutate(comb=paste(contrast, MCls, sep="_"))
#%>% drop_na(p.value)
head(resDE)
nrow(resDE)
max(resDE$estimate)
min(resDE$estimate)


#resDE$Significance <- NULL
#head(resDE)

resDE[!is.na(resDE$p.adjusted) & resDE$p.adjusted<0.1 & resDE$estimate > 0.5, "Significance"] <- "Up"
resDE[!is.na(resDE$p.adjusted) & resDE$p.adjusted<0.1 & resDE$estimate < -0.5, "Significance"] <- "Down"
resDE[!is.na(resDE$p.adjusted) & resDE$p.adjusted>=0.1, "Significance"] <- "Not Significant"
resDE[is.na(resDE$p.adjusted),"Significance"] <- "NA"
head(resDE)
nrow(resDE)
table(resDE$Significance)


resDE[!is.na(resDE$p.adjusted) & resDE$p.adjusted<0.1 & resDE$estimate > 0.25, "Significance.25"] <- "Up"
resDE[!is.na(resDE$p.adjusted) & resDE$p.adjusted<0.1 & resDE$estimate < -0.25, "Significance.25"] <- "Down"
resDE[!is.na(resDE$p.adjusted) & resDE$p.adjusted>=0.1, "Significance.25"] <- "Not Significant"
resDE[is.na(resDE$p.adjusted),"Significance.25"] <- "NA"
#head(resDE)
table(resDE$Significance.25)


resDE[!is.na(resDE$p.adjusted) & resDE$p.adjusted<0.1 & resDE$estimate > 0, "Significance.0"] <- "Up"
resDE[!is.na(resDE$p.adjusted) & resDE$p.adjusted<0.1 & resDE$estimate < 0, "Significance.0"] <- "Down"
resDE[!is.na(resDE$p.adjusted) & resDE$p.adjusted>=0.1, "Significance.0"] <- "Not Significant"
resDE[is.na(resDE$p.adjusted),"Significance.0"] <- "NA"
#head(resDE)
table(resDE$Significance.0)


resDE <- resDE %>% mutate(comb_2=paste(comb, Significance, sep="_"))
head(resDE)
table(resDE$Significance)
length(unique(resDE$gene))


resDE <- resDE %>% mutate(contrast.gene=paste0(contrast, ".", gene))
head(resDE)
nrow(resDE)


#dplyr::filter does not work for vector
#mit.genes <- resDE$gene %>% dplyr::filter(substr(resDE$gene, 1, 2)=="MT")
#the foll works but dont need as doing it other way
#mit.genes <- resDE$gene[substr(resDE$gene, 1, 3)=="MT-"]
#mit.genes

#this select all the rows with "MT" present anywhere in string
#library(data.table)
#mit.resDE <- resDE[resDE$gene %like% "MT", ]
#mit.resDE$gene

resDE <- resDE %>% mutate(mito=ifelse(substr(resDE$gene, 1, 3)=="MT-", 1, 0))
#head(resDE)
nrow(resDE)
table(resDE$mito)


resDE <- resDE %>% mutate(z_score=estimate/stderror)
head(resDE)
nrow(resDE)


#----mitofilt resDE----
resDE <- resDE %>% dplyr::filter(resDE$mito==0)

head(resDE)
nrow(resDE)
length(unique(resDE$gene))
table(resDE$mito)
table(resDE$Significance)
table(resDE$Significance.25)
table(resDE$Significance.0)



#opfn <- paste0("./3_treatEffect_2_2.outs/resDE.mitofilt.rds")

#opfn <- paste0("./3_treatEffect_2_2.outs/resDE_contrastetOH.mitofilt.rds")

opfn <- paste0("./3_treatEffect_2_2.outs/resDE_control.mitofilt.rds")

write_rds(resDE, opfn)


resDE <- read_rds(opfn)

head(resDE)
nrow(resDE)

length(unique(resDE$gene))
table(resDE$mito)
table(resDE$Significance)
table(resDE$Significance.25)
table(resDE$Significance.0)

## resDE.etoh <- read_rds(opfn)
## #head(resDE.etoh)
## nrow(resDE.etoh)
## length(unique(resDE.etoh$gene))
## table(resDE.etoh$mito)
## table(resDE.etoh$Significance)
## table(resDE.etoh$Significance.25)
## table(resDE.etoh$Significance.0)
## resDE <- resDE.etoh

#----all.DE----
all.DE  <- resDE %>% dplyr::pull(gene)
all.DE <- as.character(unique(all.DE))
head(all.DE)
length(all.DE)



## #----sig_resDE----
sig_resDE.5 <- resDE %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)
sig_resDE.25 <- resDE %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<0.1, abs(estimate)>0.25)
sig_resDE.0 <- resDE %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<0.1, abs(estimate)>0)

sig <- table(sig_resDE.5$MCls, sig_resDE.5$contrast)
#sig
sig.matrix <- as.matrix(sig)
sig.melt <- melt(sig.matrix)
sig.p <- ggplot(sig.melt, aes(x = Var2, y = Var1)) +
    geom_tile(aes(fill=value), color="black")+
    scale_fill_gradient(low = "white", high = "white")+
    geom_text(aes(label = value), color = "black", size = 5) +
    ggtitle("Significant_0.1_0.1-0.5") +
    theme(axis.text.x = element_text(angle = 90),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position="none")

figure <- ggarrange(sig.p +
                    font("x.text", size = 14))
                   # ncol = 2, nrow = 2)
png(paste0("./3_treatEffect_2_2.outs/plot_sig_treatandind_filt_all_0.1_0.1_0.5_cn_mitofilt_control.png"), width=500, height=500, pointsize=16, res=125)
print(figure)
dev.off()


sig <- table(sig_resDE.25$MCls, sig_resDE.25$contrast)
sig.matrix <- as.matrix(sig)
sig.melt <- melt(sig.matrix)
sig.p <- ggplot(sig.melt, aes(x = Var2, y = Var1)) +
      geom_tile(aes(fill=value), color="black")+
    scale_fill_gradient(low = "white", high = "white")+
      geom_text(aes(label = value), color = "black", size = 5) +
      ggtitle("Significant_0.1_0.1-0.25") +
      theme(axis.text.x = element_text(angle = 90),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    legend.position="none")
figure <- ggarrange(sig.p +
                    font("x.text", size = 14))
# ncol = 2, nrow = 2)
png(paste0("./3_treatEffect_2_2.outs/plot_sig_treatandind_filt_all_0.1_0.1_0.25_cn_mitofilt.png"), width=500, height=500, pointsize=16, res=125)
print(figure)
dev.off()


sig <- table(sig_resDE.0$MCls, sig_resDE.0$contrast)
sig.matrix <- as.matrix(sig)
sig.melt <- melt(sig.matrix)
sig.p <- ggplot(sig.melt, aes(x = Var2, y = Var1)) +
      geom_tile(aes(fill=value), color="black")+
      scale_fill_gradient(low = "white", high = "white")+
      geom_text(aes(label = value), color = "black", size = 5) +
      ggtitle("Significant_0.1_0.1-0") +
      theme(axis.text.x = element_text(angle = 90),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    legend.position="none")
figure <- ggarrange(sig.p +
                                          font("x.text", size = 14))
# ncol = 2, nrow = 2)
png(paste0("./3_treatEffect_2_2.outs/plot_sig_treatandind_filt_all_0.1_0.1_0_cn_mitofilt.png"), width=500, height=500, pointsize=16, res=125)
print(figure)
dev.off()


sig_resDE <- sig_resDE.5


#----to check----
nrow(sig_resDE %>% dplyr::filter(abs(estimate)>5))


## #----DE----
DE  <- sig_resDE %>% dplyr::pull(gene)
DE <- as.character(unique(DE))
head(DE)
length(DE)




#######################################
#----mat data for significant----
#######################################
head(resDE)

comb <- unique(resDE$comb)
sigdgs.qval <- unique(sig_resDE$gene)

##
mat.sigdgs.qval.2 <- map_dfc(comb, function(ii){
      res2 <- resDE %>% dplyr::filter(comb==ii)
      b <- res2$statistic
      ##b <- res2$estimate
      names(b) <- res2$gene
      b[sigdgs.qval]
})
print("Number of significant in mat")
nrow(mat.sigdgs.qval.2)
mat <- mat.sigdgs.qval.2


mat <- as.matrix(mat)
rownames(mat) <- DE
colnames(mat) <- comb
#head(mat)
nrow(mat)
max(mat)
min(mat)


##
## ii <- rowSums(is.na(mat))
## ii > 0
## head(ii)
## length(ii)
## ##
## #ii <- ii < 56
## #head(ii)
## #length(ii)
## #table (ii)

## ##
## #EITHER
## mat2 <- mat[ii==0,]

## #OR
## mat2 <- mat[ii<length(comb),]
## #mat2[!is.na(mat2)] = 1
## mat2[is.na(mat2)] = 0

#OR
mat[is.na(mat)] = 0
mat2 <- mat

comb.2 <- sort(comb)
mat3 <- mat2[, comb.2]


#----plot----
#breaks and color scale
y <- as.numeric(mat2)
max(y)
abs(min(y))

if(max(y)<abs(min(y))){
    scale <- max(y)
    } else{
        scale <- abs(min(y))
        }
scale

#length(y)

scale <- quantile(abs(y),probs=0.99)
scale


##y0 <- y[abs(y)< quantile(abs(y),probs=0.99)] #99% percent quantile(abs(y),probs=0.99)
#length(y0)

#mybreaks <- c(min(y),quantile(y0,probs=seq(0,1,length.out=98)),max(y))
#length(mybreaks)

mybreaks <- seq(-scale, scale, length.out=100)
mybreaks

names(mybreaks) <- NULL

mycol <- colorRampPalette(rev(brewer.pal(n=11, name="RdBu")))(100)
col1 <- c("caffeine"="red", "nicotine"="tan",
                    "vitA"="tan4", "vitD"="seagreen4",
                    "vitE"="salmon3", "water"="grey", "zinc"="maroon3")
col2 <- c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
                    "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00", "1-TCM"="pink",
                    "3-TEM"="blue", "5-CD8Naive"="green", "7-dnT"="black")
tmp_colors <- list(treatment=col1, celltype=col2) #brewer.pal(4,"Set1")

##
x <- str_split(comb, "_", simplify=T)
tmp_column <- data.frame(treatment=x[,1], celltype=x[,2])
rownames(tmp_column) <- comb
fig1 <- pheatmap(mat2,  col=mycol, breaks=mybreaks, border_color="NA",
                 cluster_rows=T, cluster_cols=F, treeheight_row = 0,
                 annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend=T,
                 show_colnames=T,
                 show_rownames=F,
                 na_col="white",
                 fontsize_col=7,
                 fontsize_row=2)
##figfn <- "./3_treatEffect_2_2.outs/Figure1.1_heatmap.beta.DG.sig.cell_control_2.png"
figfn <- "./3_treatEffect_2_2.outs/Figure1.1_heatmap.beta.DG.sig.25.cell_zscore_control_2_quantile_0.99.png"
png(figfn, width=1800, height=2100,res=225)
print(fig1)
dev.off()
##
x <- str_split(comb.2, "_", simplify=T)
tmp_column <- data.frame(treatment=x[,1], celltype=x[,2])
rownames(tmp_column) <- comb.2
fig1 <- pheatmap(mat3,  col=mycol, breaks=mybreaks, border_color="NA",
                 cluster_rows=T, cluster_cols=F, treeheight_row = 0,
                 annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend=T,
                 show_colnames=T,
                 show_rownames=F,
                 na_col="white",
                 fontsize_col=7,
                 fontsize_row=2)
##figfn <- "./3_treatEffect_2_2.outs/Figure1.1_heatmap.beta.DG.sig.treat_control_2.png"
figfn <- "./3_treatEffect_2_2.outs/Figure1.1_heatmap.beta.DG.sig.treat_zscore_control_2_quantile_0.99.png"
png(figfn, width=1800, height=2100,res=225)
print(fig1)
dev.off()


## #----plotting as done in motif file----
## #breaks and color scale
## #breaks <- (seq(-8, 8, length.out=17))
## #col_fun <- colorRamp2(breaks, colorRampPalette(rev(brewer.pal(n=11, name="RdBu")))(17))
## #----annotation as cell type----
## x <- str_split(comb, "_", simplify=T)
## column_ha_ct <- HeatmapAnnotation(
##       treatment=x[,1],
##       celltype=x[,2],
##       col=list(treatment=c("caffeine"="red", "nicotine"="tan",
##                            "vitA"="tan4", "vitD"="seagreen4",
##                            "vitE"="salmon3", "water"="grey", "zinc"="maroon3"),
##                celltype=c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
##                           "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00",
##                           "1-TCM"="pink", "3-TEM"="blue",
##                           "5-CD8Naive"="green", "7-dnT"="black")))
## fig <- Heatmap(mat2,
## #               col=col_fun,
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
## figfn <- "./3_treatEffect_2_2.outs/Figure1.1.2_heatmap.beta.DG.sig.cell_control.png"
## ##figfn <- "./3_treatEffect_2_2.outs/Figure1.1.2_heatmap.beta.DG.sig.cell_zscore_control.png"
## png(figfn, width=1800, height=2100,res=225)
## print(fig)
## dev.off()
## #----annotation as treatment----
## x <- str_split(comb.2, "_", simplify=T)
## column_ha_treat <- HeatmapAnnotation(
##       treatment=x[,1],
##       celltype=x[,2],
##       col=list(treatment=c("caffeine"="red", "nicotine"="tan",
##                            "vitA"="tan4", "vitD"="seagreen4",
##                            "vitE"="salmon3", "water"="grey", "zinc"="maroon3"),
##                celltype=c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
##                           "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00",
##                           "1-TCM"="pink", "3-TEM"="blue",
##                           "5-CD8Naive"="green", "7-dnT"="black")))
## fig <- Heatmap(mat3,
## #               col=col_fun,
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
## figfn <- "./3_treatEffect_2_2.outs/Figure1.1.2_heatmap.beta.DG.sig.treat_control.png"
## ##figfn <- "./3_treatEffect_2_2.outs/Figure1.1.2_heatmap.beta.DG.sig.treat_zscore_control.png"
## png(figfn, width=1800, height=2100,res=225)
## print(fig)
## dev.off()

## end
## #--------






###########################
### correlation heatmap ###
###########################
#Neworder <- c("Monocyte_LPS", "Monocyte_PHA", "Bcell_LPS", "Bcell_PHA",
   "NKcell_LPS", "NKcell_PHA", "Tcell_LPS", "Tcell_PHA",
   "Monocyte_LPS-DEX", "Monocyte_PHA-DEX", "Bcell_LPS-DEX", "Bcell_PHA-DEX",
   "NKcell_LPS-DEX", "NKcell_PHA-DEX", "Tcell_LPS-DEX", "Tcell_PHA-DEX")

head(comb.2)
head(mat3) #mat3 treat is ordered per treatment

Neworder <- c(comb.2)
Neworder

corr <- cor(mat3, method = "spearman")[Neworder, Neworder]
head(corr)
length(corr)

##breaks and color scale
y <- as.numeric(corr)
##min(y)
##max(y)

##quantile(abs(y),probs=0.99)
##quantile
##y0 <- y[abs(y)< quantile(abs(y),probs=0.99)] #99% percent quantile(abs(y),probs=0.99)
##y0
##seq(0,1,length.out=98)
##mybreaks <- c(min(y),quantile(y0,probs=seq(0,1,length.out=98)),max(y))
##mybreaks <- c(-1,quantile(y0,probs=seq(0,1,length.out=98)),max(y))

mybreaks <- seq(-1, 1, length.out = 100)
mybreaks

names(mybreaks) <- NULL

##
## mycol <- colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100)
## col1 <- c("caffeine"="red", "nicotine"="tan",
##                     "vitA"="tan4", "vitD"="seagreen4",
##                     "vitE"="salmon3", "water"="grey", "zinc"="maroon3")
## col2 <- c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
##                     "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00", "1-TCM"="pink",
##                     "3-TEM"="blue", "5-CD8Naive"="green", "7-dnT"="black")
## tmp_colors <- list(treatment=col1, celltype=col2) #brewer.pal(4,"Set1")
mycol <- colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100)


##
x <- str_split(colnames(corr), "_", simplify=T) %>% as.data.frame()
#head(x)
#colnames(x)
#x <- x[order(x$V2),]
x

tmp_column <- data.frame(celltype=x[,2], treatment=x[,1])
rownames(tmp_column) <- colnames(corr)

tmp_colors <- list(celltype=col2, treatment=col1)

fig2 <- pheatmap(corr, col=mycol,
                 breaks=mybreaks,
                 scale="none", border_color="NA",
                 cluster_rows=F, cluster_cols=F,
                 annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend =T,
                 show_colnames=T, show_rownames=F,
                 fontsize_row=5,
                 fontsize_col=5,
                 na_col="white")

###
##figfn <- "./3_treatEffect_2_2.outs/Figure1.2_heatmap.corr.DG.sig.treat_control_2.png"
figfn <- "./3_treatEffect_2_2.outs/Figure1.2_heatmap.corr.DG.sig.treat_zscore_control_2.png"
png(figfn, width=600, height=600,res=120)
print(fig2)
dev.off()



##----as celltype----
comb.cell <- str_split(comb.2, "_", simplify=T) %>% as.data.frame()
comb.cell <- comb.cell[order(comb.cell$V2),] %>% mutate(V3 = paste0(V1, "_", V2))
head(comb.cell)

head(mat2)

Neworder <- c(comb.cell$V3)
Neworder

corr <- cor(mat2)[Neworder, Neworder]
head(corr)
length(corr)

##breaks and color scale
y <- as.numeric(corr)
##min(y)
##max(y)

##quantile(abs(y),probs=0.99)
##quantile
##y0 <- y[abs(y)< quantile(abs(y),probs=0.99)] #99% percent quantile(abs(y),probs=0.99)
##y0
##seq(0,1,length.out=98)
##mybreaks <- c(min(y),quantile(y0,probs=seq(0,1,length.out=98)),max(y))
##mybreaks <- c(-1,quantile(y0,probs=seq(0,1,length.out=98)),max(y))

mybreaks <- seq(-1, 1, length.out = 100)
mybreaks

names(mybreaks) <- NULL

##
## mycol <- colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100)
## col1 <- c("caffeine"="red", "nicotine"="tan",
##                     "vitA"="tan4", "vitD"="seagreen4",
##                     "vitE"="salmon3", "water"="grey", "zinc"="maroon3")
## col2 <- c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
##                     "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00", "1-TCM"="pink",
##                     "3-TEM"="blue", "5-CD8Naive"="green", "7-dnT"="black")
## tmp_colors <- list(treatment=col1, celltype=col2) #brewer.pal(4,"Set1")
mycol <- colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100)


##
x <- str_split(colnames(corr), "_", simplify=T) %>% as.data.frame()
#head(x)
#colnames(x)
#x <- x[order(x$V2),]
x

tmp_column <- data.frame(celltype=x[,2], treatment=x[,1])
rownames(tmp_column) <- colnames(corr)

tmp_colors <- list(celltype=col2, treatment=col1)

fig2 <- pheatmap(corr, col=mycol,
                 breaks=mybreaks,
                 scale="none", border_color="NA",
                 cluster_rows=F, cluster_cols=F,
                 annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend =T,
                 show_colnames=T, show_rownames=F,
                 fontsize_row=5,
                 fontsize_col=5,
                 na_col="white")

###
##figfn <- "./3_treatEffect_2_2.outs/Figure1.2_heatmap.corr.DG.sig.cell_control_2.png"
figfn <- "./3_treatEffect_2_2.outs/Figure1.2_heatmap.corr.DG.sig.cell_zscore_control_2.png"
png(figfn, width=600, height=600,res=120)
print(fig2)
dev.off()




## #----plotting as done in motif file----
## #breaks and color scale
## #breaks <- (seq(-8, 8, length.out=17))
## #col_fun <- colorRamp2(breaks, colorRampPalette(rev(brewer.pal(n=11, name="RdBu")))(17))
## #annotation as treatment
## x <- str_split(comb.2, "_", simplify=T)
## column_ha_treat <- HeatmapAnnotation(
##       treatment=x[,1],
##       celltype=x[,2],
##       col=list(treatment=c("caffeine"="red", "nicotine"="tan",
##                            "vitA"="tan4", "vitD"="seagreen4",
##                            "vitE"="salmon3", "water"="grey", "zinc"="maroon3"),
##                celltype=c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
##                           "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00",
##                           "1-TCM"="pink", "3-TEM"="blue",
##                           "5-CD8Naive"="green", "7-dnT"="black")))
## fig <- Heatmap(corr,
## #               col=col_fun,
##                cluster_rows=F, cluster_columns=F, show_row_dend = FALSE,
##                top_annotation=column_ha_treat,
##                heatmap_legend_param=list(title="fold.enrichment",
##                                          title_gp=gpar(fontsize=10),
##                                          labels_gp=gpar(fontsize=10)),
##                show_row_names=T,
##                show_column_names=T,
##                row_names_gp=gpar(fontsize=2),
##                column_names_gp=gpar(fontsize=7),
##                raster_device="png")
## figfn <- "./3_treatEffect_2_2.outs/Figure1.2.2_heatmap.corr.DG.sig.treat_control.png"
## ##figfn <- "./3_treatEffect_2_2.outs/Figure1.2.2_heatmap.corr.DG.sig.treat_zscore_control.png"
## png(figfn, width=1800, height=2100,res=225)
## print(fig)
## dev.off()

## end
## #--------
















#################################################
#----mat data for for top 10 significant only----
#################################################
mat.10topdgs.qval <- map_dfr(comb, function(ii){
      res2 <- sig_resDE %>% dplyr::filter(comb==ii)
      res2 <- res2 %>% arrange(p.adjusted)
      res2 <- res2 [1:10, ]
      })
print("Number of top 10 significant")
nrow(mat.10topdgs.qval)

head(mat.10topdgs.qval)

write.csv(mat.10topdgs.qval, "./2.2_compareRNAandATAC_3.outs/DG.top10_control.csv")

topdgs.qval <- unique(mat.10topdgs.qval$gene)

mat.10topdgs.qval.2 <- map_dfc(comb, function(ii){
    res2 <- resDE %>% dplyr::filter(comb==ii)
    b <- res2$statistic
    ##b <- res2$estimate
    names(b) <- res2$gene
    b[topdgs.qval]
      })
print("Number of top 10 significant in mat")
nrow(mat.10topdgs.qval.2)
mat <- mat.10topdgs.qval.2

mat <- as.matrix(mat)
rownames(mat) <- topdgs.qval
colnames(mat) <- comb
head(mat)
max(mat)
min(mat)


##
## ii <- rowSums(is.na(mat))
## ii > 0
## head(ii)
## length(ii)
## ##
## #ii <- ii < 56
## #head(ii)
## #length(ii)
## #table (ii)
## ##

## #EITHER
## mat2 <- mat[ii==0,]

## #OR
## mat2 <- mat[ii<length(comb),]
## #mat2[!is.na(mat2)] = 1
## mat2[is.na(mat2)] = 0

#OR
mat[is.na(mat)] = 0
mat2 <- mat

comb.2 <- sort(comb)
mat3 <- mat2[, comb.2]





## #----plot in loop----
## #breaks and color scale
## y <- as.numeric(mat2)
## y

## abs(max(y))
## abs(min(y))

## length(y[y < -9])

## scale <- max(y)
## #scale <- abs(min(y))

## ##quantile(abs(y),probs=0.99)
## ##y0 <- y[abs(y)< quantile(abs(y),probs=0.99)] #99% percent quantile(abs(y),probs=0.99)

## ##mybreaks <- c(min(y),quantile(y0,probs=seq(0,1,length.out=98)),max(y))


## for(i in c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)){
##     mybreaks <- seq(-scale, scale, length.out=i)
##     ##mybreaks
##     names(mybreaks) <- NULL
##     mycol <- colorRampPalette(rev(brewer.pal(n=11, name="RdBu")))(50)
##     col1 <- c("caffeine"="red", "nicotine"="tan",
##               "vitA"="tan4", "vitD"="seagreen4",
##               "vitE"="salmon3", "water"="grey", "zinc"="maroon3")
##     col2 <- c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
##               "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00", "1-TCM"="pink",
##               "3-TEM"="blue", "5-CD8Naive"="green", "7-dnT"="black")
##     tmp_colors <- list(treatment=col1, celltype=col2) #brewer.pal(4,"Set1")
##     ##
##     x <- str_split(comb, "_", simplify=T)
##     tmp_column <- data.frame(treatment=x[,1], celltype=x[,2])
##     rownames(tmp_column) <- comb
##     fig1 <- pheatmap(mat2,  col=mycol, breaks=mybreaks, border_color="NA",
##                      cluster_rows=T, cluster_cols=F, treeheight_row = 0,
##                      annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend=T,
##                      show_colnames=T,
##                      show_rownames=T,
##                      na_col="white",
##                      fontsize_col=7,
##                      fontsize_row=2)
##     figfn <- paste0("./3_treatEffect_2_2.outs/control3/Figure1.1_heatmap.beta.DG.top10.cell_", i, "_cg-50_control_3.png")
##     ##figfn <- "./3_treatEffect_2_2.outs/control3/Figure1.1_heatmap.beta.DG.top10.cell_zscore_control_3.png"
##     png(figfn, width=1800, height=2100,res=225)
##     print(fig1)
##     dev.off()
## }





#----plot----
#breaks and color scale
y <- as.numeric(mat2)
length(y)
max(y)
abs(min(y))


scale <- quantile(abs(y),probs=0.99)
scale
length(y[abs(y)>scale])



## if(max(y)<abs(min(y))){
##     scale <- max(y)
##     } else{
##         scale <- abs(min(y))
##         }
## scale



##quantile(abs(y),probs=0.99)
####y <- y[abs(y)< scale]
##y0 <- y[abs(y)< quantile(abs(y),probs=0.99)] #99% percent quantile(abs(y),probs=0.99)
##mybreaks <- c(min(y),quantile(y0,probs=seq(0,1,length.out=98)),max(y))
##mybreaks <- c(-scale,quantile(y0,probs=seq(0,1,length.out=98)),scale)
##mybreaks

##mybreaks <- seq(min(y), max(y), length.out=100)
mybreaks <- seq(-scale, scale, length.out=100)
mybreaks

names(mybreaks) <- NULL

mycol <- colorRampPalette(rev(brewer.pal(n=11, name="RdBu")))(100)
col1 <- c("caffeine"="red", "nicotine"="tan",
                    "vitA"="tan4", "vitD"="seagreen4",
                    "vitE"="salmon3", "water"="grey", "zinc"="maroon3")
col2 <- c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
                    "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00", "1-TCM"="pink",
                    "3-TEM"="blue", "5-CD8Naive"="green", "7-dnT"="black")
tmp_colors <- list(treatment=col1, celltype=col2) #brewer.pal(4,"Set1")
##
x <- str_split(comb, "_", simplify=T)
tmp_column <- data.frame(treatment=x[,1], celltype=x[,2])
rownames(tmp_column) <- comb
fig1 <- pheatmap(mat2,  col=mycol, breaks=mybreaks, border_color="NA",
                 cluster_rows=T, cluster_cols=F, treeheight_row = 0,
                 annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend=T,
                 show_colnames=T,
                 show_rownames=T,
                 na_col="white",
                 fontsize_col=7,
                 fontsize_row=2)
##figfn <- "./3_treatEffect_2_2.outs/Figure1.1_heatmap.beta.DG.top10.cell_control_2_quantile_0.99.png"
figfn <- "./3_treatEffect_2_2.outs/Figure1.1_heatmap.beta.DG.25.top10.cell_zscore_control_2_quantile_0.99.png"
png(figfn, width=1800, height=2100,res=225)
print(fig1)
dev.off()
##
x <- str_split(comb.2, "_", simplify=T)
tmp_column <- data.frame(treatment=x[,1], celltype=x[,2])
rownames(tmp_column) <- comb.2
fig1 <- pheatmap(mat3,  col=mycol, breaks=mybreaks, border_color="NA",
                 cluster_rows=T, cluster_cols=F, treeheight_row = 0,
                 annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend=T,
                 show_colnames=T,
                 show_rownames=T,
                 na_col="white",
                 fontsize_col=7,
                 fontsize_row=2)
##figfn <- "./3_treatEffect_2_2.outs/Figure1.1_heatmap.beta.DG.top10.treat_control_2_quantile_0.99.png"
figfn <- "./3_treatEffect_2_2.outs/Figure1.1_heatmap.beta.DG.25.top10.treat_zscore_control_2_quantile_0.99.png"
png(figfn, width=1800, height=2100,res=225)
print(fig1)
dev.off()


#----plotting as done in motif file----
#breaks and color scale
#breaks <- (seq(-8, 8, length.out=17))
#col_fun <- colorRamp2(breaks, colorRampPalette(rev(brewer.pal(n=11, name="RdBu")))(17))
#annotation as cell type
x <- str_split(comb, "_", simplify=T)
column_ha_ct <- HeatmapAnnotation(
      treatment=x[,1],
      celltype=x[,2],
      col=list(treatment=c("caffeine"="red", "nicotine"="tan",
                           "vitA"="tan4", "vitD"="seagreen4",
                           "vitE"="salmon3", "water"="grey", "zinc"="maroon3"),
               celltype=c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
                          "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00",
                          "1-TCM"="pink", "3-TEM"="blue",
                          "5-CD8Naive"="green", "7-dnT"="black")))
fig <- Heatmap(mat2,
#               col=col_fun,
               cluster_rows=T, cluster_columns=F, show_row_dend = FALSE,
               top_annotation=column_ha_ct,
               heatmap_legend_param=list(title="fold.enrichment",
                                         title_gp=gpar(fontsize=10),
                                         labels_gp=gpar(fontsize=10)),
               show_row_names=T,
               show_column_names=T,
               row_names_gp=gpar(fontsize=2),
               column_names_gp=gpar(fontsize=7),
               raster_device="png")
figfn <- "./3_treatEffect_2_2.outs/Figure1.1.2_heatmap.beta.DG.top10.cell_control.png"
##figfn <- "./3_treatEffect_2_2.outs/Figure1.1.2_heatmap.beta.DG.top10.cell_zscore_control.png"
png(figfn, width=1800, height=2100,res=225)
print(fig)
dev.off()
#annotation as treatment
x <- str_split(comb.2, "_", simplify=T)
column_ha_treat <- HeatmapAnnotation(
      treatment=x[,1],
      celltype=x[,2],
      col=list(treatment=c("caffeine"="red", "nicotine"="tan",
                           "vitA"="tan4", "vitD"="seagreen4",
                           "vitE"="salmon3", "water"="grey", "zinc"="maroon3"),
               celltype=c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
                          "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00",
                          "1-TCM"="pink", "3-TEM"="blue",
                          "5-CD8Naive"="green", "7-dnT"="black")))
fig <- Heatmap(mat3,
#               col=col_fun,
               cluster_rows=T, cluster_columns=F, show_row_dend = FALSE,
               top_annotation=column_ha_treat,
               heatmap_legend_param=list(title="fold.enrichment",
                                         title_gp=gpar(fontsize=10),
                                         labels_gp=gpar(fontsize=10)),
               show_row_names=T,
               show_column_names=T,
               row_names_gp=gpar(fontsize=2),
               column_names_gp=gpar(fontsize=7),
               raster_device="png")
figfn <- "./3_treatEffect_2_2.outs/Figure1.1.2_heatmap.beta.DG.top10.treat_control.png"
##figfn <- "./3_treatEffect_2_2.outs/Figure1.1.2_heatmap.beta.DG.top10.treat_zscore_control.png"
png(figfn, width=1800, height=2100,res=225)
print(fig)
dev.off()

end
#--------





