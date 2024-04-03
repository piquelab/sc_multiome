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
outdir <- "./3_treatEffect_4_countmatrix.outs/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=F)


## #####################################################################
## ###            heatmap for resDP                                  ###
## #####################################################################
## #----resDP----
## resDP <- read_rds("../2_Differential/1.2_DiffPeak.outs/2.0_DESeq.results_0.1_cn.rds") %>% as.data.frame()
## nrow(resDP)
## nrow(resDP %>% dplyr::filter(is.na(p.adjusted)))

## #resDP <- resDP %>% mutate(MCls=gsub("_", "-", MCls)) %>% mutate(comb=paste(MCls, contrast, sep="_"))%>% dplyr::rename("peak"="gene")
## #%>% drop_na(p.value)
## resDP <- resDP %>% mutate(MCls=gsub("_", "-", MCls)) %>% mutate(comb=paste(contrast, MCls, sep="_"))%>% dplyr::rename("peak"="gene")
## #%>% drop_na(p.value)
## head(resDP)
## nrow(resDP)
## max(resDP$estimate)
## min(resDP$estimate)


## resDP[!is.na(resDP$p.adjusted) & resDP$p.adjusted<0.1 & resDP$estimate > 0.5, "Significance"] <- "Up"
## resDP[!is.na(resDP$p.adjusted) & resDP$p.adjusted<0.1 & resDP$estimate < -0.5, "Significance"] <- "Down"
## resDP[!is.na(resDP$p.adjusted) & resDP$p.adjusted>=0.1, "Significance"] <- "Not Significant"
## resDP[is.na(resDP$p.adjusted),"Significance"] <- "NA"
## resDP <- resDP %>% mutate(comb_2=paste(comb, Significance, sep="_"))
## head(resDP)
## table(resDP$Significance)
## table(resDP$MCls, resDP$contrast)


## resDP[!is.na(resDP$p.adjusted) & resDP$p.adjusted<0.1 & resDP$estimate > 0.25, "Significance.25"] <- "Up"
## resDP[!is.na(resDP$p.adjusted) & resDP$p.adjusted<0.1 & resDP$estimate < -0.25, "Significance.25"] <- "Down"
## resDP[!is.na(resDP$p.adjusted) & resDP$p.adjusted>=0.1, "Significance.25"] <- "Not Significant"
## resDP[is.na(resDP$p.adjusted),"Significance.25"] <- "NA"
## #head(resDP)
## table(resDP$Significance.25)


## resDP[!is.na(resDP$p.adjusted) & resDP$p.adjusted<0.1 & resDP$estimate > 0, "Significance.0"] <- "Up"
## resDP[!is.na(resDP$p.adjusted) & resDP$p.adjusted<0.1 & resDP$estimate < 0, "Significance.0"] <- "Down"
## resDP[!is.na(resDP$p.adjusted) & resDP$p.adjusted>=0.1, "Significance.0"] <- "Not Significant"
## resDP[is.na(resDP$p.adjusted),"Significance.0"] <- "NA"
## #head(resDP)
## table(resDP$Significance.0)

## opfn <- paste0("./3_treatEffect_2_2.outs/resDP.rds")

## write_rds(resDP, opfn)

## resDP <- read_rds(opfn)


## #----allDP----
## all.DP  <- resDP %>% dplyr::pull(peak)
## all.DP <- as.character(unique(all.DP))
## head(all.DP)
## length(all.DP)


## #----sig_resDP----
## sig_resDP.5 <- resDP %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)
## sig_resDP.25 <- resDP %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<0.1, abs(estimate)>0.25)
## sig_resDP.0 <- resDP %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<0.1, abs(estimate)>0)

## sig <- table(sig_resDP.5$MCls, sig_resDP.5$contrast)
## sig.matrix <- as.matrix(sig)
## sig.melt <- melt(sig.matrix)
## sig.p <- ggplot(sig.melt, aes(x = Var2, y = Var1)) +
##       geom_tile(aes(fill=value), color="black")+
##       scale_fill_gradient(low = "white", high = "white")+
##       geom_text(aes(label = value), color = "black", size = 5) +
##       ggtitle("Significant_0.1_0.1-0.5") +
##       theme(axis.text.x = element_text(angle = 90),
##                     axis.title.x = element_blank(),
##                     axis.title.y = element_blank(),
##                     legend.position="none")
## figure <- ggarrange(sig.p +
##                     font("x.text", size = 14))
## # ncol = 2, nrow = 2)
## png(paste0("./3_treatEffect_2_2.outs/plot_sig_treatandind_filt_DP_all_0.1_0.1_0.5_cn_mitofilt.png"), width=500, height=500, pointsize=16, res=125)
## print(figure)
## dev.off()




## sig <- table(sig_resDP.25$MCls, sig_resDP.25$contrast)
## sig.matrix <- as.matrix(sig)
## sig.melt <- melt(sig.matrix)
## sig.p <- ggplot(sig.melt, aes(x = Var2, y = Var1)) +
##       geom_tile(aes(fill=value), color="black")+
##       scale_fill_gradient(low = "white", high = "white")+
##       geom_text(aes(label = value), color = "black", size = 5) +
##       ggtitle("Significant_0.1_0.1-0.25") +
##       theme(axis.text.x = element_text(angle = 90),
##                     axis.title.x = element_blank(),
##                     axis.title.y = element_blank(),
##                     legend.position="none")
## figure <- ggarrange(sig.p +
##                                           font("x.text", size = 14))
## # ncol = 2, nrow = 2)
## png(paste0("./3_treatEffect_2_2.outs/plot_sig_treatandind_filt_DP_all_0.1_0.1_0.25_cn_mitofilt.png"), width=500, height=500, pointsize=16, res=125)
## print(figure)
## dev.off()




## sig <- table(sig_resDP.0$MCls, sig_resDP.0$contrast)
## sig.matrix <- as.matrix(sig)
## sig.melt <- melt(sig.matrix)
## sig.p <- ggplot(sig.melt, aes(x = Var2, y = Var1)) +
##       geom_tile(aes(fill=value), color="black")+
##       scale_fill_gradient(low = "white", high = "white")+
##       geom_text(aes(label = value), color = "black", size = 5) +
##       ggtitle("Significant_0.1_0.1-0") +
##       theme(axis.text.x = element_text(angle = 90),
##                     axis.title.x = element_blank(),
##                     axis.title.y = element_blank(),
##                     legend.position="none")
## figure <- ggarrange(sig.p +
##                                           font("x.text", size = 14))
## # ncol = 2, nrow = 2)
## png(paste0("./3_treatEffect_2_2.outs/plot_sig_treatandind_filt_DP_all_0.1_0.1_0_cn_mitofilt.png"), width=500, height=500, pointsize=16, res=125)
## print(figure)
## dev.off()


## sig_resDP <- sig_resDP.5


## #----to check----
## nrow(sig_resDP %>% dplyr::filter(abs(estimate)>5))


## #----DP-----
## DP  <- sig_resDP %>% dplyr::pull(peak)
## DP <- as.character(unique(DP))
## head(DP)
## length(DP)




########################################
#----count matrices for ATAC and RNA----
########################################
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(Signac)

pbmc <- read_rds("../1_processing/7.1_Clustering_macs2.outs/1_seurat.cluster.combined.mitofilt.rds")

head(pbmc@meta.data)
colnames(pbmc@meta.data)

YtX_rna <- read_rds("./1_DiffRNA_2.outs/1_YtX.sel_clusters_0.1_cn.rds")
head(YtX_rna)


bti2 <- colnames(YtX_rna)
#length(bti2)

cvt0 <- str_split(bti2, "_", simplify=T)
#head(cvt0)

cvt <- data.frame(bti=bti2, MCls=paste(cvt0[,1], cvt0[,2], sep="_"), treat=cvt0[,3], sampleID=cvt0[,4])
head(cvt)

#comb <- table(cvt$MCls, cvt$treat)
#write.csv(comb, paste0("./1_DiffRNA_2.outs/combinations.csv"), quote = FALSE, row.names = TRUE)
#comb <- read.csv("./1_DiffRNA_2.outs/combinations.csv")
#comb
#dd2 <- dd %>% dplyr::filter(bti %in% bti2)
#sum(dd2$ncell)


cvt <- cvt%>%mutate(MCls_IDs=paste(cvt$MCls, ".", cvt$sampleID, sep=""))
#head(cvt)

MCls.IDs <- unique(cvt$MCls_IDs)
#MCls.IDs

cvt.filtered <- map_dfr(MCls.IDs, function(oneX){
          time0 <- Sys.time()
                  ###
                  cvt.2 <- cvt%>%dplyr::filter(MCls==str_split(oneX, "\\.", simplify=T)[1], sampleID==str_split(oneX, "\\.", simplify=T)[2])
                  ##
                  cvt.2 <- cvt.2%>%mutate(count=ifelse(nrow(cvt.2)>=8,1,0))
                  ##
                  time1 <- Sys.time()
                  elapsed <- difftime(time1, time0, units="mins")
                  cat(oneX, elapsed, "Done\n")
                  cvt.2
                })
#cvt.filtered

cvt.2 <- cvt.filtered %>% dplyr::filter(count==1)
#head(cvt.2)

comb.2 <- table(cvt.2$MCls, cvt.2$treat)
comb.2
#write.csv(comb.2, paste0("./1_DiffRNA_2.outs/combinations_2_", reso, "_cn.csv"), quote = FALSE, row.names = TRUE)
#comb.2 <- read.csv("./1_DiffRNA_2.outs/combinations_2.csv")
#comb.2

bti3 <- cvt.2$bti
#bti3

#ncol(YtX)
YtX_rna2 <- YtX_rna[,bti3]
#ncol(YtX_rna2)
head(YtX_rna2)

head(cvt.2)

cvt.3 <- cvt.2 %>% mutate(MCls = gsub("_", "-", MCls)) %>% mutate(MCls_treat = paste0(MCls, "_", treat))
head(cvt.3)

bti4 <- factor(cvt.3$MCls_treat)
sort(unique(bti4))

X <- model.matrix(~0+bti4)
#head(X)

YtX_rna3 <- YtX_rna2 %*% X
YtX_rna3 <- as.matrix(YtX_rna3)
colnames(YtX_rna3) <- gsub("^bti4", "", colnames(YtX_rna3))
head(YtX_rna3)
nrow(YtX_rna3)
ncol(YtX_rna3)


mat_rna <- YtX_rna3

comb_rna <- colnames(mat_rna)

x <- str_split(comb_rna, "_", simplify=T) %>% as.data.frame()
x <- x %>% arrange(x$V2) %>% mutate(V3 = paste0(V1, "_", V2))
head(x)

comb_rna.2 <- x$V3

mat3_rna <- mat_rna[, comb_rna.2]
head(mat3_rna)



#----corr----
Neworder <- c(comb_rna.2)
Neworder

corr <- cor(mat3_rna, method =c("spearman"))[Neworder, Neworder]

corr <- cor(mat3_rna)[Neworder, Neworder]

head(corr)

mycol <- colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100)
##

###
x <- str_split(colnames(corr), "_", simplify=T)
head(x)

tmp_column <- data.frame(celltype=x[,1], treatment=x[,2])
rownames(tmp_column) <- colnames(corr)

col1 <- c("caffeine"="red", "nicotine"="tan",
          "vitA"="tan4", "vitD"="seagreen4",
          "vitE"="salmon3", "water"="grey", "zinc"="maroon3", "etOH"="grey")
col2 <- c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
          "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00", "1-TCM"="pink",
          "3-TEM"="blue", "5-CD8Naive"="green", "7-dnT"="black")

tmp_colors <- list(celltype=col2, treatment=col1)

fig2 <- pheatmap(corr, col=mycol, scale="none", border_color="NA",
                 cluster_rows=F, cluster_cols=F,
                 annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend =T,
                 show_colnames=T, show_rownames=T,
                 fontsize_row=5,
                 fontsize_col=5,
                 na_col="white")

###
figfn <- "./3_treatEffect_4_countmatrix.outs/Figure1.2_heatmap.corr.DG.treat.count.pearson.png"
png(figfn, width=600, height=600,res=120)
print(fig2)
dev.off()




end
















#----ATAC----
YtX_atac <- read_rds("./1.2_DiffPeak.outs/1_YtX.sel_0.1_cn.rds")
nrow(YtX_atac)
ncol(YtX_atac)

bti2 <- colnames(YtX_atac)
#bti2

cvt0 <- str_split(bti2, "_", simplify=T)
#cvt0

cvt <- data.frame(bti=bti2, MCls=paste(cvt0[,1], cvt0[,2], sep="_"), treat=cvt0[,3], sampleID=cvt0[,4])
#cvt
# table(cvt$treat)
# table(cvt$MCls)
# table(cvt$MCls, cvt$treat)

# comb <- table(cvt$MCls, cvt$treat)
# write.csv(comb, paste0("./1_DiffPeak.outs/combinations.csv"), quote = FALSE, row.names = TRUE)
# comb <- read.csv("./1_DiffPeak.outs/combinations.csv")
# comb
# dd
# dd2 <- dd %>% dplyr::filter(bti %in% bti2)
# dd2
# sum(dd2$ncell)

cvt <- cvt%>%mutate(MCls_IDs=paste(cvt$MCls, ".", cvt$sampleID, sep=""))
#head(cvt)

MCls.IDs <- unique(cvt$MCls_IDs)
#MCls.IDs

cvt.filtered <- map_dfr(MCls.IDs, function(oneX){
          time0 <- Sys.time()
                  ###
                  cvt.2 <- cvt%>%dplyr::filter(MCls==str_split(oneX, "\\.", simplify=T)[1], sampleID==str_split(oneX, "\\.", simplify=T)[2])
                  ##
                  cvt.2 <- cvt.2%>%mutate(count=ifelse(nrow(cvt.2)>=8,1,0))
                  ##
                  time1 <- Sys.time()
                  elapsed <- difftime(time1, time0, units="mins")
                  cat(oneX, elapsed, "Done\n")
                  cvt.2
                })
#cvt.filtered

cvt.2 <- cvt.filtered %>% dplyr::filter(count==1)

comb.2 <- table(cvt.2$MCls, cvt.2$treat)
comb.2
#write.csv(comb.2, paste0("./1.2_DiffPeak.outs/combinations_2_", reso, "_cn.csv"), quote = FALSE, row.names = TRUE)
#comb.2 <- read.csv("./1.2_DiffPeak.outs/combinations_2.csv")
#comb.2

bti3 <- cvt.2$bti
#bti3

#ncol(YtX)
YtX_atac2 <- YtX_atac[,bti3]
ncol(YtX_atac2)
head(YtX_atac2)

cvt.3 <- cvt.2 %>% mutate(MCls = gsub("_", "-", MCls))  %>% mutate(MCls_treat = paste0(MCls, "_", treat))
head(cvt.3)

bti4 <- factor(cvt.3$MCls_treat)
sort(unique(bti4))

X <- model.matrix(~0+bti4)
#head(X)

YtX_atac3 <- YtX_atac2 %*% X
YtX_atac3 <- as.matrix(YtX_atac3)
colnames(YtX_atac3) <- gsub("^bti4", "", colnames(YtX_atac3))
head(YtX_atac3)
nrow(YtX_atac3)
ncol(YtX_atac3)

mat_atac <- YtX_atac3

comb_atac <- colnames(mat_atac)

x <- str_split(comb_atac, "_", simplify=T) %>% as.data.frame()
x <- x %>% arrange(x$V2) %>% mutate(V3 = paste0(V1, "_", V2))
head(x)

comb_atac.2 <- x$V3

mat3_atac <- mat_atac[, comb_atac.2]
head(mat3_atac)



#----corr----
Neworder <- c(comb_atac.2)
Neworder

corr <- cor(mat3_atac, method =c("spearman"))[Neworder, Neworder]

corr <- cor(mat3_atac)[Neworder, Neworder]

head(corr)

mycol <- colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100)
##

###
x <- str_split(colnames(corr), "_", simplify=T)
head(x)

tmp_column <- data.frame(celltype=x[,1], treatment=x[,2])
rownames(tmp_column) <- colnames(corr)

col1 <- c("caffeine"="red", "nicotine"="tan",
          "vitA"="tan4", "vitD"="seagreen4",
          "vitE"="salmon3", "water"="grey", "zinc"="maroon3", "etOH"="grey")
col2 <- c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
          "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00", "1-TCM"="pink",
          "3-TEM"="blue", "5-CD8Naive"="green", "7-dnT"="black")

tmp_colors <- list(celltype=col2, treatment=col1)

fig2 <- pheatmap(corr, col=mycol, scale="none", border_color="NA",
                 cluster_rows=F, cluster_cols=F,
                 annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend =T,
                 show_colnames=T, show_rownames=T,
                 fontsize_row=5,
                 fontsize_col=5,
                 na_col="white")

###
figfn <- "./3_treatEffect_4_countmatrix.outs/Figure1.2_heatmap.corr.DP.treat.count.pearson.png"
png(figfn, width=600, height=600,res=120)
print(fig2)
dev.off()




end















