###
library(Matrix)
library(tidyverse)
library(qqman)
library(qvalue)
##
library(DESeq2)
library(biobroom)
library(ashr)
library(GenomicRanges)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(Signac)
### annotation required package
library(clusterProfiler)
library(org.Hs.eg.db)
library(annotables)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(EnsDb.Hsapiens.v86)
library(ChIPseeker)
#library(ChIPseeker, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
###
library(ggplot2)
library(cowplot, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(grid)
library(gridExtra)
library(ggExtra)
library(RColorBrewer)
library(ggsci)
library(RColorBrewer)
library(viridis)
theme_set(theme_grey())
##
library(ComplexHeatmap)
library(circlize)

outdir <- "./2.2_compareRNAandATAC_4_linkpeaks_sep.outs/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)


## ####################################################
## ### peak annotattion using linkpeakstogenes - link
## ####################################################
#link <- read_rds("../1_processing/6.1_linkPeakstoGenes.outs/links_dist100000.rds")
#head(link)

#dat.link <- link %>% as.data.frame()
#head(dat.link)
#nrow(dat.link)

dat.caffeine <- readRDS(paste0("../1_processing/6.1_linkPeakstoGenes_2.outs/links_caffeine_dist100000.rds"))
dat.nicotine <- readRDS(paste0("../1_processing/6.1_linkPeakstoGenes_2.outs/links_nicotine_dist100000.rds"))
dat.zinc <- readRDS(paste0("../1_processing/6.1_linkPeakstoGenes_2.outs/links_zinc_dist100000.rds"))
dat.vitA <- readRDS(paste0("../1_processing/6.1_linkPeakstoGenes_2.outs/links_vitA_dist100000.rds"))
dat.vitD <- readRDS(paste0("../1_processing/6.1_linkPeakstoGenes_2.outs/links_vitD_dist100000.rds"))
dat.vitE <- readRDS(paste0("../1_processing/6.1_linkPeakstoGenes_2.outs/links_vitE_dist100000.rds"))
#dat.water <- readRDS(paste0("../1_processing/6.1_linkPeakstoGenes_2.outs/links_water_dist100000.rds"))
#dat.etOH <- readRDS(paste0("../1_processing/6.1_linkPeakstoGenes_2.outs/links_etOH_dist100000.rds"))


dat.caffeine$treats <- "caffeine"
dat.nicotine$treats <- "nicotine"
dat.zinc$treats <- "zinc"
dat.vitA$treats <- "vitA"
dat.vitD$treats <- "vitD"
dat.vitE$treats <- "vitE"
#dat.water$treats <- "water"
#dat.etOH$treats <- "etOH" 


dat <- rbind(dat.caffeine, dat.nicotine, dat.vitA, dat.vitD, dat.vitE, dat.zinc)
head(dat)
nrow(dat)


dat2.caffeine <- dat.caffeine %>% arrange(desc(zscore)) %>% distinct(peak, .keep_all = TRUE)
dat2.nicotine <- dat.nicotine %>% arrange(desc(zscore)) %>% distinct(peak, .keep_all = TRUE)
dat2.zinc <- dat.zinc %>% arrange(desc(zscore)) %>% distinct(peak, .keep_all = TRUE)
dat2.vitA <- dat.vitA %>% arrange(desc(zscore)) %>% distinct(peak, .keep_all = TRUE)
dat2.vitD <- dat.vitD %>% arrange(desc(zscore)) %>% distinct(peak, .keep_all = TRUE)
dat2.vitE <- dat.vitE %>% arrange(desc(zscore)) %>% distinct(peak, .keep_all = TRUE)
#dat2.water <- dat.water %>% arrange(desc(zscore)) %>% distinct(peak, .keep_all = TRUE)
#dat2.etOH <- dat.etOH %>% arrange(desc(zscore)) %>% distinct(peak, .keep_all = TRUE)


head(dat2.caffeine)
head(dat2.nicotine)

nrow(dat2.caffeine)
nrow(dat2.nicotine)
nrow(dat2.zinc)
nrow(dat2.vitA)
nrow(dat2.vitD)
nrow(dat2.vitE)
#nrow(dat2.water)
#nrow(dat2.etOH)

dat2 <- rbind(dat2.caffeine, dat2.nicotine, dat2.vitA, dat2.vitD, dat2.vitE, dat2.zinc)
head(dat2)
nrow(dat2)


#################################################
### resDP ###
#################################################
resDP <- read_rds("../2_Differential/3_treatEffect_2_2.outs/resDP_control.rds") %>% as.data.frame()

#resDP <- resDP %>% drop_na(p.adjusted)
##%>% drop_na(p.value)
#head(resDP)
#nrow(resDP)

#table(resDP$contrast)
#table(resDP$MCls)

sig_resDP <- resDP %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)
head(sig_resDP)
nrow(sig_resDP)

peakAll <- unique(sig_resDP$peak)

end


## #################################################
## ### resDP.treat ###
## #################################################
treats <- unique(resDP$contrast)

for (i in treats){
    df <- resDP %>% dplyr::filter(contrast==i)
    assign(paste0("resDP.",i), df)
    }

head(resDP.caffeine)



## ## #################################################
## ## ### dat.com ###
## ## #################################################
## common.peak <- intersect(dat.caffeine$peak, dat.nicotine$peak)
## length(common.peak)
## common.peak <- intersect(common.peak, dat.vitA$peak)
## length(common.peak)
## common.peak <- intersect(common.peak, dat.vitD$peak)
## length(common.peak)
## common.peak <- intersect(common.peak, dat.vitE$peak)
## length(common.peak)
## common.peak <- intersect(common.peak, dat.zinc$peak)
## length(common.peak)


## dat.com.caffeine <- dat.caffeine %>% dplyr::filter(peak %in% common.peak) %>% arrange(desc(zscore)) %>% distinct(peak, .keep_all = TRUE)
## dat.com.nicotine <- dat.nicotine %>% dplyr::filter(peak %in% common.peak) %>% arrange(desc(zscore)) %>% distinct(peak, .keep_all = TRUE)
## dat.com.zinc <- dat.zinc %>% dplyr::filter(peak %in% common.peak) %>% arrange(desc(zscore)) %>% distinct(peak, .keep_all = TRUE)
## dat.com.vitA <- dat.vitA %>% dplyr::filter(peak %in% common.peak) %>% arrange(desc(zscore)) %>% distinct(peak, .keep_all = TRUE)
## dat.com.vitD <- dat.vitD %>% dplyr::filter(peak %in% common.peak) %>% arrange(desc(zscore)) %>% distinct(peak, .keep_all = TRUE)
## dat.com.vitE <- dat.vitE %>% dplyr::filter(peak %in% common.peak) %>% arrange(desc(zscore)) %>% distinct(peak, .keep_all = TRUE)


## nrow(dat.com.caffeine)
## nrow(dat.com.nicotine)
## nrow(dat.com.vitA)
## nrow(dat.com.vitD)
## nrow(dat.com.vitE)
## nrow(dat.com.zinc)


## dat.com <- rbind(dat.com.caffeine, dat.com.nicotine, dat.com.vitA, dat.com.vitD, dat.com.vitE, dat.com.zinc)
## head(dat.com)
## nrow(dat.com)



## ## #################################################
## ## ### resDP.com ###
## ## #################################################
## head(resDP)
## head(dat.com)

## table(resDP$contrast)

## table(resDP$MCls)

## resDP.com <- resDP %>% dplyr::filter(peak %in% common.peak)

## head(resDP.com)
## nrow(resDP.com %>% dplyr::filter(contrast=="caffeine", MCls=="0-CD4Naive"))













## #################################################
## ### dat correlations ###
## #################################################
## comb <- unique(dat2$treats)
## comb

## scorr_treat <- c()
## i_treat <- c()
## j_treat <- c()

## for(i in comb){
##     for(j in comb){
##         cvt0 <- dat.com%>%dplyr::filter(treats==i)
##         cvt1 <- dat.com%>%dplyr::filter(treats==j)
##         ##
##         scorr <- as.numeric(cor.test(cvt0$zscore, cvt1$zscore, method = "spearman")$estimate)
##         print(scorr)
##         scorr_treat <- c(scorr_treat, scorr)
##         i_treat <- c(i_treat, i)
##         j_treat <- c(j_treat, j)
##         ##data1 <- plotPCA(vd1, intgroup=c("new_treat"), returnData=TRUE)
##         ##print(head(data1))
##         ##data1 <- data1 %>% mutate(sampleID=str_split(data$name, "_", simplify=T)[,4])
##         ##head(data1)
##     }
## }

## head(scorr_treat)
## head(i_treat)
## head(j_treat)

    


## peaks <- unique(dat2$peak)



## #################################################
## ### heatmap of dat2 scores ###
## #################################################
## ##----mat_estimate----
## mat <- map_dfc(comb, function(ii){
##     res2 <- dat2 %>% dplyr::filter(treats==ii)
##     b <- res2$zscore
##     names(b) <- res2$peak
##     b[sigdps.qval]
## })
## print("Number of significant in mat")
## nrow(mat.sigdps.qval.2)
## mat <- mat.sigdps.qval.2

## mat <- as.matrix(mat)
## rownames(mat) <- sigdps.qval #DP
## colnames(mat) <- comb
## head(mat)
## max(mat)
## min(mat)

## ##
## mat[is.na(mat)] = 0
## mat2 <- mat








## #################################################
## ### resDP.dat2.treat ###
## #################################################
resDP.dat2.caffeine <- left_join(resDP.caffeine, dat2.caffeine, by="peak") %>% mutate(contrast.gene=paste0(comb, ".", gene)) %>% drop_na(gene)
resDP.dat2.nicotine <- left_join(resDP.nicotine, dat2.nicotine, by="peak") %>% mutate(contrast.gene=paste0(comb, ".", gene)) %>% drop_na(gene)
resDP.dat2.zinc <- left_join(resDP.zinc, dat2.zinc, by="peak") %>% mutate(contrast.gene=paste0(comb, ".", gene)) %>% drop_na(gene)
resDP.dat2.vitA <- left_join(resDP.vitA, dat2.vitA, by="peak") %>% mutate(contrast.gene=paste0(comb, ".", gene)) %>% drop_na(gene)
resDP.dat2.vitD <- left_join(resDP.vitD, dat2.vitD, by="peak") %>% mutate(contrast.gene=paste0(comb, ".", gene)) %>% drop_na(gene)
resDP.dat2.vitE <- left_join(resDP.vitE, dat2.vitE, by="peak") %>% mutate(contrast.gene=paste0(comb, ".", gene)) %>% drop_na(gene)

head(resDP.dat2.caffeine)

sig_resDP.dat2.caffeine <- resDP.dat2.caffeine %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)
sig_resDP.dat2.nicotine <- resDP.dat2.nicotine %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)
sig_resDP.dat2.zinc <- resDP.dat2.zinc %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)
sig_resDP.dat2.vitA <- resDP.dat2.vitA %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)
sig_resDP.dat2.vitD <- resDP.dat2.vitD %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)
sig_resDP.dat2.vitE <- resDP.dat2.vitE %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)

nrow(resDP.dat2.caffeine)
nrow(resDP.dat2.nicotine)
nrow(resDP.dat2.zinc)
nrow(resDP.dat2.vitA)
nrow(resDP.dat2.vitD)
nrow(resDP.dat2.vitE)

nrow(sig_resDP.dat2.caffeine)
nrow(sig_resDP.dat2.nicotine)
nrow(sig_resDP.dat2.zinc)
nrow(sig_resDP.dat2.vitA)
nrow(sig_resDP.dat2.vitD)
nrow(sig_resDP.dat2.vitE)

resDP.dat2 <- rbind(resDP.dat2.caffeine, resDP.dat2.nicotine, resDP.dat2.vitA, resDP.dat2.vitD, resDP.dat2.vitE, resDP.dat2.zinc)
head(resDP.dat2)
nrow(resDP.dat2)

#sig_resDP.dat2 <- resDP.dat2 %>% dplyr::filter(p.adjusted < 0.1, abs(estimate)>0.5)
#nrow(sig_resDP.dat2)

sig_resDP.dat2 <- rbind(sig_resDP.dat2.caffeine, sig_resDP.dat2.nicotine, sig_resDP.dat2.vitA, sig_resDP.dat2.vitD, sig_resDP.dat2.vitE, sig_resDP.dat2.zinc)
head(sig_resDP.dat2)

nrow(sig_resDP.dat2)

head(sig_resDP.dat2)










#################################################
### resDE ###
#################################################
resDE <- read_rds("../2_Differential/3_treatEffect_2_2.outs/resDE_control.mitofilt.rds") %>% as.data.frame()

head(resDE)

#resDE <- resDE %>% drop_na(p.adjusted) %>% mutate(contrast.gene=paste0(comb, ".", gene))
#%>% drop_na(p.value)
#head(resDE)
#nrow(resDE)

table(resDE$contrast)
table(resDE$MCls)

sig_resDE <- resDE %>% dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)
head(sig_resDE)
nrow(sig_resDE)

geneAll <- unique(sig_resDE$gene)




#########################################
## heatmap ##
#########################################
head(sig_resDP.dat2.caffeine)

##--------##
comb <- unique(resDP.dat2$comb) ## this can be significant data
comb

#OR
comb.2 <- unique(resDP.dat2$contrast) ## this can be significant data
comb.2

sigdps.qval <- unique(sig_resDP.dat2$peak)

##----mat_estimate----
mat.sigdps.qval.2 <- map_dfc(comb, function(ii){
    res2 <- resDP.dat2 %>% dplyr::filter(contrast==ii)
    b <- res2$estimate
    names(b) <- res2$peak
    b[sigdps.qval]
})
print("Number of significant in mat")
nrow(mat.sigdps.qval.2)
mat <- mat.sigdps.qval.2


mat <- as.matrix(mat)
rownames(mat) <- sigdps.qval #DP
colnames(mat) <- comb.2
head(mat)
max(mat)
min(mat)

##
mat[is.na(mat)] = 0
mat2 <- mat
head(mat2)

#comb.2 <- sort(comb)
#mat3 <- mat2[, comb.2]
#head(mat3)

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
names(mybreaks) <- NULL
mybreaks

mycol <- colorRampPalette(rev(brewer.pal(n=11, name="RdBu")))(100)

col1 <- c("caffeine"="red", "nicotine"="tan",
          "vitA"="tan4", "vitD"="seagreen4",
          "vitE"="salmon3", "water"="grey", "zinc"="maroon3")

#col2 <- c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
#          "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00", "1-TCM"="pink",
#          "3-TEM"="blue", "5-CD8Naive"="green", "7-dnT"="black")

#tmp_colors <- list(treatment=col1, celltype=col2) #brewer.pal(4,"Set1")
tmp_colors <- list(treatment=col1) #brewer.pal(4,"Set1")


#x <- str_split(comb, "_", simplify=T)
#tmp_column <- data.frame(treatment=x[,1], celltype=x[,2])
tmp_column <- data.frame(treatment=comb)
rownames(tmp_column) <- comb
tmp_column

fig1 <- pheatmap(mat2,  col=mycol, breaks=mybreaks, border_color="NA",
                 cluster_rows=T, cluster_cols=F,
                 annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend=T,
                 show_colnames=T,
                 show_rownames=F,
                 treeheight_row = 0,
                 na_col="white",
                 fontsize_col=7,
                 fontsize_row=3)

figfn <- "./2.2_compareRNAandATAC_4_linkpeaks_sep.outs/Figure1.1_heatmap.beta.DP.sig_control_quantile_0.99.png"
png(figfn, width=1800, height=2100,res=225)
print(fig1)
dev.off()



end



#################################################
#----mat data for for top 10 significant only----
#################################################
nrow(sig_resDP.dat2)

comb <- unique(sig_resDP.dat2$comb)
comb
mat.10topdps.qval <- map_dfr(comb, function(ii){
    res2 <- sig_resDP.dat2 %>% dplyr::filter(comb==ii)
    res2 <- res2 %>% arrange(p.adjusted)
    res2 <- res2[1:10,]
})
print("Number of top 10 significant")
head(mat.10topdps.qval)
nrow(mat.10topdps.qval)

topdps.qval <- unique(mat.10topdps.qval$gene)
topdps.qval <- topdps.qval[!is.na(topdps.qval)]
topdps.qval

comb <- unique(resDP.dat2$comb)
comb
mat.10topdps.qval.2 <- map_dfc(comb, function(ii){
    res2 <- resDP.dat2 %>% dplyr::filter(comb==ii)
    b <- res2$statistic
    names(b) <- res2$gene
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
head(mat2)

#comb.2 <- sort(comb)
#comb.2
#mat3 <- mat2[, comb.2]


#----top10 DP plot Fig1.1----
#breaks and color scale
y <- as.numeric(mat2)
max(y)
abs(min(y))

scale <- quantile(abs(y),probs=0.99)
scale

##y0 <- y[abs(y)< quantile(abs(y),probs=0.99)] #99% percent quantile(abs(y),probs=0.99)
##mybreaks <- c(min(y),quantile(y0,probs=seq(0,1,length.out=98)),max(y))

mybreaks <- seq(-scale, scale, length.out=100)
names(mybreaks) <- NULL
mybreaks

mycol <- colorRampPalette(rev(brewer.pal(n=11, name="RdBu")))(100)

col1 <- c("caffeine"="red", "nicotine"="tan",
          "vitA"="tan4", "vitD"="seagreen4",
          "vitE"="salmon3", "water"="grey", "zinc"="maroon3")
col2 <- c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
          "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00", "1-TCM"="pink",
          "3-TEM"="blue", "5-CD8Naive"="green", "7-dnT"="black")

tmp_colors <- list(treatment=col1, celltype=col2) #brewer.pal(4,"Set1")
#tmp_colors <- list(treatment=col1) #brewer.pal(4,"Set1")


##
x <- str_split(comb, "_", simplify=T)
tmp_column <- data.frame(treatment=x[,1], celltype=x[,2])
#tmp_column <- data.frame(treatment=comb)
rownames(tmp_column) <- comb
tmp_column

fig1 <- pheatmap(mat2,  col=mycol, breaks=mybreaks, border_color="NA",
                 cluster_rows=T, cluster_cols=F, treeheight_row = 0,
                 annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend=T,
                 show_colnames=T,
                 show_rownames=T,
                 na_col="white",
                 fontsize_col=7,
                 fontsize_row=3)

figfn <- "./2.2_compareRNAandATAC_4_linkpeaks_sep.outs/Figure1.1_heatmap.beta.DP.top10_zscore_control_quantile_0.99_comb.png"
png(figfn, width=1800, height=2400,res=200)
print(fig1)
dev.off()
##
#x <- str_split(comb.2, "_", simplify=T)
# tmp_column <- data.frame(treatment=x[,1], celltype=x[,2])
# rownames(tmp_column) <- comb.2
# fig1 <- pheatmap(mat3,  col=mycol, breaks=mybreaks, border_color="NA",
#                                   cluster_rows=T, cluster_cols=F,treeheight_row = 0,
#                                   annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend=T,
#                                   show_colnames=T,
#                                   show_rownames=T,
#                                   na_col="white",
#                                   fontsize_col=7,
#                                   fontsize_row=2)
# figfn <- "./3_treatEffect_2_2.outs/Figure1.1_heatmap.beta.DP.top10.treat_control_quantile_0.99.png"
# png(figfn, width=1800, height=2100,res=225)
# print(fig1)
# dev.off()





#################################################
### CORRELATION
#################################################
 MCls <- unique(resDE$MCls)
MCls
comb.2 <- MCls

treats <- unique(resDE$contrast)
#treats <- c("caffeine", "nicotine", "vitA", "vitD", "vitE", "zinc")
treats


#----all----
logFC <- "NA"
include <- "all"


#----Significance.0----
logFC <- 0
conditions <- c("both", "sc", "GXE")
include <- "union"


## alt_func <- function(j){
##         if (j=="caffeine"){
##             #scorr_caff <- c(scorr_caff, "NA")
##             scorr_caff_z <- c(scorr_caff_z, 0)
##         } else if (j=="nicotine"){
##             #scorr_nic <- c(scorr_nic, "NA")
##             scorr_nic_z <- c(scorr_nic_z, 0)
##         } else if (j=="vitA"){
##             #scorr_vitA <- c(scorr_vitA, "NA")
##             scorr_vitA_z <- c(scorr_vitA_z, 0)
##         } else if (j=="vitD"){
##             #scorr_vitD <- c(scorr_vitD, "NA")
##             scorr_vitD_z <- c(scorr_vitD_z, 0)
##         } else if (j=="vitE"){
##             #scorr_vitE <- c(scorr_vitE, "NA")
##             scorr_vitE_z <- c(scorr_vitE_z, 0)
##         } else if (j=="water"){
##             #scorr_water <- c(scorr_water, "NA")
##             scorr_water_z <- c(scorr_water, 0)
##         } else if (j=="zinc"){
##             #scorr_zinc <- c(scorr_zinc, "NA")
##             scorr_zinc_z <- c(scorr_zinc_z, 0)
##         }
## }


scorr_caff <- c()
scorr_nic <- c()
scorr_vitA <- c()
scorr_vitD <- c()
scorr_vitE <- c()
scorr_zinc <- c()
scorr_caff_z <- c()
scorr_nic_z <- c()
scorr_vitA_z <- c()
scorr_vitD_z <- c()
scorr_vitE_z <- c()
scorr_water_z <- c()
scorr_zinc_z <- c()

head(resDP.link)

head(resDP.peakAnno)


for (j in treats){
    resDE2 <- resDE %>% dplyr::filter(contrast==j)
    ##
    resDP2 <- resDP.peakAnno %>% dplyr::filter(contrast==j)
    ##
    ##resDP2 <- resDP.link %>% dplyr::filter(contrast==j)
    ##
    resDP2.sum <- resDP2 %>% arrange(p.adjusted)
    resDP2.sum <- resDP2.sum %>% distinct(symbol, .keep_all = TRUE) #for peakAnno
    ##resDP2.sum <- resDP2.sum %>% distinct(gene, .keep_all = TRUE) #for link
    ##
    resDE.DP <- left_join(resDE2, resDP2.sum, by="contrast.gene")
    resDE.DP <- resDE.DP %>% dplyr::filter(!is.na(estimate.y))
    ##
    resDE.DP <- resDE.DP %>% mutate(Significance=ifelse(Significance.y %in% c("Up", "Down") &  Significance.x %in% c("Up", "Down"), "both", ifelse(Significance.y %in% c("Up", "Down") & (Significance.x %in% c("Not Significant", "NA") | is.na(Significance.x)), "DP", ifelse((Significance.y %in% c("Not Significant", "NA") | is.na(Significance.y)) & Significance.x %in% c("Up", "Down"), "DG", "NA"))))
    ##
    resDE.DP <- resDE.DP %>% mutate(Significance.25=ifelse(Significance.25.y %in% c("Up", "Down") &  Significance.25.x %in% c("Up", "Down"), "both", ifelse(Significance.25.y %in% c("Up", "Down") & (Significance.25.x %in% c("Not Significant", "NA") | is.na(Significance.25.x)), "DP", ifelse((Significance.25.y %in% c("Not Significant", "NA") | is.na(Significance.25.y)) & Significance.25.x %in% c("Up", "Down"), "DG", "NA"))))
    ##
    resDE.DP <- resDE.DP %>% mutate(Significance.0=ifelse(Significance.0.y %in% c("Up", "Down") &  Significance.0.x %in% c("Up", "Down"), "both", ifelse(Significance.0.y %in% c("Up", "Down") & (Significance.0.x %in% c("Not Significant", "NA") | is.na(Significance.0.x)), "DP", ifelse((Significance.0.y %in% c("Not Significant", "NA") | is.na(Significance.0.y)) & Significance.0.x %in% c("Up", "Down"), "DG", "NA"))))
    ####################
    ## scorr_values ##
    ####################
    for (i in comb.2){
        tryCatch({
            res1 <- resDE.DP %>% dplyr::filter(contrast.x==j, MCls.x==i)       #                          %>% dplyr::filter(Significance.0 %in% conditions)
            ##scorr <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$estimate)
            scorr_z <- as.numeric(cor.test(res1$statistic.x, res1$statistic.y, method = "spearman")$estimate)
            ##print(scorr)
            print(scorr_z)
            if (is.numeric(scorr_z) & j=="caffeine"){
                ##scorr_caff <- c(scorr_caff, scorr)
                scorr_caff_z <- c(scorr_caff_z, scorr_z)
            } else if (is.numeric(scorr_z) & j=="nicotine"){
                ##scorr_nic <- c(scorr_nic, scorr)
                scorr_nic_z <- c(scorr_nic_z, scorr_z)
            } else if (is.numeric(scorr_z) & j=="vitA"){
                ##scorr_vitA <- c(scorr_vitA, scorr)
                scorr_vitA_z <- c(scorr_vitA_z, scorr_z)
            } else if (is.numeric(scorr_z) & j=="vitD"){
                ##scorr_vitD <- c(scorr_vitD, scorr)
                scorr_vitD_z <- c(scorr_vitD_z, scorr_z)
            } else if (is.numeric(scorr_z) & j=="vitE"){
                ##scorr_vitE <- c(scorr_vitE, scorr)
                scorr_vitE_z <- c(scorr_vitE_z, scorr_z)
            } else if (is.numeric(scorr_z) & j=="water"){
                ##scorr_water <- c(scorr_water, scorr)
                scorr_water_z <- c(scorr_water, scorr_z)
            } else if (is.numeric(scorr_z) & j=="zinc"){
                ##scorr_zinc <- c(scorr_zinc, scorr)
                scorr_zinc_z <- c(scorr_zinc_z, scorr_z)
            }
            ##
        },  error=function(j){
            if (j=="caffeine"){
                ##scorr_caff <- c(scorr_caff, "NA")
                scorr_caff_z <- c(scorr_caff_z, 0)
        } else if (j=="nicotine"){
            #scorr_nic <- c(scorr_nic, "NA")
            scorr_nic_z <- c(scorr_nic_z, 0)
        } else if (j=="vitA"){
            #scorr_vitA <- c(scorr_vitA, "NA")
            scorr_vitA_z <- c(scorr_vitA_z, 0)
        } else if (j=="vitD"){
            #scorr_vitD <- c(scorr_vitD, "NA")
            scorr_vitD_z <- c(scorr_vitD_z, 0)
        } else if (j=="vitE"){
            #scorr_vitE <- c(scorr_vitE, "NA")
            scorr_vitE_z <- c(scorr_vitE_z, 0)
        } else if (j=="water"){
            #scorr_water <- c(scorr_water, "NA")
            scorr_water_z <- c(scorr_water, 0)
        } else if (j=="zinc"){
            #scorr_zinc <- c(scorr_zinc, "NA")
            scorr_zinc_z <- c(scorr_zinc_z, 0)
        }
        })    
    }
}


#cat("ERROR", "\n")

scorr_caff_z <- c(scorr_caff_z, 0)
scorr_nic_z <- c(scorr_nic_z, 0)
scorr_zinc_z <- c(scorr_zinc_z, 0)
scorr_vitD_z <- c(scorr_vitD_z, 0)


scorr_caff_z
scorr_nic_z
scorr_vitA_z
scorr_vitD_z
scorr_vitE_z
scorr_water_z
scorr_zinc_z






library(ggsci)
library(circlize)
library(reshape2)
library("ggpubr")


scorr_data <- data.frame(MCls, scorr_caff_z, scorr_nic_z,
                         scorr_vitA_z, scorr_vitD_z,
                         scorr_vitE_z, scorr_zinc_z)
is.num <- sapply(scorr_data, is.numeric)
scorr_data[is.num] <- lapply(scorr_data[is.num], round, 2)
##
sig.melt <- melt(scorr_data)
##
sig.p <- ggplot(sig.melt, aes(x = variable, y = MCls)) +
    geom_tile(aes(fill=value), color="black")+
    scale_fill_gradient(low = "white", high = "white")+
    geom_text(aes(label = value), color = "black", size = 5) +
    ggtitle(paste0("logFC_", logFC, "_", include, "_zscorevszscore_corr_DG_vs_DP")) +
    theme(axis.text.x = element_text(angle = 90),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position="none")
##
figure <- ggarrange(sig.p + font("x.text", size = 14))
png(paste0("./2.2_compareRNAandATAC_3.outs/logFC_", logFC, "_", include, "_zscorevszscore_corr_link.png"), width=500, height=500, pointsize=16, res=125)
print(figure)
dev.off()


end

## checking if gene_name and symbol are same
for (i in (sig_resDP.peakAnno$gene_name==sig_resDP.peakAnno$symbol)){
    if(i==FALSE){
    print(sig_resDP.peakAnno$gene_name)
    print(sig_resDP.peakAnno$symbol)
    } 
}








#################################################
### SCATTER PLOT 
#################################################
library(ggplot2)
library(ggpmisc)

comb <- unique(sig_resDP.dat2$comb)
comb

#comb.2 <- unique(sig_resDP.dat2$contrast)
#comb.2

#comb.3 <- unique(sig_resDP.dat2$MCls)
#comb.3

head(sig_resDP.dat2)

sig_comb <- map_dfr(comb, function(ii){
    resDP2 <- sig_resDP.dat2%>%dplyr::filter(comb==ii) #link
    resDE2 <- sig_resDE%>%dplyr::filter(comb==ii)
    ##DP <- unique(resDP2$peak)
    DEG <- unique(resDE2$gene)
    ##x <- peakAnno2%>%dplyr::filter(peak%in%DP, gene_name%in%DEG)
    x <- resDP2 %>% dplyr::filter(gene%in%DEG)
    if(nrow(x)==0){
        res1 <- NULL
        ##next
    }else{
        ##res1 <- unique(x$gene_name)
        ##res1 <- x
        #res1 <- x%>%mutate(gene=gene_name)
        res1 <- x
        res1 <- res1%>%left_join(resDE2, by="gene")
    }
    ##res <- list(DEG=res1)
    res1
})

## columns <- colnames(sig_resDP.dat2)
## columns

## res <- data.frame(matrix(nrow = 0, ncol = length(columns)))
## res

## colnames(res) <- cols
## res

## head(res)

## for (i in comb){
##     for (j in MCls){
##     resDP2 <- sig_resDP.dat2%>%dplyr::filter(contrast==i, MCls==j) #link
##     resDE2 <- sig_resDE%>%dplyr::filter(contrast==i, MCls==j)
##     ##DP <- unique(resDP2$peak)
##     DEG <- unique(resDE2$gene)
##     ##x <- peakAnno2%>%dplyr::filter(peak%in%DP, gene_name%in%DEG)
##     x <- resDP2 %>% dplyr::filter(gene%in%DEG)
##     if(nrow(x)==0){
##         res1 <- NULL
##         ##next
##     }else{
##         ##res1 <- unique(x$gene_name)
##         ##res1 <- x
##         #res1 <- x%>%mutate(gene=gene_name)
##         res1 <- x
##         res1 <- res1%>%left_join(resDE2, by="gene")
##     }
##     ##res <- list(DEG=res1)
##     ##res1
##     ##assign(paste0("pbmc.",i), df2)
##     res <- rbind(res, res1)
## }}


## head(res)

head(sig_comb)

nrow(sig_comb)

colnames(sig_comb)
length(unique(sig_comb$gene))
table(sig_comb$comb.x)

max(sig_comb$estimate.y)
min(sig_comb$estimate.y)
max(sig_comb$estimate.x)
min(sig_comb$estimate.x)






##----plot loop----
#install.packages("ggpmisc")
library(ggpmisc)

##----celltype----
comb<- unique(resDP.dat2$contrast)
comb



head(sig_comb)

## figfn <- paste0("./2.2_compareRNAandATAC_3.outs/Figure1.1_significant_scatter_macs2_0.1_cn_peakAnno_control_celltype.png")
## png(figfn, width=6000, height=8000, pointsize=33, res=250)
## par(mar=c(4,4,2,2),mgp=c(2,1,0))
## x <- matrix(1:8, 2, 4, byrow=T)
## layout(x)
## ##for (i in c(1)){
res1 <- sig_comb %>% dplyr::filter(MCls.x==comb[1])
p1 <- ggplot(res1, aes(x=estimate.y, y=estimate.x, color=MCls.x), show.legend = FALSE) +
    geom_point()+
    xlab("sig_genes") + ylab("sig_peaks")+
    xlim(-4, 10) + ylim(-4, 5)+
    ##stat_poly_line() +
    ##stat_poly_eq() +
    ggtitle(comb.3[1])
p1 + theme(legend.position = "none")
scorr_est_1 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$estimate)
scorr_pval_1 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$p.value)
##print(p)
##
res1 <- sig_comb %>% dplyr::filter(MCls.x==comb[2])
p2 <- ggplot(res1, aes(x=estimate.y, y=estimate.x, color=MCls.x), show.legend = FALSE) +
    geom_point()+
    xlab("sig_genes") + ylab("sig_peaks")+
    xlim(-4, 10) + ylim(-4, 5)+
    ##stat_poly_line() +
    ##stat_poly_eq() +
    ggtitle(comb.3[2])
p2 + theme(legend.position = "none")
scorr_est_2 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$estimate)
scorr_pval_2 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$p.value)
##print(p)
##
res1 <- sig_comb %>% dplyr::filter(MCls.x==comb[3])
p3 <- ggplot(res1, aes(x=estimate.y, y=estimate.x, color=MCls.x), show.legend = FALSE) +
    geom_point()+
    xlab("sig_genes") + ylab("sig_peaks")+
    xlim(-4, 10) + ylim(-4, 5)+
    ##stat_poly_line() +
    ##stat_poly_eq() +
    ggtitle(comb.3[3])
p3 + theme(legend.position = "none")
scorr_est_3 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$estimate)
scorr_pval_3 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$p.value)
##print(p)
##
res1 <- sig_comb %>% dplyr::filter(MCls.x==comb[4])
p4 <- ggplot(res1, aes(x=estimate.y, y=estimate.x, color=MCls.x), show.legend = FALSE) +
    geom_point()+
    xlab("sig_genes") + ylab("sig_peaks")+
    xlim(-4, 10) + ylim(-4, 5)+
    ##stat_poly_line() +
    ##stat_poly_eq() +
    ggtitle(comb.3[4])
p4 + theme(legend.position = "none")
scorr_est_4 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$estimate)
scorr_pval_4 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$p.value)
##print(p)
##
res1 <- sig_comb %>% dplyr::filter(MCls.x==comb[5])
p5 <- ggplot(res1, aes(x=estimate.y, y=estimate.x, color=MCls.x), show.legend = FALSE) +
    geom_point()+
    xlab("sig_genes") + ylab("sig_peaks")+
    xlim(-4, 10) + ylim(-4, 5)+
    ##stat_poly_line() +
    ##stat_poly_eq() +
    ggtitle(comb.3[5])
p5 + theme(legend.position = "none")
scorr_est_5 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$estimate)
scorr_pval_5 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$p.value)
##print(p)
##
res1 <- sig_comb %>% dplyr::filter(MCls.x==comb[6])
p6 <- ggplot(res1, aes(x=estimate.y, y=estimate.x, color=MCls.x), show.legend = FALSE) +
    geom_point()+
    xlab("sig_genes") + ylab("sig_peaks")+
    xlim(-4, 10) + ylim(-4, 5)+
    ##stat_poly_line() +
    ##stat_poly_eq() +
    ggtitle(comb.3[6])
p6 + theme(legend.position = "none")
scorr_est_6 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$estimate)
scorr_pval_6 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$p.value)
##print(p)


## ##
## res1 <- sig_comb %>% dplyr::filter(MCls.x==comb[7])
## p7 <- ggplot(res1, aes(x=estimate.y, y=estimate.x, color=MCls.x), show.legend = FALSE) +
##     geom_point()+
##     xlab("sig_genes") + ylab("sig_peaks")+
##     xlim(-4, 10) + ylim(-4, 5)+
##     ##stat_poly_line() +
##     ##stat_poly_eq() +
##     ggtitle(comb.3[7])
## p7 + theme(legend.position = "none")
## scorr_est_7 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$estimate)
## scorr_pval_7 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$p.value)
## ##print(p)
## ##
## res1 <- sig_comb %>% dplyr::filter(MCls.x==comb[8])
## p8 <- ggplot(res1, aes(x=estimate.y, y=estimate.x, color=MCls.x), show.legend = FALSE) +
##     geom_point()+
##     xlab("sig_genes") + ylab("sig_peaks")+
##     xlim(-4, 10) + ylim(-4, 5)+
##     ##stat_poly_line() +
##     ##stat_poly_eq() +
##     ggtitle(comb.3[8])
## p8 + theme(legend.position = "none")
## scorr_est_8 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$estimate)
## scorr_pval_8 <- as.numeric(cor.test(res1$estimate.x, res1$statistic.y, method = "spearman")$p.value)
## ##print(p)
## ##
## ##print(mtext(1, side=4, line=0.5, cex=1, col="blue"))
## ##}
## ##dev.off()



print(scorr_est_1)
print(scorr_est_2)
print(scorr_est_3)
print(scorr_est_4)
print(scorr_est_5)
print(scorr_est_6)
print(scorr_est_7)
print(scorr_est_8)
print(scorr_pval_1)
print(scorr_pval_2)
print(scorr_pval_3)
print(scorr_pval_4)
print(scorr_pval_5)
print(scorr_pval_6)
print(scorr_pval_7)
print(scorr_pval_8)



library("ggpubr")
figure <- ggarrange(p1, p2, p3, p4, p5, p6, p7, p8 +
                                             font("x.text", size = 18),
                    ncol = 4, nrow = 2)
png(paste0("./2.2_compareRNAandATAC_3.outs/Figure1.1_significant_scatter_macs2_0.1_cn_peakAnno_control_celltype_noline.png"), width=4000, height=2000, pointsize=20, res=175)
print(figure)
dev.off()




##----treatment----
unique(sig_comb$contrast.x)


## figfn <- paste0("./2.2_compareRNAandATAC_3.outs/Figure1.1_significant_scatter_macs2_0.1_cn_peakAnno_control_treat.png")
## png(figfn, width=6000, height=8000, pointsize=33, res=250)
## par(mar=c(4,4,2,2),mgp=c(2,1,0))
## x <- matrix(1:5, 1, 5, byrow=T)
## layout(x)
##print(layout(x))
##for (i in comb.2){
res1 <- sig_comb %>% dplyr::filter(contrast.x=="caffeine")
p1 <- ggplot(res1, aes(x=estimate.y, y=estimate.x, color=contrast.x), show.legend = FALSE) +
    geom_point()+
    xlab("sig_genes") + ylab("sig_peaks")+
    xlim(-4, 10) + ylim(-4, 5)+
    ##stat_poly_line() +
    ##stat_poly_eq() +
    ggtitle("caffeine")
p1 + theme(legend.position = "none")
p1 + guides(color = FALSE)
scorr_est_1 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$estimate)
scorr_pval_1 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$p.value)
##print(p)
##
res1 <- sig_comb %>% dplyr::filter(contrast.x=="nicotine")
p2 <- ggplot(res1, aes(x=estimate.y, y=estimate.x, color=contrast.x), show.legend = FALSE) +
    geom_point()+
    xlab("sig_genes") + ylab("sig_peaks")+
    xlim(-4, 10) + ylim(-4, 5)+
    ##stat_poly_line() +
    ##stat_poly_eq() +
    ggtitle("nicotine")
p2 + theme(legend.position = "none")
p2 + guides(color = FALSE)
scorr_est_2 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$estimate)
scorr_pval_2 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$p.value)
##print(p)
##
res1 <- sig_comb %>% dplyr::filter(contrast.x=="vitA")
p3 <- ggplot(res1, aes(x=estimate.y, y=estimate.x, color=contrast.x), show.legend = FALSE) +
    geom_point()+
    xlab("sig_genes") + ylab("sig_peaks")+
    xlim(-4, 10) + ylim(-4, 5)+
    ##stat_poly_line() +
    ##stat_poly_eq() +
    ggtitle("vitA")
p3 + theme(legend.position = "none")
p3 + guides(color = FALSE)
scorr_est_3 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$estimate)
scorr_pval_3 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$p.value)
##print(p)
##
res1 <- sig_comb %>% dplyr::filter(contrast.x=="vitD")
p4 <- ggplot(res1, aes(x=estimate.y, y=estimate.x, color=contrast.x), show.legend = FALSE) +
    geom_point()+
    xlab("sig_genes") + ylab("sig_peaks")+
    xlim(-4, 10) + ylim(-4, 5)+
    ##stat_poly_line() +
    ##stat_poly_eq() +
    ggtitle("vitD")
p4 + theme(legend.position = "none")
p4 + guides(color = FALSE)
scorr_est_4 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$estimate)
scorr_pval_4 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$p.value)
##print(p)
##
## res1 <- sig_comb %>% dplyr::filter(contrast.x=="vitE")
## p5 <- ggplot(res1, aes(x=estimate.y, y=estimate.x, color=contrast.x), show.legend = FALSE) +
##     geom_point()+
##     xlab("sig_genes") + ylab("sig_peaks")+
##     xlim(-4, 10) + ylim(-4, 5)+
##     ##stat_poly_line() +
##     ##stat_poly_eq() +
##     ggtitle("vitE")
## p5 + theme(legend.position = "none")
## p5 + guides(color = FALSE)
## scorr_est_5 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$estimate)
## scorr_pval_5 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$p.value)
##
res1 <- sig_comb %>% dplyr::filter(contrast.x=="zinc")
p6 <- ggplot(res1, aes(x=estimate.y, y=estimate.x, color=contrast.x), show.legend = FALSE) +
    geom_point()+
    xlab("sig_genes") + ylab("sig_peaks")+
    xlim(-4, 10) + ylim(-4, 5)+
    ##stat_poly_line() +
    ##stat_poly_eq() +
    ggtitle("zinc")
p6 + theme(legend.position = "none")
p6 + guides(color = FALSE)
scorr_est_6 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$estimate)
scorr_pval_6 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$p.value)
##print(p)

##
##print(mtext(1, side=4, line=0.5, cex=1, col="blue"))
#}
##dev.off()

nrow(sig_comb %>% dplyr::filter(contrast.x=="vitE"))

print(scorr_est_1)
print(scorr_est_2)
print(scorr_est_3)
print(scorr_est_4)
#print(scorr_est_5)
print(scorr_est_6)
print(scorr_pval_1)
print(scorr_pval_2)
print(scorr_pval_3)
print(scorr_pval_4)
#print(scorr_pval_5)
print(scorr_pval_6)


library("ggpubr")

figure <- ggarrange(p1, p2, p3, p4, p6 +
                                    font("x.text", size = 16),
                    ncol = 3, nrow = 2)

png(paste0("./2.2_compareRNAandATAC_4_linkpeaks_sep.outs/Figure1.1_significant_scatter_macs2_0.1_cn_linkpeaks_control_treat_noline.png"), width=3000, height=2000, pointsize=18, res=175)
print(figure)
dev.off()




##----celltype and treatments----
figfn <- paste0("./2.2_compareRNAandATAC_3.outs/Figure1.1_significant_scatter_macs2_0.1_cn_peakAnno_control_allcomb.png")
png(figfn, width=6000, height=8000, pointsize=33, res=250)
par(mar=c(4,4,2,2),mgp=c(2,1,0))
x <- matrix(1:40, 8, 5, byrow=T)
layout(x)
for (i in comb.3){
    res1 <- sig_comb %>% dplyr::filter(MCls.x==i, contrast.x=="caffeine")
    p <- ggplot(res1, aes(x=estimate.y, y=estimate.x, color=contrast.x)) +
        geom_point()+
        xlab("sig_genes") + ylab("sig_peaks")+
        xlim(-4, 10) + ylim(-4, 5)+
        stat_poly_line() +
        stat_poly_eq() +
        ggtitle("fold.enrichment")
    print(p)
    ##
    res1 <- sig_comb %>% dplyr::filter(MCls.x==i, contrast.x=="nicotine")
    p <- ggplot(res1, aes(x=estimate.y, y=estimate.x, color=contrast.x)) +
        geom_point()+
        xlab("sig_genes") + ylab("sig_peaks")+
        xlim(-4, 10) + ylim(-4, 5)+
        stat_poly_line() +
        stat_poly_eq() +
        ggtitle("fold.enrichment")
    print(p)
    ##
    res1 <- sig_comb %>% dplyr::filter(MCls.x==i, contrast.x=="vitA")
    p <- ggplot(res1, aes(x=estimate.y, y=estimate.x, color=contrast.x)) +
        geom_point()+
        xlab("sig_genes") + ylab("sig_peaks")+
        xlim(-4, 10) + ylim(-4, 5)+
        stat_poly_line() +
        stat_poly_eq() +
        ggtitle("fold.enrichment")
    print(p)
    ##
    res1 <- sig_comb %>% dplyr::filter(MCls.x==i, contrast.x=="vitD")
    p <- ggplot(res1, aes(x=estimate.y, y=estimate.x, color=contrast.x)) +
        geom_point()+
        xlab("sig_genes") + ylab("sig_peaks")+
        xlim(-4, 10) + ylim(-4, 5)+
        stat_poly_line() +
        stat_poly_eq() +
        ggtitle("fold.enrichment")
    print(p)
    ##
    res1 <- sig_comb %>% dplyr::filter(MCls.x==i, contrast.x=="zinc")
    p <- ggplot(res1, aes(x=estimate.y, y=estimate.x, color=contrast.x)) +
        geom_point()+
        xlab("sig_genes") + ylab("sig_peaks")+
        xlim(-4, 10) + ylim(-4, 5)+
        stat_poly_line() +
        stat_poly_eq() +
        ggtitle("fold.enrichment")
    print(p)
    ##
    print(mtext(i, side=4, line=0.5, cex=1, col="blue"))
}
dev.off()
    


end


## ##----plot----
## lm_eqn <- function(df){
##     m <- lm(y ~ x, df);
##     eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
##                      list(a = format(unname(coef(m)[1]), digits = 2),
##                           b = format(unname(coef(m)[2]), digits = 2),
##                           r2 = format(summary(m)$r.squared, digits = 3)))
##     as.character(as.expression(eq));
## }


## png("./2.2_compareRNAandATAC_3.outs/Figure1.1_significant_scatter_macs2_0.1_cn_peakAnno_control_scatter.png")

## library("ggpubr")
## p <- ggscatter(sig_comb, x = "estimate.y", y = "estimate.x",
##           #color=sig_comb$comb.x,
##           add = "reg.line", conf.int = TRUE,
##           cor.coef = TRUE, cor.method = "spearman",
##           xlab = "sig_genes", ylab = "sig_peaks")


## p <- ggplot(sig_comb, aes(x=estimate.y, y=estimate.x, color=comb.x)) +
##     geom_point()+
##     xlab("sig_genes") + ylab("sig_peaks")+
##     xlim(-4, 10) + ylim(-4, 5)+
##     ##stat_poly_line() +
##     ##stat_poly_eq() +
##     ##geom_text(label = lm_eqn(df), parse = TRUE)+
##     ggtitle("fold.enrichment")

## print(p)
## dev.off()

## scorr <- as.numeric(cor.test(sig_comb$estimate.y, sig_comb$estimate.x, method = "spearman")$estimate)
## print(scorr)

## scorr.pval <- as.numeric(cor.test(sig_comb$estimate.y, sig_comb$estimate.x, method = "spearman")$p.value)
## print(scorr.pval)






#################################################
### QQ PLOT 
#################################################
head(resDE)
nrow(resDE)

head(resDP.peakAnno)
nrow(resDP.peakAnno)

head(resDP.pa.sum)
nrow(resDP.pa.sum)

comb <- sort(unique(resDE$comb))
comb

dfNew <- map_dfr(comb, function(ii){
  res2 <- sig_resDP.peakAnno%>%dplyr::filter(comb==ii)
  DEG <- sig_resDE%>%dplyr::filter(comb==ii)%>%dplyr::pull(gene)
  res2 <- res2%>%mutate(is_DEG=ifelse(symbol%in%DEG,1,0))
  ###
  dx <- map_dfr(c(0,1),function(i){
     di <- res2%>%dplyr::filter(is_DEG==i)%>%arrange(p.value)
     ngene <- nrow(di)
     di <- di%>%mutate(observed=-log10(p.value), expected=-log10(ppoints(ngene)))
     di
  })
  dx
})

head(dfNew)
nrow(dfNew)

###
###
## lab1 <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX", "PHA"="PHA", "PHA-DEX"="PHA+DEX")

## p1 <- ggplot(dfNew, aes(x=expected,y=observed, colour=factor(is_DEG)))+
##    ggrastr::rasterise(geom_point(size=0.3),dpi=300)+
##    geom_abline(colour="red")+
##    scale_colour_manual(values=c("0"="grey40","1"="green"),
##       labels=c("0"="Not DEG", "1"="DEG"),
##       guide=guide_legend(override.aes=list(size=1)))+
##    facet_grid(MCls~contrast, scales="free", labeller=labeller(contrast=lab1))+
##    xlab(bquote("Expected"~ -log[10]~"("~italic(plain(P))~")"))+
##    ylab(bquote("Observed"~-log[10]~"("~italic(plain(P))~")"))+
##    theme_bw()+
##    theme(legend.title=element_blank(),strip.text=element_text(size=12))

## figfn <- "./2.2_compareRNAandATAC.outs/Figure1.1_qq.signac.png"
## png(figfn, width=750, height=750, res=120)
## print(p1)
## dev.off()


lab1 <- c("caffeine"="caffeine", "nicotine"="nicotine", "zinc"="zinc",
          "vitA"="vitA", "vitD"="vitD", "vitE"="vitE")
          ##"water"="water")

lab2 <- c("4-Bcell"="4-Bcell", "6-Monocyte"="6-Monocyte",
          "2-NKcell"="2-NKcell", "0-CD4Naive"="0-CD4Naive", "1-TCM"="1-TCM",
          "3-TEM"="3-TEM", "5-CD8Naive"="5-CD8Naive", "7-dnT"="7-dnT")

p <- ggplot(dfNew, aes(x=expected,y=observed, colour=factor(is_DEG)))+
    ggrastr::rasterise(geom_point(size=0.3),dpi=300)+
    geom_abline(colour="red")+
    scale_colour_manual(values=c("1"="green", "0"="grey40"),
                        labels=c("1"="sig DAR, sig DEG", "0"="sig DAR, non-sig DEG"),
                        guide=guide_legend(override.aes=list(size=1)))+
    facet_grid(MCls~contrast, scales="free",
               labeller=labeller(contrast=lab1, MCls=lab2))+
    xlab(bquote("Expected"~ -log[10]~"("~italic(plain(P))~")"))+
    ylab(bquote("Observed"~-log[10]~"("~italic(plain(P))~")"))+
    theme_bw()+
    theme(strip.text=element_text(size=12))

figfn <- "./2.2_compareRNAandATAC_3.outs/Figure1.1_qq.signac_control.png"
png(figfn, width=1600, height=1600, res=175)
print(p)
dev.off()



## #######################################
## ### peak annotation from ChIPseeker ###
## #######################################
## peakAnno <- read_rds("./2.2_compareRNAandATAC.outs/2_annot.ChIPseeker.rds")

## ##
## figfn <- "./2.2_compareRNAandATAC.outs/Figure2.1_annot.pie.png"
## png(figfn, width=480, height=320)
## plotAnnoPie(peakAnno)
## dev.off()

## ###
## figfn <- "./2.2_compareRNAandATAC.outs/Figure2.2_annot.ven.png"
## png(figfn, width=700, height=500)
## vennpie(peakAnno)
## dev.off()

## ###
## peakAnno <- as.data.frame(peakAnno)
## peakAnno <- peakAnno%>%
##    mutate(chr=gsub("chr","",seqnames),
##           peak_region=paste(chr,start,end,sep="-"))%>%
##    dplyr::select(peak_region, geneId, flank_geneIds) 


## ###
## res <- read_rds("./1.2_DiffPeak.outs/2.0_DESeq.results.rds")%>%as.data.frame()
## res <- res%>%dplyr::rename("peak_region"="gene")%>%drop_na(p.value)

## res <- res%>%
##    left_join(peakAnno, by="peak_region")%>%
##    mutate(comb=paste(MCls, contrast, sep="_"))


## ### previous identified DEGs
## fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
## resDE <- read_rds(fn)%>%dplyr::filter(qval<0.1, abs(beta)>0.5)
## resDE <- resDE%>%
##    mutate(comb=paste(MCls, contrast, sep="_"))


## ###
## ### qq for cloest genes
## comb <- sort(unique(resDE$comb))
## dfNew <- map_dfr(comb, function(ii){
##   res2 <- res%>%dplyr::filter(comb==ii)
##   DEG <- resDE%>%dplyr::filter(comb==ii)%>%dplyr::pull(gene)
##   res2 <- res2%>%mutate(is_DEG=ifelse(geneId%in%DEG,1,0))  
##   ###
##   dx <- map_dfr(c(0,1),function(i){
##      di <- res2%>%dplyr::filter(is_DEG==i)%>%arrange(p.value)
##      ngene <- nrow(di)
##      di <- di%>%mutate(observed=-log10(p.value), expected=-log10(ppoints(ngene)))
##      di
##   })
##   dx
## })

## ###
## ###
## lab1 <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX", "PHA"="PHA", "PHA-DEX"="PHA+DEX")
## p1 <- ggplot(dfNew, aes(x=expected,y=observed, colour=factor(is_DEG)))+
##    ggrastr::rasterise(geom_point(size=0.3),dpi=300)+
##    geom_abline(colour="red")+
##    scale_colour_manual(values=c("0"="grey40","1"="green"),
##       labels=c("0"="Not DEG", "1"="DEG"),
##       guide=guide_legend(override.aes=list(size=1)))+
##    facet_grid(MCls~contrast, scales="free", labeller=labeller(contrast=lab1))+
##    xlab(bquote("Expected"~ -log[10]~"("~italic(plain(P))~")"))+
##    ylab(bquote("Observed"~-log[10]~"("~italic(plain(P))~")"))+
##    theme_bw()+
##    theme(legend.title=element_blank(),strip.text=element_text(size=12))

## figfn <- "./2.2_compareRNAandATAC.outs/Figure2.3_qq.DEG.cloest.png"
## png(figfn, width=750, height=750, res=120)
## print(p1)
## dev.off()



## ###
## ### flanking region
## ### generate data for qq plot
## comb <- sort(unique(resDE$comb))
## dfNew <- map_dfr(comb, function(ii){
##   res2 <- res%>%dplyr::filter(comb==ii)
##   DEG <- resDE%>%dplyr::filter(comb==ii)%>%dplyr::pull(gene)
##   geneList <- str_split(res2$flank_geneIds, ";")
##   is_DEG <- map_dbl(geneList, function(i){
##      x <- ifelse(any(i%in%DEG),1,0)
##      x
##   })
##   res2$is_DEG <- is_DEG
  
##   ###
##   dx <- map_dfr(c(0,1),function(i){
##      di <- res2%>%dplyr::filter(is_DEG==i)%>%arrange(p.value)
##      ngene <- nrow(di)
##      di <- di%>%mutate(observed=-log10(p.value), expected=-log10(ppoints(ngene)))
##      di
##   })
##   dx
## })


## ###
## lab1 <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX", "PHA"="PHA", "PHA-DEX"="PHA+DEX")
## p2 <- ggplot(dfNew, aes(x=expected,y=observed, colour=factor(is_DEG)))+
##    ggrastr::rasterise(geom_point(size=0.3),dpi=300)+
##    geom_abline(colour="red")+
##    scale_colour_manual(values=c("0"="grey40","1"="green"),
##       labels=c("0"="Not DEG", "1"="DEG"),
##       guide=guide_legend(override.aes=list(size=1)))+
##    facet_grid(MCls~contrast, scales="free", labeller=labeller(contrast=lab1))+
##    xlab(bquote("Expected"~ -log[10]~"("~italic(plain(P))~")"))+
##    ylab(bquote("Observed"~-log[10]~"("~italic(plain(P))~")"))+
##    theme_bw()+
##    theme(legend.title=element_blank(),strip.text=element_text(size=12))

## figfn <- "./2.2_compareRNAandATAC.outs/Figure2.4_qq.DEG.png"
## png(figfn, width=750, height=750, res=120)
## print(p2)
## dev.off()




## ############################
## ### if enrich are in DVG ###
## ############################
## ### previous identified DEGs
## fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/10_RNA.Variance_output/tmp9/3_phiNew.meta"
## resDE <- read.table(fn,header=T)%>%dplyr::filter(qval<0.1, abs(beta)>0.5)
## resDE <- resDE%>%
##    mutate(comb=paste(MCls, contrast, sep="_"))


## ###
## ### qq for cloest genes
## comb <- sort(unique(resDE$comb))
## dfNew <- map_dfr(comb, function(ii){
##   res2 <- res%>%dplyr::filter(comb==ii)
##   DEG <- resDE%>%dplyr::filter(comb==ii)%>%dplyr::pull(gene)
##   res2 <- res2%>%mutate(is_DEG=ifelse(geneId%in%DEG,1,0))  
##   ###
##   dx <- map_dfr(c(0,1),function(i){
##      di <- res2%>%dplyr::filter(is_DEG==i)%>%arrange(p.value)
##      ngene <- nrow(di)
##      di <- di%>%mutate(observed=-log10(p.value), expected=-log10(ppoints(ngene)))
##      di
##   })
##   dx
## })

## ###
## ###
## lab1 <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX", "PHA"="PHA", "PHA-DEX"="PHA+DEX")
## p1 <- ggplot(dfNew, aes(x=expected,y=observed, colour=factor(is_DEG)))+
##    ggrastr::rasterise(geom_point(size=0.3),dpi=300)+
##    geom_abline(colour="red")+
##    scale_colour_manual(values=c("0"="grey40","1"="green"),
##       labels=c("0"="Not DVG", "1"="DVG"),
##       guide=guide_legend(override.aes=list(size=1)))+
##    facet_grid(MCls~contrast, scales="free", labeller=labeller(contrast=lab1))+
##    xlab(bquote("Expected"~ -log[10]~"("~italic(plain(P))~")"))+
##    ylab(bquote("Observed"~-log[10]~"("~italic(plain(P))~")"))+
##    theme_bw()+
##    theme(legend.title=element_blank(),strip.text=element_text(size=12))

## figfn <- "./2.2_compareRNAandATAC.outs/Figure2.5_qq.DVG.cloest.png"
## png(figfn, width=750, height=750, res=120)
## print(p1)
## dev.off()



###############################################
### forest plots DEG if is enriched in DARs ###
###############################################
### fisher test function
cal.fisher <- function(df){
   ###
 resfisher <- map_dfr(1:nrow(df),function(i){
      dmat <- matrix(as.numeric(df[i,]), 2, 2)
      colnames(dmat) <- c("interest", "not.interest")
      rownames(dmat) <- c("in.DAR", "not.DAR")
      ##rownames(dmat) <- c("in.DEG", "not.DEG")
      res <- fisher.test(dmat)
      res2 <- data.frame(odds=res$estimate, pval.fisher=res$p.value,
         lower=res$conf.int[1], upper=res$conf.int[2])
      res2
   })
   resfisher
}

## ##
## df.fisher <- function(resDE, resDP, peakAnno){
##    ### 
##    comb <- sort(unique(resDE$comb))
##    df.ls <- lapply(1:length(comb), function(i){
##       ii <- comb[i] 
##       cell <- gsub("_.*", "", ii)
##       contrast <- gsub(".*_", "", ii)
##       resDE2 <- resDE%>%dplyr::filter(comb==ii) ## dplyr::filter(qval<0.1, abs(beta)>0.5)
##       resDP2 <- resDP%>%
##          dplyr::filter(comb==ii, p.adjusted<0.1, abs(estimate)>0.5)%>%
##          dplyr::left_join(peakAnno, by="peak")
##       if( nrow(resDP2)>5){
##          sigGene <- resDE2%>%
##             dplyr::filter(qval<0.1, abs(beta)>0.5)%>%
##             dplyr::pull(gene)
##          ###
##          interest.in.DARs <- sum(sigGene%in%resDP2$geneID)
##          interest.not.DARs <- length(sigGene)-interest.in.DARs
##       ###
##          notSig <- setdiff(res2$gene, sigGene)
##          not.interest.in.DARs <- sum(notSig%in%resDP2$geneId)
##          not.interest.not.DARs <- length(notSig)-not.interest.in.DARs
##          df2 <- data.frame("interest.in.DARs"=interest.in.DARs,
##             "interest.not.DARs"=interest.not.DARs,
##             "not.interest.in.DARs"=not.interest.in.DARs,
##             "not.interest.not.DARs"=not.interest.not.DARs)
##          df2$cell <- cell
##          df2$contrast <- contrast
##          df2$comb <- ii
##       }else{
##         df2 <- NA    
##       }
##       df2
##    })    
##    df.ls <- df.ls[!is.na(df.ls)]
##    df <- do.call(rbind, df.ls)    
##    res <- cal.fisher(df[,1:4])
##    res2 <- cbind(df, res)
##    ### 
##    as.data.frame(res2)
## }    

head(resDP.peakAnno)
nrow(resDP.peakAnno)
colnames(resDP.peakAnno)

df.fisher <- function(resDE, resDP.peakAnno){
   ### 
   comb <- sort(unique(resDE$comb))
   df.ls <- lapply(1:length(comb), function(i){
      ii <- comb[i] 
      contrast <- gsub("_.*", "", ii)
      cell <- gsub(".*_", "", ii)
      ##
      resDE2 <- resDE%>%dplyr::filter(comb==ii) ## dplyr::filter(qval<0.1, abs(beta)>0.5)
      ##resDP2 <- resDP%>%
      ##   dplyr::filter(comb==ii, p.adjusted<0.1, abs(estimate)>0.5)%>%
      ##   dplyr::left_join(peakAnno, by="peak")
      resDP2 <- resDP.peakAnno%>%dplyr::filter(comb==ii)
      ##
      ##
      if(nrow(resDP2)>5){
          ##
          ##
          sigGene <- resDE2%>%dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)%>%
              dplyr::pull(gene)
          sigPeak <- unique(resDP2%>%dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)%>%dplyr::pull(symbol))
          ##nonsigPeak <- unique(resDP2%>%dplyr::filter(p.adjusted>0.1 | (p.adjusted<0.1 & abs(estimate)>0.5))%>%dplyr::pull(symbol))
          notSig <- setdiff(resDE2$gene, sigGene)
          ##
          ##
          interest.in.DARs <- sum(sigGene%in%sigPeak)  #%in%resDP2$symbol)
          interest.not.DARs <- length(sigGene)-interest.in.DARs
          ##interest.in.DEGs <- sum(sigPeak%in%sigGene)
          ##interest.not.DEGs <- sum(sigPeak%in%notSig)
          ##interest.in.DEGs <- sum(sigPeak%in%sigGene)
          ##interest.not.DEGs <- length(sigPeak)-interest.in.DEGs
          ##
          ##
          not.interest.in.DARs <- sum(notSig%in%sigPeak) #%in%resDP2$symbol)
          not.interest.not.DARs <- length(notSig)-not.interest.in.DARs
          ##not.interest.in.DEGs <- sum(nonsigPeak%in%sigGene)
          ##not.interest.not.DEGs <- sum(nonsigPeak%in%notSig)
          ##not.interest.in.DEGs <- sum(nonsigPeak%in%sigGene)
          ##not.interest.not.DEGs <- length(nonsigPeak)-not.interest.in.DEGs
          ##
          ##
          df2 <- data.frame("interest.in.DARs"=interest.in.DARs,
                            "interest.not.DARs"=interest.not.DARs,
                            "not.interest.in.DARs"=not.interest.in.DARs,
                            "not.interest.not.DARs"=not.interest.not.DARs)
          #df2 <- data.frame("interest.in.DEGs"=interest.in.DEGs,
          #                  "interest.not.DEGs"=interest.not.DEGs,
          #                  "not.interest.in.DEGs"=not.interest.in.DEGs,
          #                  "not.interest.not.DEGs"=not.interest.not.DEGs)
          ##
          ##
          df2$cell <- cell
          df2$contrast <- contrast
          df2$comb <- ii
      }else{
          df2 <- NA    
      }
      df2
   })    
    df.ls <- df.ls[!is.na(df.ls)]
    df <- do.call(rbind, df.ls)    
    res <- cal.fisher(df[,1:4])
    res2 <- cbind(df, res)
    ##
    as.data.frame(res2)
}
    

resFisher <- df.fisher(resDE, resDP.peakAnno)
head(resFisher)

#opfn <- "./2.2_compareRNAandATAC_3.outs/3.1_enrich.DEG_control.csv"
#write.csv(resFisher, opfn, row.names=F)

opfn <- "./2.2_compareRNAandATAC_3.outs/3.1_enrich.DEG_howmanygenes_control.csv"
write.csv(resFisher, opfn, row.names=F)



#################
### read data ###
#################
## fn <- "./1.2_DiffPeak.outs/2.0_DESeq.results_control.rds"
## resDP <- read_rds(fn)%>%drop_na(p.value)%>%
##    mutate(comb=paste(MCls, contrast, sep="_"))%>%
##    dplyr::rename("peak"="gene") 


## fn <- "./2.2_compareRNAandATAC.outs/2_annot.ChIPseeker.rds"
## peakAnno <- read_rds(fn)%>%as.data.frame()%>%
##    mutate(peak=paste(gsub("chr", "", seqnames), start, end, sep="-")) 
## peakAnno2 <- peakAnno%>%dplyr::select(peak, geneId, SYMBOL)


## ### test if DEG is DARs 
## fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
## res <- read_rds(fn)%>%drop_na(pval)%>%
##    mutate(comb=paste(MCls, contrast, sep="_"))%>%
##    as.data.frame()
## ##
## resFisher <- df.fisher(res, resDP, peakAnno2)
## opfn <- "./2.2_compareRNAandATAC.outs/3.1_enrich.DEG.csv"
## write.csv(resFisher, opfn, row.names=F)


## ### test if DVGs is DARs
## fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/10_RNA.Variance_output/tmp9/3_phiNew.meta"
## res <- read.table(fn, header=T) %>%drop_na(pval)%>%
##    mutate(comb=paste(MCls, contrast, sep="_"))

## resFisher2 <- df.fisher(res, resDP, peakAnno2)
## opfn <- "./2.2_compareRNAandATAC.outs/3.2_enrich.DVG.csv"
## write.csv(resFisher2, opfn, row.names=F)







#######################
### error bar plots ###
#######################
## col.treat <- c("LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
##              "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
## col.MCl <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
##              "NKcell"="#aa4b56", "Tcell"="#ffaa00")




col.treat <- c("caffeine"="red", "nicotine"="tan",
          "vitA"="tan4", "vitD"="seagreen4",
          "vitE"="salmon3",
          ##"water"="grey",
          "zinc"="maroon3")
col.MCl <- c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
          "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00", "1-TCM"="pink",
          "3-TEM"="blue", "5-CD8Naive"="green", "7-dnT"="black")



## lab1 <- c("caffeine"="caffeine", "nicotine"="nicotine", "zinc"="zinc",
##           "vitA"="vitA", "vitD"="vitD", "vitE"="vitE", "water"="water")

## lab2 <- c("4_Bcell"="4_Bcell", "6_Monocyte"="6_Monocyte",
##           "2_NKcell"="2_NKcell", "0_CD4Naive"="0_CD4Naive", "1_TCM"="1_TCM",
##           "3_TEM"="3_TEM", "5_CD8Naive"="5_CD8Naive", "7_dnT"="7_dnT")




fn <- "./2.2_compareRNAandATAC_3.outs/3.1_enrich.DEG_control.csv"
df <- read.csv(fn, header=T)
head(df)
nrow(df)

df <- resFisher
head(df)

## no need to calculate qval as only one row for each contrast-cell combination
#df <- df%>%group_by(cell, contrast)%>%
#       mutate(qvalue.fisher=p.adjust(pval.fisher, "BH"))%>%
#       ungroup()%>%
#       as.data.frame()
#head(df)
#nrow(df)

#df <- df %>% dplyr::filter(pval.fisher<0.1)
#nrow(df)


df%>%dplyr::filter(contrast=="caffeine")
df%>%dplyr::filter(contrast=="zinc")


df1 <- data.frame(odds=df$odds,
   CI.low=df$lower, CI.high=df$upper,
   comb=gsub("_", ".", df$comb),
   MCls=df$cell, contrast=df$contrast, gene="DEG")


## df1.notinf <- df1 %>% dplyr::filter(CI.high!=Inf)
## head(df1.notinf)
## error.max <- max(df1.notinf$CI.high)
## error.max
## high.quant <- quantile(abs(df1.notinf$CI.high),probs=0.75)
## high.quant

head(df1)

## fn <- "./2.2_compareRNAandATAC.outs/3.2_enrich.DVG.csv"
## df <- read.csv(fn, header=T)
## df2 <- data.frame(odds=df$odds,
##    CI.low=df$lower, CI.high=df$upper,
##    comb=gsub("_", ".", df$comb),
##    MCls=df$cell, contrast=df$contrast, gene="DVG")

## df2 <- rbind(df1, df2)

df2 <- df1

comb2 <- sort(unique(df2$comb))
comb2

##ylab2 <- gsub("-", "+", gsub(".*\\.", "", df2$comb))
ylab2 <- gsub(".*\\.", "", df2$comb)
names(ylab2) <- comb2
ylab2

p <- ggplot(df2, aes(x=odds, y=comb))+
   ##geom_errorbarh(aes(xmax=CI.high, xmin=CI.low, colour=MCls),
   ##    size=0.5, height=0.2)+
   geom_point(aes(colour=contrast), shape=19, size=1.5)+
   scale_colour_manual(values=col.treat)+
   geom_vline(aes(xintercept=1), size=0.25, linetype="dashed")+ 
   xlab("Odds ratio")+
   scale_y_discrete(labels=ylab2)+
   facet_wrap(~factor(gene), ncol=2)+ 
   theme_bw()+
   theme(plot.title=element_text(hjust=0.5, size=9),
         axis.title.y=element_blank(),
         axis.title.x=element_text(size=9),
         ## axis.text=element_text(size=8),
         legend.position="right")

###
figfn <- "./2.2_compareRNAandATAC_3.outs/Figure3.0_enrich.DAR_howmanygenes_control.png"
png(figfn, width=1500, height=1500, res=175)
p
dev.off()



##
### poster
p <- ggplot(df2%>%dplyr::filter(gene=="DEG"), aes(y=odds, x=comb))+
   geom_errorbar(aes(ymax=CI.high, ymin=CI.low, colour=MCls),
       size=0.5, width=0.2)+
   geom_point(aes(colour=MCls), shape=19, size=1.5)+
   scale_colour_manual(values=col.MCl)+
   geom_hline(aes(yintercept=1), size=0.25, linetype="dashed")+ 
   ylab("Odds ratio")+
   scale_x_discrete(labels=ylab2)+
   ## facet_wrap(~factor(gene), ncol=2)+ 
   theme_bw()+
   theme(plot.title=element_text(hjust=0.5, size=9),
         axis.title.x=element_blank(),
         axis.title.y=element_text(size=12),
         axis.text.x=element_text(size=10, angle=-90, hjust=0, vjust=0.5),
         axis.text.y=element_text(size=10),
         legend.position="none")

###
figfn <- "./2.2_compareRNAandATAC.outs/Figure3.1_enrich.DAR.png"
png(figfn, width=620, height=480, res=120)
p
dev.off()


###
p2 <- ggplot(df2%>%filter(gene=="DEG"), aes(x=odds, y=comb))+
   geom_errorbarh(aes(xmax=CI.high, xmin=CI.low, colour=MCls),
       size=0.5, height=0.2)+
   geom_point(aes(colour=MCls), shape=19, size=1.5)+
   scale_colour_manual(values=col.MCl)+
   geom_vline(aes(xintercept=1), size=0.25, linetype="dashed")+ 
   xlab("Odds ratio")+
   scale_y_discrete(labels=ylab2)+
   ## facet_wrap(~factor(gene), ncol=2)+ 
   theme_bw()+
   theme(plot.title=element_text(hjust=0.5, size=9),
         axis.title.y=element_blank(),
         axis.title.x=element_text(size=9),
         ## axis.text=element_text(size=8),
         legend.position="none")
###
figfn <- "./2.2_compareRNAandATAC.outs/Figure3.4_enrich.DAR.png"
png(figfn, width=380, height=480, res=120)
p2
dev.off()

##                  
## p1 <- ggplot(df2, aes(x=odds, y=comb))+
##    geom_errorbarh(aes(xmax=CI.high, xmin=CI.low, colour=MCls),
##        size=0.5, height=0.2)+
##    geom_point(aes(colour=MCls), shape=19, size=1.5)+
##    scale_colour_manual(values=col.MCl)+
##    geom_vline(aes(xintercept=1), size=0.25, linetype="dashed")+ 
##    xlab("Odds ratio")+
##    ## scale_y_discrete(labels=ylab2)+ 
##    ggtitle("DEGs")+    
##    theme_bw()+
##    theme(plot.title=element_text(hjust=0.5, size=9),
##          axis.title.y=element_blank(),
##          axis.title.x=element_text(size=9),
##          axis.text=element_text(size=8),
##          legend.position="none")

## ###
## fn <- "./2.2_compareRNAandATAC.outs/3.2_enrich.DVG.csv"
## df <- read.csv(fn, header=T)
## df2 <- data.frame(odds=df$odds,
##    CI.low=df$lower, CI.high=df$upper,
##    comb=gsub("_", ".", df$comb),
##    MCls=df$cell, contrast=df$contrast)

## p2 <- ggplot(df2, aes(x=odds, y=comb))+
##    geom_errorbarh(aes(xmax=CI.high, xmin=CI.low, colour=MCls),
##        size=0.5, height=0.2)+
##    geom_point(aes(colour=MCls), shape=19, size=1.5)+
##    scale_colour_manual(values=col.MCl)+
##    geom_vline(aes(xintercept=1), size=0.25, linetype="dashed")+ 
##    xlab("Odds ratio")+
##    ## scale_y_discrete(labels=ylab2)+   
##    ggtitle("DVGs")+    
##    theme_bw()+
##    theme(plot.title=element_text(hjust=0.5, size=9),
##          axis.title.y=element_blank(),
##          axis.title.x=element_text(size=9),
##          axis.text=element_text(size=8),
##          legend.position="none")

## ###
## figfn <- "./2.2_compareRNAandATAC.outs/Figure3_enrich.DAR.png"
## png(figfn, width=680, height=480, res=120)
## plot_grid(p1, p2, ncol=2)
## dev.off()



#######################
### Distance to TSS ###
#######################

###
### DEGs or DVGs are Distance to peaks

###
### differential peaks
fn <- "./1.2_DiffPeak.outs/2.0_DESeq.results.rds"
resDP <- read_rds(fn)%>%drop_na(p.value)%>%
   mutate(comb=paste(MCls, contrast, sep="_"))%>%
   dplyr::rename("peak"="gene")%>%
   dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)%>%
   as.data.frame()

peakAll <- unique(resDP$peak)


###
### annotation
fn <- "./2.2_compareRNAandATAC.outs/2_annot.ChIPseeker.rds"
peakAnno <- read_rds(fn)%>%as.data.frame()%>%
   mutate(peak=paste(gsub("chr", "", seqnames), start, end, sep="-")) 
peakAnno2 <- peakAnno%>%dplyr::filter(peak%in%peakAll)%>%
   dplyr::select(peak, geneId, SYMBOL, dtss=distanceToTSS)%>%mutate(dtss=abs(dtss))


 ## df2gene <- peakAnno2%>%
 ##     group_by(geneId)%>%summarise(npeaks=n(), .groups="drop")%>%
 ##     ungroup()%>%as.data.frame()


###
### DEG is DARs 
fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
resDEG <- read_rds(fn)%>%drop_na(pval)%>%
   mutate(comb=paste(MCls, contrast, sep="_"))%>%
   dplyr::filter(qval<0.1, abs(beta)>0.5)%>% 
   as.data.frame()

## DEG <- res%>%dplyr::filter(qval<0.1, abs(beta)>0.5)%>%
##     dplyr::pull(gene)%>%unique()
##

###
### DVG
fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/10_RNA.Variance_output/tmp9/3_phiNew.meta"
resDVG <- read.table(fn, header=T) %>%drop_na(pval)%>%
   mutate(comb=paste(MCls, contrast, sep="_"))%>%
   dplyr::filter(qval<0.1, abs(beta)>0.5)%>%
   as.data.frame()


###
###
comb <- sort(unique(resDEG$comb))
geneList <- lapply(comb, function(ii){
###
   resDP2 <- resDP%>%dplyr::filter(comb==ii)
   resDE2 <- resDEG%>%dplyr::filter(comb==ii)
   resDV2 <- resDVG%>%dplyr::filter(comb==ii)
###
   DP <- unique(resDP2$peak)
   DEG <- unique(resDE2$gene)
   DVG <- unique(resDV2$gene)
    
   DEGunq <- setdiff(DEG, DVG)
   shared <- intersect(DEG, DVG)
   DVGunq <- setdiff(DVG, DEG)
   DGall <- union(DEG, DVG)

   ## DEG 
   x <- peakAnno2%>%dplyr::filter(peak%in%DP, geneId%in%DEG)
   if(nrow(x)==0){
      res1 <- NA
    }else{
       res1 <- unique(x$geneId)
    }
    ##
    x <- peakAnno2%>%dplyr::filter(peak%in%DP, geneId%in%DVG)
    if(nrow(x)==0){
      res2 <- NA
    }else{
       res2 <- unique(x$geneId)
    }
    res <- list(DEG=res1, DVG=res2)
})

gene1 <- lapply(geneList, function(x) x[[1]])
gene1 <- unique(unlist(gene1[!is.na(gene1)]))

gene2 <- lapply(geneList, function(x) x[[2]])
gene2 <- unique(unlist(gene2[!is.na(gene2)]))

DEG <- unique(resDEG$gene)
DVG <- unique(resDVG$gene)
## DEG, 3115/6571, 47%
## DVG, 579/1409, 41%




###
### compute the number of differential peaks nearby DEG or DVG genes and distance to DPs of DEG or DVG 
comb <- sort(unique(resDEG$comb))
dd <- map_dfr(comb, function(ii){
###
   resDP2 <- resDP%>%dplyr::filter(comb==ii)
   resDE2 <- resDEG%>%dplyr::filter(comb==ii)
   resDV2 <- resDVG%>%dplyr::filter(comb==ii)
###
   DP <- unique(resDP2$peak)
   DEG <- unique(resDE2$gene)
   DVG <- unique(resDV2$gene)
    
   DEGunq <- setdiff(DEG, DVG)
   shared <- intersect(DEG, DVG)
   DVGunq <- setdiff(DVG, DEG)
   DGall <- union(DEG, DVG)
   
###
  
   ## d0 <- peakAnno2%>%dplyr::filter(peak%in%DP, !geneId%in%DGall)
   ## d0 <- d0%>%group_by(geneId)%>%
   ##    summarise(npeaks=n(), dtss=median(dtss), .groups="drop")%>%
   ##    ungroup()%>%as.data.frame()%>%
   ##    mutate(conditions=ii, grp="0")
    
  ###  
  d1 <- peakAnno2%>%dplyr::filter(peak%in%DP, geneId%in%DEG)  
  d1 <- d1%>%group_by(geneId)%>%
     summarise(npeaks=n(), dtss=min(dtss), .groups="drop")%>%
     ungroup()%>%as.data.frame()%>%
     mutate(conditions=ii, grp="1")
  ##
  ## d2 <- peakAnno2%>%dplyr::filter(peak%in%DP, geneId%in%shared)  
  ## d2 <- d2%>%group_by(geneId)%>%
  ##    summarise(npeaks=n(), dtss.min=min(dtss),.groups="drop")%>%
  ##    ungroup()%>%as.data.frame()%>%
  ##    mutate(conditions=ii, grp="2")
  ##
  d2 <- peakAnno2%>%dplyr::filter(peak%in%DP, geneId%in%DVG)  
  d2 <- d2%>%group_by(geneId)%>%
     summarise(npeaks=n(), dtss=min(dtss), .groups="drop")%>%
     ungroup()%>%as.data.frame()%>%
     mutate(conditions=ii, grp="2")
  ##
  dd <- rbind(d1, d2)  
#
  dd    
})    


          
###
p1 <- ggplot(dd, aes(x=factor(grp), y=log2(npeaks), fill=factor(grp)))+
   geom_violin()+
   scale_x_discrete(labels=c("0"="background", "1"="DEG", "3"="DVG-only"))+ 
   ## scale_fill_manual(values=c("DEG"="#d95f02", "DVG"="#1b9e77"))+
   ## scale_color_manual(values=c("DEG"="#d95f02", "DVG"="#1b9e77"))+
   ylab(bquote("Complexity ("~log[2]~" peaks per gene)"))+ 
   theme_bw()+
   theme(legend.position="none",
         axis.title.x=element_blank(),
         legend.background=element_blank()) 

figfn <- "./2.2_compareRNAandATAC.outs/Figure4.1_complexity_violin.png"
png(figfn, width=480, height=320, res=120)
print(p1)
dev.off()

          
###
p2 <- ggplot(dd, aes(x=factor(grp), y=log2(dtss+0.1), fill=factor(grp)))+
   geom_violin()+
   scale_x_discrete(labels=c("0"="background", "1"="DEG", "3"="DVG-only"))+ 
   ## scale_fill_manual(values=c("DEG"="#d95f02", "DVG"="#1b9e77"))+
   ## scale_color_manual(values=c("DEG"="#d95f02", "DVG"="#1b9e77"))+
   ylab(bquote("Distance to neary peaks"~log[2]))+ 
   theme_bw()+
   theme(legend.position="none",
         axis.title.x=element_blank(),
         legend.background=element_blank()) 

figfn <- "./2.2_compareRNAandATAC.outs/Figure4.2_dtss_violin.png"
png(figfn, width=480, height=320, res=120)
print(p2)
dev.off()


###
### summary number of peaks
dd1 <- dd%>%dplyr::filter(grp==1)
sum(dd1$npeaks<=10)/nrow(dd1)
sum(abs(dd1$dtss)<3000)/nrow(dd1)
## dd1 <- dd1%>%group_by(npeaks)%>%summarise(freq=n(),.groups="drop")%>%
##    mutate(prop=freq/sum(freq))%>%
##    ungroup()%>%as.data.frame()

## x1 <- data.frame(order=1, freq=sum(dd1[dd1$npeaks>=1&dd1$npeaks<8,2]), prop=sum(dd1[dd1$npeaks>=1&dd1$npeaks<8,3]))
## ## x2 <- data.frame(order=2, freq=sum(dd1[dd1$npeaks>=5&dd1$npeaks<11,2]), prop=sum(dd1[dd1$npeaks>=5&dd1$npeaks<11,3]))
## x2 <- data.frame(order=2, freq=sum(dd1[dd1$npeaks>=8,2]), prop=sum(dd1[dd1$npeaks>=8,3]))                            
## a <- rbind(x1,x2)
## a$grp <- 1

dd2 <- dd%>%dplyr::filter(grp==2)
sum(dd2$npeaks<=10)/nrow(dd2)
sum(abs(dd2$dtss)<3000)/nrow(dd2)


###
