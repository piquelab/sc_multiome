###
## library("rhdf5")
## library("corpcor")
library(Matrix)
## library(MASS)
## library(scales)
library(tidyverse)
## library(parallel)
## library(data.table)
## library(future)
## library(purrr)
## library(furrr)
## library("BiocParallel")
## library(Rcpp)
## library(reshape)
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
library(clusterProfiler)
library(org.Hs.eg.db)
library(annotables)  
###
library(ggplot2)
#library(cowplot,lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(cowplot)
library(grid)
library(gridExtra)
library(ggExtra)
library(RColorBrewer)
library(ggsci)
library(RColorBrewer)
library(viridis)
theme_set(theme_grey())

#outdir <- "./1_DiffRNA_2.outs/"
#if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)

### Last modified by Mohammed Husain Bharmal

####################################
### sys_args ###
####################################
reso = 0.1

####################################
### 1. Generate pseudo-bulk data ###
####################################
combined <- read_rds("../1_processing/3_Clustering.outs/1_seurat.cluster.combined.mitofilt.rds")
#colnames(combined@meta.data)

#combined <- subset(x = combined, subset = percent.mt.RNA < 20)
#colnames(combined@meta.data)
#max(combined@meta.data$percent.mt.RNA)

#--------------------------------------------------------------------------
## meta <- atac@meta.data%>%mutate(treat=gsub(".*-ATAC-|_.*", "", NEW_BARCODE)) 
meta <- combined@meta.data

meta$treat <- meta$treats

table(meta$wsnn_res.0.1)

#colnames(meta)
#table(meta$wsnn_res.0.05)

#----cell labeling----
# cell <- c("0"="Tcell", "1"="Tcell", "2"="NKcell", "3"="Bcell",
# "4"="Tcell", "5"="Tcell", "6"="Monocyte", "7"="Tcell", "8"="Tcell",
# "9"="Tcell", "10"="Tcell", "11"="DC", "12"="DC")

## cell <- c("0"="CD4Naive", "1"="TCM", "2"="NKcell", "3"="Bcell",
##           "4"="CD8Naive", "5"="TEM", "6"="Monocyte", "7"="CTL", "8"="dnT",
##           "9"="Treg", "10"="Platelet", "11"="DC", "12"="HSPC", "13"="Bcell")


## cell <- c("0"="CD4Naive", "1"="TCM", "2"="NKcell", "3"="Bcell",
##           "4"="CD8Naive", "5"="TEM", "6"="Monocyte", "7"="CTL", "8"="dnT",
##           "9"="Treg", "10"="Platelet", "11"="DC")

cell <- c("0"="CD4Naive", "1"="TCM", "2"="NKcell", "3"="TEM",
           "4"="Bcell", "5"="CD8Naive", "6"="Monocyte", "7"="dnT", "8"="MAIT",
           "9"="Platelet", "10"="DC")

#cell <- c("0"="Tcell", "1"="Tcell", "2"="Tcell", "3"="NKcell",
#         "4"="Bcell", "5"="Monocyte", "6"="Platelet", "7"="DC")


#meta$MCls <-NULL
meta$MCls <- cell[as.character(meta$wsnn_res.0.1)]
colnames(meta)
#table(meta$MCls)
#----cell labeling----


meta <- meta%>%
   mutate(bti=paste(wsnn_res.0.1, MCls, treat, SNG.BEST.GUESS, sep="_"))%>%
   dplyr::select(NEW_BARCODE, bti)
#colnames(meta)
head(meta)

dd <- meta%>%group_by(bti)%>%summarise(ncell=n(),.groups="drop")
dd

#write_rds(dd, "./1_DiffRNA_2.outs/0_ncell_clusters.rds")

#dd <- read_rds("./1_DiffRNA_2.outs/0_ncell_clusters.rds")
#dd


#----------------------------------------------------------------------------
count <- combined@assays$RNA@counts
#head(count)
#nrow(count)

count_genes <- rownames(count)
#head(count_genes)
#length(count_genes)





# filtering X and Y
#-------------------------------------------------------------------
library(biomaRt)

mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
#mart

genes.table <- try(biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description", "gene_biotype", "chromosome_name", "start_position"), mart = mart, useCache = F))
#colnames(genes.table)
#table(genes.table$chromosome_name)
#nrow(genes.table)

genes.table.2 <- genes.table%>%dplyr::filter(chromosome_name %in% c("X", "Y"), external_gene_name != "", external_gene_name %in% count_genes)
#table(genes.table.2$chromosome_name)
#colnames(genes.table.2)
#nrow(genes.table.2)

no_genes <- genes.table.2$external_gene_name
length(no_genes)

#check <- no_genes %in% count_genes
#check

#--------------------------------------------------------------
yes_genes <- count_genes[!(count_genes %in% no_genes)]

#length(yes_genes)

count.2 <- count[yes_genes,]
#head(count.2)
#head(rownames(count.2))
nrow(count.2)
ncol(count.2)

#---------------------------------------------------------------
anno <- data.frame(rn=rownames(count.2), rnz=rowSums(count.2))
#str(anno)
head(anno)

#autosome <- paste("chr", as.character(1:22), sep="")
#autosome

annoSel <- anno%>%dplyr::filter(rnz>0)
head(annoSel)

Y <- count[annoSel$rn,]
head(Y)


#------------------------------------------------------------------
##pseudo-bulk peak data
#length(meta$bti)
bti <- factor(meta$bti)
#head(bti)
#length(bti)
sort(unique(bti))

X <- model.matrix(~0+bti)
head(X)
nrow(X)
ncol(X)

YtX <- Y %*% X
YtX <- as.matrix(YtX)

colnames(YtX) <- gsub("^bti", "", colnames(YtX))
#colnames(YtX)
head(YtX)
nrow(YtX)
ncol(YtX)

###rnz>0,chr:1-22, 111,750*400
#opfn <- "./1_DiffRNA_2.outs/1_YtX.comb_clusters.rds" 
#write_rds(YtX, file=opfn)



#-------------------------------------------------------------------
###rnz>100, ncell>20, 111,746 * 331
#dd <- read_rds("./1_DiffRNA_2.outs/0_ncell_clusters.rds")
#dd <- read_rds("./1_DiffRNA_2.outs/0_ncell_clusters.rds")%>%filter(ncell>20)

dd <- dd %>%filter(ncell>20)
#dd
#table(dd$bti)

## YtX <- read_rds("./1_DiffRNA_2.outs/1_YtX.comb_clusters.rds")
## head(YtX)
## ncol(YtX)
## nrow(YtX)

anno <- data.frame(rn=rownames(YtX), rnz=rowSums(YtX))
head(anno)

annoSel <- anno%>%filter(rnz>100)
#head(annoSel)

YtX_sel <- YtX[annoSel$rn, dd$bti]
head(YtX_sel)

nrow(YtX_sel)
ncol(YtX_sel)

opfn <- paste0("./1_DiffRNA_2.outs/1_YtX.sel_clusters_", reso, "_cn.rds")
write_rds(YtX_sel, file=opfn)

YtX_sel <- read_rds("./1_DiffRNA_2.outs/1_YtX.sel_clusters_0.1_cn.rds")

#####################################################
### 2. Differential analysis for CRP, High vs Low ###
#####################################################
### function, adjusted p value
## myqval <- function(pval){
##   qval <- pval
##   ii0 <- !is.na(pval)
##   qval[ii0] <- qvalue(pval[ii0],pi0=1)$qvalues
##   qval
## }


#### Read data
#YtX <- read_rds("./1_DiffRNA_2.outs/1_YtX.sel_clusters.rds")

YtX  <- YtX_sel
ncol(YtX)
nrow(YtX)
#str(YtX)
#colnames(YtX)

bti2 <- colnames(YtX)
#length(bti2)

cvt0 <- str_split(bti2, "_", simplify=T)
#head(cvt0)

cvt <- data.frame(bti=bti2, MCls=paste(cvt0[,1], cvt0[,2], sep="_"), treat=cvt0[,3], sampleID=cvt0[,4])
#cvt

#comb <- table(cvt$MCls, cvt$treat)
#write.csv(comb, paste0("./1_DiffRNA_2.outs/combinations.csv"), quote = FALSE, row.names = TRUE)
#comb <- read.csv("./1_DiffRNA_2.outs/combinations.csv")
#comb
#dd2 <- dd %>% dplyr::filter(bti %in% bti2)
#sum(dd2$ncell)




# filtered
#-------------------------------------------------------------------------
cvt <- cvt%>%mutate(MCls_IDs=paste(cvt$MCls, ".", cvt$sampleID, sep="")) 
head(cvt)

cvt <- cvt%>%mutate(treat_IDs=paste(cvt$treat, ".", cvt$sampleID, sep="")) 
head(cvt)

cvt <- cvt%>%mutate(treat_MCls=paste(cvt$treat, ".", cvt$MCls, sep="")) 
head(cvt)

MCls.IDs <- unique(cvt$MCls_IDs)
#MCls.IDs
treat.IDs <- unique(cvt$treat_IDs)
#treat.IDs 
treat.MCls <- unique(cvt$treat_MCls)
#treat.MCls

####
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
head(cvt.2)

cvt.3 <- cvt.2 %>% mutate(vehicle=treat)
#cvt.3 <- cvt.3 %>% mutate(vehicle=gsub("etOH", "etOH", vehicle))
#cvt.3 <- cvt.3 %>% mutate(vehicle=gsub("water", "etOH", vehicle))
cvt.3 <- cvt.3 %>% mutate(vehicle=gsub("caffeine", "water", vehicle))
cvt.3 <- cvt.3 %>% mutate(vehicle=gsub("nicotine", "water", vehicle))
cvt.3 <- cvt.3 %>% mutate(vehicle=gsub("zinc", "water", vehicle))
cvt.3 <- cvt.3 %>% mutate(vehicle=gsub("vitA", "etOH", vehicle))
cvt.3 <- cvt.3 %>% mutate(vehicle=gsub("vitD", "etOH", vehicle))
cvt.3 <- cvt.3 %>% mutate(vehicle=gsub("vitE", "etOH", vehicle))

cvt.3 <- cvt.3 %>% mutate(new_treat=treat)
cvt.3 <- cvt.3 %>% mutate(new_treat=gsub("etOH", "control", new_treat))
cvt.3 <- cvt.3 %>% mutate(new_treat=gsub("water", "control", new_treat))

head(cvt.3)

comb.2 <- table(cvt.2$MCls, cvt.2$treat)
comb.2
 
##write.csv(comb.2, paste0("./1_DiffRNA_2.outs/combinations_2_", reso, "_cn.csv"), quote = FALSE, row.names = TRUE)
##comb.2 <- read.csv("./1_DiffRNA_2.outs/combinations_2.csv")
##comb.2

bti3 <- cvt.2$bti 
##bti3

##ncol(YtX)
YtX <- YtX[,bti3]
##ncol(YtX)

dd2 <- dd %>% dplyr::filter(bti %in% bti3)
##dd2
##sum(dd2$ncell)



## ##----------------------------------------------------------------------------
## ##-----------------------------trial------------------------------------------
## ##----------------------------------------------------------------------------
## head(YtX_sel)


## bti <- colnames(YtX_sel)
## ##sort(unique(bti))
## ##length(bti)
## ##head(bti)

## cvt.trial  <- str_split(bti, "_", simplify=T) %>% as.data.frame() %>% mutate(comb=paste(V1, V2, V3, V4, sep="_"))
## head(cvt.trial)
## ##nrow(cvt.trial)


## cvt.trial.2 <- cvt.trial %>% dplyr::filter(V2=="Monocyte")

## head(cvt.trial.2)

## cell <- factor(cvt.trial.2$V2)
## ##cell
## treat <- factor(cvt.trial.2$V3)
## ##treat
## ind <- factor(cvt.trial.2$V4)
## ##ind


## YtX_sel.2 <- YtX_sel[, cvt.trial.2$comb]
## head(YtX_sel.2)


## X <- model.matrix(~0+treat+ind, data=cvt.trial.2)
## head(X)
## nrow(X)
## ncol(X)

## YtX_sel.3 <- YtX_sel.2 %*% X
## YtX_sel.3 <- as.matrix(YtX_sel.3)

## colnames(YtX_sel.3) <- gsub("^cell", "", colnames(YtX_sel.3))
## colnames(YtX_sel.3) <- gsub("^treat", "", colnames(YtX_sel.3))
## colnames(YtX_sel.3) <- gsub("^ind", "", colnames(YtX_sel.3))

## head(YtX_sel.3)
## nrow(YtX_sel.3)
## ncol(YtX_sel.3)



#################################################################
### 2.1 call DESeq ###
#################################################################
contrast.list <- list("caffeine"=c("treat", "caffeine", "water"),
                      "nicotine"=c("treat", "nicotine", "water"),
                      "vitA"=c("treat", "vitA", "etOH"),
                      "vitD"=c("treat", "vitD", "etOH"),
                      "vitE"=c("treat", "vitE", "etOH"),
                      "water"=c("treat", "water", "etOH"),
                      "zinc"=c("treat", "zinc", "water"))


contrast.list <- list("caffeine"=c("treat", "caffeine", "etOH"),
                      "nicotine"=c("treat", "nicotine", "etOH"),
                      "vitA"=c("treat", "vitA", "etOH"),
                      "vitD"=c("treat", "vitD", "etOH"),
                      "vitE"=c("treat", "vitE", "etOH"),
                      "water"=c("treat", "water", "etOH"),
                      "zinc"=c("treat", "zinc", "etOH"))


contrast.list <- list("caffeine"=c("new_treat", "caffeine", "control"),
                      "nicotine"=c("new_treat", "nicotine", "control"),
                      "vitA"=c("new_treat", "vitA", "control"),
                      "vitD"=c("new_treat", "vitD", "control"),
                      "vitE"=c("new_treat", "vitE", "control"),
                      "zinc"=c("new_treat", "zinc", "control"))


cat("2.1", "run DESeq batch separately", "\n")


#MCls <- c("CD4Naive", "TCM", "NKcell", "Bcell", "CD8Naive", "TEM", "Monocyte", "CTL", "dnT")

MCls <- unique(cvt.2$MCls)
MCls


## #----INDIVIDUAL STEP----
## head(YtX)
## head(cvt.2)

## MCls <- "0_CD4Naive"

## treatment=c("caffeine"="red", "nicotine"="tan",
##             "vitA"="tan4", "vitD"="seagreen4",
##             "vitE"="salmon3", "water"="grey", "zinc"="maroon3")


## cvt0 <- cvt.3%>%dplyr::filter(MCls=="0_CD4Naive")
## #head(cvt0)
## #cvt0$treat <- relevel(factor(cvt0$treat), ref = "water")

## YtX0 <- YtX[,cvt0$bti]
## #head(YtX0)

## ##
## dds <- DESeqDataSetFromMatrix(YtX0, cvt0, ~sampleID+vehicle+new_treat) #also controlling individual effect

## #dds <- DESeqDataSetFromMatrix(YtX0, cvt0, ~treat) #also controlling individual effect

## ##head(cvt.trial.2)
## ##head(YtX_sel.3)
## ##dds <- DESeqDataSetFromMatrix(YtX_sel.3, cvt.trial.2, ~V4+V3) #also controlling individual effect


## vd <- vst(dds)

## treatment=c("caffeine"="red", "nicotine"="tan",
##             "vitA"="tan4", "vitD"="seagreen4",
##             "vitE"="salmon3", "water"="grey",
##             "zinc"="maroon3", "etOH"="grey")

## cols <- c("red", "grey", "tan", "tan4", "seagreen4", "salmon3", "grey", "maroon3")

## color_key <- c("red", "tan", "tan4", "seagreen4", "salmon3", "grey", "maroon3")
## names(color_key) <- c("caffeine", "nicotine", "vitA", "vitD", "vitE", "water", "zinc")

## data <- plotPCA(vd, intgroup=c("new_treat"), returnData=TRUE)
## head(data)

## data <- data %>% mutate(sampleID=str_split(data$name, "_", simplify=T)[,4])
## head(data)

## #p <- qplot(PC1,PC2,color=treat,data=data,size=I(10))

## #p <- ggplot(data, aes(x=PC1, y=PC2)) + geom_point(aes(colour = factor(treat)))
## #p + scale_colour_manual(values = treatment)
## #p + scale_fill_manual(values = color_key[names(color_key) %in% data$treat])

## p <- ggplot(data, aes(x=PC1, y=PC2, colour = factor(treat))) +
##     geom_point() +
##     scale_colour_manual(values = treatment)+
##     ggtitle("")
## ggsave("./1_DiffRNA_2.outs/PCA_condition.png",p)

## png(paste0("./1_DiffRNA_2.outs/plotpca_treatandind_filt_", reso, "_", oneX, "_cn.png"), width=2000, height=2000, pointsize=14, res=175)
## print(p)
## dev.off()



## dds <- DESeq(dds)

## colnames(coef(dds))

## head(coef(dds))

## head(results(dds))
## length(results(dds))

## res <- contrast.list%>%map(~results(dds, contrast=.x))
## head(res)

## res2 <- tibble(contrast=names(res),MCls=oneX, data=map(res,tidy))%>%unnest(data)
## ##
## time1 <- Sys.time()
## elapsed <- difftime(time1, time0, units="mins")
## cat(oneX, elapsed, "Done\n")
## res2





## #----data test----
## head(YtX)

## #relevel(cvt.3$vehicle)
## table(cvt.3$vehicle, cvt.3$treat)
## table(cvt.3$vehicle, cvt.3$new_treat)
## table(cvt.3$treat, cvt.3$new_treat)

## head(cvt.3)




#----FUNCTION----
res.DE <- map_dfr(MCls, function(oneX){
  time0 <- Sys.time()
###
  ##cvt0 <- cvt.2%>%dplyr::filter(MCls==oneX)
  cvt0 <- cvt.3%>%dplyr::filter(MCls==oneX)
  YtX0 <- YtX[,cvt0$bti]
##
  ##dds <- DESeqDataSetFromMatrix(YtX0, cvt0, ~ sampleID + vehicle + new_treat) #also controlling individual effect
  ##dds <- DESeqDataSetFromMatrix(YtX0, cvt0, ~ sampleID + treat) #also controlling individual effect
  dds <- DESeqDataSetFromMatrix(YtX0, cvt0, ~ sampleID + new_treat) #also controlling individual effect
  ##
  ##
  ##
  ##
  ##vd <- vst(dds)
  ###p <- plotPCA(vd, intgroup=c("treat"))
  ###png(paste0("./1_DiffRNA_2.outs/plotpca_treatandind_filt_", reso, "_", oneX, "_cn.png"), width=2000, height=2000, pointsize=14, res=175)
  ###print(p)
  ###dev.off()
  ##
  ##
  ##correlations - but not in loop
  ##
  ##
  ##data <- plotPCA(vd, intgroup=c("new_treat"), returnData=TRUE)
  ##data <- data %>% mutate(sampleID=str_split(data$name, "_", simplify=T)[,4])
  ###p <- qplot(PC1,PC2,color=treatment,data=data,size=I(10))
  ##p <- ggplot(data, aes(x=PC1, y=PC2, colour = factor(sampleID))) +
   ##geom_point() +
      ###scale_colour_manual(values = treatment)+
      ##ggtitle(oneX)
  ##ggsave(paste0("./1_DiffRNA_2.outs/plotpca_treatandind_filt_", reso, "_", oneX, "_cn_3_include_vehicle.png"),p)
  ##
  ##
  ##
  ##
  dds <- DESeq(dds)
  res <- contrast.list%>%map(~results(dds, contrast=.x))
  res2 <- tibble(contrast=names(res),MCls=oneX, data=map(res,tidy))%>%unnest(data)
##
  time1 <- Sys.time()
  elapsed <- difftime(time1, time0, units="mins")
  cat(oneX, sum(res2$p.adjusted<0.1), elapsed, "Done\n")
  res2
})

head(res.DE)
table(res.DE$contrast)
table(res.DE$MCls)


opfn <- paste0("./1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_", reso, "_cn.rds") #~sampleID+treat

opfn <- paste0("./1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_", reso, "_cn_sampleID+treat_allvsetOH.rds")

opfn <- paste0("./1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_", reso, "_cn_include_vehicle.rds")

opfn <- paste0("./1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_", reso, "_cn_sampleID+new_treat.rds")

write_rds(res.DE, opfn)

##res.DE <- readRDS("./1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn.rds")

##res.DE.contrastetOH <- readRDS("./1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn_sampleID+treat_allvsetOH.rds")

res.DE <- readRDS("./1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn_sampleID+new_treat.rds")

head(res.DE)
#table(res.DE$MCls)
#nrow(res.DE)
#nrow(res.DE %>% drop_na(p.adjusted))
#nrow(res.DE %>% drop_na(p.value))

#res.DE.vehicle <- readRDS("./1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn_include_vehicle.rds")
#res.DE.vehicle <- res.DE
#head(res.DE.vehicle)
#table(res.DE.vehicle$contrast)





## ######################################################
## # significant genes
## ######################################################
res <- res.DE
#res.etoh <- res.DE.contrastetOH
padj_cutoff <- 0.1




#------------------------------
sig_res <- dplyr::filter(res, p.adjusted < padj_cutoff, abs(estimate) > 0.5) %>% dplyr::arrange(p.adjusted)

#write.csv(sig_res, paste0("./1_DiffRNA_2.outs/sig_genes.csv"), quote = FALSE, row.names = FALSE)

#sig_res <- read.csv("./1_DiffRNA_2.outs/sig_genes.csv")

head(sig_res)
table(sig_res$MCls, sig_res$contrast)





#------------------------------
##sig_res.etoh <- dplyr::filter(res.etoh, p.adjusted < padj_cutoff, abs(estimate) >0.5) %>% dplyr::arrange(p.adjusted)
##head(sig_res.etoh)
##table(sig_res.etoh$MCls, sig_res.etoh$contrast)




#------------------------------
sig_res_split <- split(sig_res, f = list(sig_res$contrast,sig_res$MCls))

head(sig_res_split)
#sig_res_split[["vitD.TEM"]]

length(sig_res_split)

names(sig_res_split)

sig_res_split2 = list()

for (i in names(sig_res_split)){
    print(i)
    sig_res_split2[[i]] <- sig_res_split[[i]]  %>% dplyr::arrange(p.adjusted) %>% dplyr::pull(gene) %>% head(n=20)
}

sig_res_split2



#######################################################
# table of up and down regulated genes #1
#######################################################
#------------------------------------------------------
# Nicole script
#------------------------------------------------------
res <- res.DE

FDR <- 0.1
pos_foldchange <- 0.5
neg_foldchange <- -0.5


#------all res------
## res.DE.0.05 <- readRDS("./1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.05_cn.rds")
## res.DE.0.1 <- readRDS("./1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn.rds")
## res.DE.0.15 <- readRDS("./1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.15_cn.rds")
## res.DE.0.2 <- readRDS("./1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.2_cn.rds")


## ## res.DE.0.05.2 <- res.DE.0.05 %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=FDR)  %>% dplyr::filter(abs(estimate) >= pos_foldchange) %>% dplyr::select(contrast, MCls, gene) %>% unique()
## ## res.DE.0.1.2 <- res.DE.0.1 %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=FDR)  %>% dplyr::filter(abs(estimate) >= pos_foldchange) %>% dplyr::select(contrast, MCls,  gene) %>% unique()
## ## res.DE.0.15.2 <- res.DE.0.15 %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=FDR)  %>% dplyr::filter(abs(estimate) >= pos_foldchange) %>% dplyr::select(contrast, MCls,  gene) %>% unique()
## ## res.DE.0.2.2 <- res.DE.0.2 %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=FDR)  %>% dplyr::filter(abs(estimate) >= pos_foldchange) %>% dplyr::select(contrast, MCls,  gene) %>% unique()


## res.DE.0.05.2 <- res.DE.0.05 %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=FDR)  %>% dplyr::filter(abs(estimate) >= pos_foldchange) %>% dplyr::select(contrast, gene) %>% unique()
## res.DE.0.1.2 <- res.DE.0.1 %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=FDR)  %>% dplyr::filter(abs(estimate) >= pos_foldchange) %>% dplyr::select(contrast, gene) %>% unique()
## res.DE.0.15.2 <- res.DE.0.15 %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=FDR)  %>% dplyr::filter(abs(estimate) >= pos_foldchange) %>% dplyr::select(contrast, gene) %>% unique()
## res.DE.0.2.2 <- res.DE.0.2 %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=FDR)  %>% dplyr::filter(abs(estimate) >= pos_foldchange) %>% dplyr::select(contrast, gene) %>% unique()


## a <-as.data.frame(table(res.DE.0.05.2$contrast))
## colnames(a) <- c('treats','0.05')
## #a
## b <-as.data.frame(table(res.DE.0.1.2$contrast))
## colnames(b) <- c('treats','0.1')
## #b
## c <- as.data.frame(table(res.DE.0.15.2$contrast))
## colnames(c) <- c('treats','0.15')
## #c
## #d <- merge(x = a, y = c(b,c), by = "Var1")
## d <- as.data.frame(table(res.DE.0.2.2$contrast))
## colnames(d) <- c('treats','0.2')
## #d
## e <- left_join(a,b, by = "treats") %>% left_join(., c, by='treats') %>% left_join(., d, by='treats')
## e
## ## colnames(res.DE.0.05.2)
#------all res------

## nrow(res %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=FDR)  %>% dplyr::filter(estimate >=pos_foldchange))
## nrow(res %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=FDR) %>% dplyr::filter(estimate <=neg_foldchange))

res[!is.na(res$p.adjusted) & res$p.adjusted<=FDR & res$estimate>=pos_foldchange, "Significance"] <- "Upregulated"
res[!is.na(res$p.adjusted) & res$p.adjusted<=FDR & res$estimate<=neg_foldchange, "Significance"] <- "Downregulated"
res[!is.na(res$p.adjusted) & res$p.adjusted>FDR,"Significance"] <- "Not Significant"
res[is.na(res$p.adjusted),"Significance"] <- "NA"
#res2 <- res %>% dplyr::filter(Significance %in% c("Upregulated", "Downregulated"))
#res




#----barplot----
caff_up <- res %>% filter(res$contrast=="caffeine" & res$Significance=="Upregulated")
nic_up <- res %>% filter(res$contrast=="nicotine" & res$Significance=="Upregulated
")
zn_up <- res %>% filter(res$contrast=="zinc" & res$Significance=="Upregulated")
vitA_up <- res %>% filter(res$contrast=="vitA" & res$Significance=="Upregulated")
vitD_up <- res %>% filter(res$contrast=="vitD" & res$Significance=="Upregulated")
vitE_up <- res %>% filter(res$contrast=="vitE" & res$Significance=="Upregulated")



caff_down <- res %>% filter(res$contrast=="caffeine" & res$Significance=="Downregulated")
nic_down <- res %>% filter(res$contrast=="nicotine" & res$Significance=="Downregulated")
zn_down <- res %>% filter(res$contrast=="zinc" & res$Significance=="Downregulated")
vitA_down <- res %>% filter(res$contrast=="vitA" & res$Significance=="Downregulated")
vitD_down <- res %>% filter(res$contrast=="vitD" & res$Significance=="Downregulated")
vitE_down <- res %>% filter(res$contrast=="vitE" & res$Significance=="Downregulated")


Treatment <- c("Caffeine", "Nicotine", "Zinc", "VitaminA", "VitaminD", "VitaminE",
               "Caffeine", "Nicotine", "Zinc", "VitaminA", "VitaminD", "VitaminE")

DEGs <- c(nrow(caff_up), nrow(nic_up),nrow(zn_up),
          nrow(vitA_up),nrow(vitD_up),nrow(vitE_up),
          nrow(caff_down), nrow(nic_down),nrow(zn_down),
          nrow(vitA_down),nrow(vitD_down),nrow(vitE_down))

df <- as.data.frame(Treatment)
df$DEGs <- DEGs
Direction <- c("Up", "Up","Up","Up", "Up","Up",
               "Down","Down","Down","Down","Down","Down")
df$Direction <- Direction

df

png("./1_DiffRNA_2.outs/RNA_barplot_control.png", width=1000, height=1000, pointsize=18, res=210)
p <- ggplot(df, aes(x=Treatment, y=DEGs, alpha=Direction, color=Treatment, fill=Treatment)) +
    ##  theme(legend.title= element_blank(), axis.title.x = element_text(size = 34), axis.title.y = element_text(size = 34), legend.text = element_text(size = 34), title = element_text(size = 34))+
    geom_bar(stat="identity") +
    scale_alpha_discrete(range=c(.5,1)) +
    theme_bw() +
    ##scale_color_manual(values=c("#00BA38","#FF6699","#C77CFF")) +
    ##scale_fill_manual(values=c("#00BA38","#FF6699","#C77CFF"))+
    scale_color_manual(values=c("red","tan", "tan4", "seagreen4", "salmon3", "maroon3")) +
    scale_fill_manual(values=c("red","tan", "tan4", "seagreen4", "salmon3","maroon3"))+
    theme(text = element_text(size = 14))+
    theme(axis.text.x = element_text(angle = 90))
print(p)
dev.off()





#---------------------------
res_sig <- res %>% dplyr::filter(Significance %in% c("Upregulated", "Downregulated"))
nrow(res_sig)

sig <- table(res_sig$MCls, res_sig$contrast)
sig

table(res_sig$Significance)

#write.csv(up, paste0("./1_DiffRNA_2.outs/up_treat&ind_filt_", reso, "_", FDR, "_", pos_foldchange, "_cn.csv"), quote = FALSE, row.names = TRUE)

#up <- read.csv(file = './1_DiffRNA_2.outs/up_treat&ind_filt.csv')
#up


## MCls <- unique(res_up$MCls)
## #MCls
## res_up.unique <- map_dfr(MCls, function(oneX){
##       time0 <- Sys.time()
##         ###
##         res_up.2 <- res_up%>%dplyr::filter(MCls==oneX)
##         ##
##         res_up.3 <- res_up%>%dplyr::filter(MCls!=oneX)
##             #mutate(count=ifelse(nrow(cvt.2)>=8,1,0)) 
##         ##
##         other.genes <- unique(res_up.3$gene)
##         res_up.2.2 <- res_up.2%>%dplyr::filter(!gene %in% other.genes)
##         time1 <- Sys.time()
##         elapsed <- difftime(time1, time0, units="mins")
##         cat(oneX, elapsed, "Done\n")
##         res_up.2.2
##       })
## #res_up.unique
## #nrow(res_up.unique)

## up.unique <- table(res_up.unique$MCls, res_up.unique$contrast)
## #up

## write.csv(up.unique, paste0("./1_DiffRNA_2.outs/up_treat&ind_filt_unique_", reso, "_", FDR, "_", pos_foldchange, ".csv"), quote = FALSE, row.names = TRUE)

library(reshape2)
sig.matrix <- as.matrix(sig)
sig.melt <- melt(sig.matrix)


sig.p <- ggplot(sig.melt, aes(x = Var2, y = Var1)) +
      geom_tile(aes(fill=value), color="black")+
      scale_fill_gradient(low = "white", high = "white")+
      geom_text(aes(label = value), color = "black", size = 5) +
      ggtitle(paste0("Significant", "_", reso, "_", FDR, "-", pos_foldchange)) +
      theme(axis.text.x = element_text(angle = 90),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position="none")

## up.unique.p <- ggplot(up.unique.melt, aes(x = Var2, y = Var1)) +
##     geom_tile(aes(fill=value))+
##     geom_text(aes(label = value), color = "white", size = 4) +
##     ggtitle("Upregulated.unique") +
##     theme(axis.text.x = element_text(angle = 90))


library("ggpubr")

figure <- ggarrange(sig.p +
                    font("x.text", size = 14))
                   # ncol = 2, nrow = 2)

png(paste0("./1_DiffRNA_2.outs/plot_sig_treatandind_filt_all_", reso, "_", FDR, "_", pos_foldchange, "_cn_control.png"), width=500, height=500, pointsize=16, res=125)
print(figure)
dev.off()








#---------------------------
res_up <- res %>% dplyr::filter(Significance %in% "Upregulated")
nrow(res_up)

up <- table(res_up$MCls, res_up$contrast)
#up

write.csv(up, paste0("./1_DiffRNA_2.outs/up_treat&ind_filt_", reso, "_", FDR, "_", pos_foldchange, "_cn.csv"), quote = FALSE, row.names = TRUE)

#up <- read.csv(file = './1_DiffRNA_2.outs/up_treat&ind_filt.csv')
#up


## MCls <- unique(res_up$MCls)
## #MCls
## res_up.unique <- map_dfr(MCls, function(oneX){
##       time0 <- Sys.time()
##         ###
##         res_up.2 <- res_up%>%dplyr::filter(MCls==oneX)
##         ##
##         res_up.3 <- res_up%>%dplyr::filter(MCls!=oneX)
##             #mutate(count=ifelse(nrow(cvt.2)>=8,1,0)) 
##         ##
##         other.genes <- unique(res_up.3$gene)
##         res_up.2.2 <- res_up.2%>%dplyr::filter(!gene %in% other.genes)
##         time1 <- Sys.time()
##         elapsed <- difftime(time1, time0, units="mins")
##         cat(oneX, elapsed, "Done\n")
##         res_up.2.2
##       })
## #res_up.unique
## #nrow(res_up.unique)

## up.unique <- table(res_up.unique$MCls, res_up.unique$contrast)
## #up

## write.csv(up.unique, paste0("./1_DiffRNA_2.outs/up_treat&ind_filt_unique_", reso, "_", FDR, "_", pos_foldchange, ".csv"), quote = FALSE, row.names = TRUE)



#---------------------------
res_down <- res %>% dplyr::filter(Significance %in% "Downregulated")

down <- table(res_down$MCls, res_down$contrast)
#down

write.csv(down, paste0("./1_DiffRNA_2.outs/down_treat&ind_filt_", reso, "_", FDR, "_", pos_foldchange, "_cn.csv"), quote = FALSE, row.names = TRUE)

#down <- read.csv(file = './1_DiffRNA_2.outs/down_treat&ind_filt.csv')
#down


## MCls <- unique(res_down$MCls)
## #MCls
## res_down.unique <- map_dfr(MCls, function(oneX){
##       time0 <- Sys.time()
##         ###
##         res_down.2 <- res_down%>%dplyr::filter(MCls==oneX)
##         ##
##         res_down.3 <- res_down%>%dplyr::filter(MCls!=oneX)
##             #mutate(count=ifelse(nrow(cvt.2)>=8,1,0)) 
##         ##
##         other.genes <- unique(res_down.3$gene)
##         res_down.2.2 <- res_down.2%>%dplyr::filter(!gene %in% other.genes)
##         time1 <- Sys.time()
##         elapsed <- difftime(time1, time0, units="mins")
##         cat(oneX, elapsed, "Done\n")
##         res_down.2.2
##       })
## #res_down.unique
## #nrow(res_down.unique)


## down.unique <- table(res_down.unique$MCls, res_down.unique$contrast)
## #down
## write.csv(down.unique, paste0("./1_DiffRNA_2.outs/down_treat&ind_filt_unique_", reso, "_", FDR, "_", pos_foldchange, ".csv"), quote = FALSE, row.names = TRUE)


# -----------number of cells------------------------------------------------
#cvt.2

dd3 <- left_join(cvt.2, dd2)
#head(dd3)

#dd3 <- dd3 %>% dplyr::filter(dd3$MCls %in% res$MCls, !dd3$treat %in% "etOH")
dd3 <- dd3 %>% dplyr::filter(dd3$MCls %in% res$MCls)
#head(dd3)
#nrow(dd3)
#sum(dd3$ncell)

dd4 <- dd3 %>% group_by(treat, MCls) %>% summarise(ncell_count = sum(ncell))
#dd4

ncell_table <- data.frame(pivot_wider(dd4, names_from=treat, values_from=ncell_count))

rownames(ncell_table) <- ncell_table$MCls
ncell_table$MCls <- NULL
ncell_table

#write.csv(ncell_table, paste0("./1_DiffRNA_2.outs/ncell_treat$ind_filt_", reso, "_", FDR, "_", pos_foldchange, "_cn.csv"), quote = FALSE, row.names = TRUE)
write.csv(ncell_table, paste0("./1_DiffRNA_2.outs/ncell_treat$ind_filt_cn.csv"), quote = FALSE, row.names = TRUE)


#ncell_table <- read.csv("./1_DiffRNA_2.outs/ncell_treat$ind_filt.csv")


#########################################################
# matrix_plots
#########################################################
library(reshape2)

up.matrix <-as.matrix(up) # create a numeric matrix object
up.melt <- melt(up.matrix)

#up.unique.matrix <-as.matrix(up.unique) # create a numeric matrix object
#up.unique.melt <- melt(up.unique.matrix)

down.matrix <- as.matrix(down)
down.melt <- melt(down.matrix)

#down.unique.matrix <- as.matrix(down.unique)
#down.unique.melt <- melt(down.unique.matrix)


#ncell
#rownames(ncell)<- ncell$X
#ncell$X <- NULL
ncell.matrix <- as.matrix(ncell_table)
ncell.melt <- melt(ncell.matrix)

#comb
#rownames(comb)<- comb$X
#comb$X <- NULL
#comb.2 <- subset(comb, select = -c(etOH))
#comb.2 <- comb.2 %>% dplyr::filter(rownames(comb.3) %in% rownames(up))
comb.3.matrix <- as.matrix(comb.2)
comb.3.matrix <- subset.matrix(comb.3.matrix, select = -c(etOH))
comb.3.melt <- melt(comb.3.matrix)

up.p <- ggplot(up.melt, aes(x = Var2, y = Var1)) +
    geom_tile(aes(fill=value), color="black")+
    scale_fill_gradient(low = "white", high = "white")+
    geom_text(aes(label = value), color = "black", size = 5) +
    ggtitle(paste0("Upregulated", "_", reso, "_", FDR, "-", pos_foldchange)) +
    theme(axis.text.x = element_text(angle = 90),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position="none")



down.p <- ggplot(down.melt, aes(x = Var2, y = Var1)) +
      geom_tile(aes(fill=value), color="black")+
      scale_fill_gradient(low = "white", high = "white")+
      geom_text(aes(label = value), color = "black", size = 5) +
      ggtitle(paste0("Downregulated", "_", reso, "_", FDR, "-", pos_foldchange)) +
      theme(axis.text.x = element_text(angle = 90),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position="none")

## up.unique.p <- ggplot(up.unique.melt, aes(x = Var2, y = Var1)) +
##     geom_tile(aes(fill=value))+
##     geom_text(aes(label = value), color = "white", size = 4) +
##     ggtitle("Upregulated.unique") +
##     theme(axis.text.x = element_text(angle = 90))

## down.unique.p <- ggplot(down.unique.melt, aes(x = Var2, y = Var1)) +
##       geom_tile(aes(fill=value))+
##       geom_text(aes(label = value), color = "white", size = 4) +
##       ggtitle("Downregulated.unique") +
##       theme(axis.text.x = element_text(angle = 90))


ncell.p <- ggplot(ncell.melt, aes(x = Var2, y = Var1)) +
      geom_tile(aes(fill=value), color="black")+
      scale_fill_gradient(low = "white", high = "white")+
      geom_text(aes(label = value), color = "black", size = 5) +
      ggtitle(paste0("ncells", "_", reso, "_", FDR, "-", pos_foldchange)) +
      theme(axis.text.x = element_text(angle = 90),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position="none")


comb.3.p <- ggplot(comb.3.melt, aes(x = Var2, y = Var1)) +
      geom_tile(aes(fill=value), color="black")+
      scale_fill_gradient(low = "white", high = "white")+
      geom_text(aes(label = value), color = "black", size = 5) +
      ggtitle(paste0("Combinations", "_", reso, "_", FDR, "-", pos_foldchange)) +
      theme(axis.text.x = element_text(angle = 90),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position="none")


## ncell.p <- ggplot(ncell.melt, aes(x = Var2, y = Var1)) +
##       geom_tile(aes(fill=value))+
##       geom_text(aes(label = value), color = "white", size = 4) +
##       ggtitle("ncells") +
##       theme(axis.text.x = element_text(angle = 90))

## comb.3.p <- ggplot(comb.3.melt, aes(x = Var2, y = Var1)) +
##       geom_tile(aes(fill=value))+
##       geom_text(aes(label = value), color = "white", size = 4) +
##       ggtitle("Combinations") +
##       theme(axis.text.x = element_text(angle = 90))

library("ggpubr")

figure <- ggarrange(comb.3.p, ncell.p, up.p, down.p +
                    font("x.text", size = 14),
                    ncol = 2, nrow = 2)

png(paste0("./1_DiffRNA_2.outs/plot_all_treatandind_filt_all_", reso, "_", FDR, "_", pos_foldchange, "_cn.png"), width=2000, height=2000, pointsize=16, res=175)
print(figure)
dev.off()


##############################################################
###
### 1. MA plots
##############################################################
table(res$MCls)
table(res$contrast)

figfn <- paste0("./1_DiffRNA_2.outs/Figure1.1_MA_clusters_treat&ind_filt_", reso, "_", FDR, "_", pos_foldchange, "_cn_control.png")
png(figfn, width=6000, height=8000, pointsize=33, res=250)
par(mar=c(4,4,2,2),mgp=c(2,1,0))
x <- matrix(1:48, 8, 6, byrow=T)
layout(x)
#fn <- "./1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered.rds"
#res <- read_rds(fn)%>%mutate(color=ifelse((p.adjusted<0.1)&(!is.na(p.adjusted)), T, F))
res.1 <- res%>%mutate(color=ifelse((p.adjusted<=FDR)&(!is.na(p.adjusted)), T, F))
# MCls <- c("CD4Naive", "TCM", "NKcell", "Bcell", "CD8Naive", "TEM", "Monocyte", "CTL", "dnT")
#MCls <- unique(cvt.2$MCls)
MCls <- unique(res$MCls)
for (oneMCl in MCls){
##1
   res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="caffeine")%>%
           dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
   print(plotMA(res.2[,1:3], colLine="NA", main="caffeine", cex.main=1, cex.axis=0.8, cex.lab=1))
###2    
   res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="nicotine")%>%
          dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
   print(plotMA(res.2[,1:3], colLine="NA", main="nicotine", cex.main=1, cex.axis=0.8, cex.lab=1))
###3
   res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="vitA")%>%
           dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
   print(plotMA(res.2[,1:3], colLine="NA", main="vitA", cex.main=1, cex.axis=0.8, cex.lab=1))
###4
   res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="vitD")%>%
           dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
   print(plotMA(res.2[,1:3], colLine="NA", main="vitD", cex.main=1, cex.axis=0.8, cex.lab=1))
###5
    res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="vitE")%>%
           dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
   print(plotMA(res.2[,1:3], colLine="NA", main="vitE", cex.main=1, cex.axis=0.8, cex.lab=1))
###6    
   ##  res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="water")%>%
   ##         dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
   ## print(plotMA(res.2[,1:3], colLine="NA", main="water", cex.main=1, cex.axis=0.8, cex.lab=1))
###7    
    res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="zinc")%>%
           dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
   print(plotMA(res.2[,1:3], colLine="NA", main="zinc", cex.main=1, cex.axis=0.8, cex.lab=1))
   print(mtext(oneMCl, side=4, line=0.5, cex=1, col="blue")) 
}
dev.off()




### 2, qq plots
figfn <- paste0("./1_DiffRNA_2.outs/Figure1.2_qq_clusters_treat&ind_filt_", reso, "_", FDR, "_", pos_foldchange, "_cn_control.png")
png(figfn, width=3500, height=4500, pointsize=18, res=225)
par(mar=c(4,4,2,2),mgp=c(2,1,0))
x <- matrix(1:48, 8, 6, byrow=T)
layout(x)
#fn <- "./1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered.rds"
#res <- read_rds(fn)%>%drop_na(p.value)
res.1 <- res%>%drop_na(p.value)
# MCls <- c("CD4Naive", "TCM", "NKcell", "Bcell", "CD8Naive", "TEM", "Monocyte", "CTL", "dnT")
#MCls <- unique(cvt.2$MCls)
MCls <- unique(res$MCls)
for (oneMCl in MCls){
   ##1
   res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="caffeine")
   print( qq(res.2$p.value, main="caffeine", cex.main=1, cex.axis=0.8, cex.lab=1))
   ##2
   res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="nicotine")
   print( qq(res.2$p.value, main="nicotine", cex.main=1, cex.axis=0.8, cex.lab=1))
   ##3 
   res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="vitA")
   print( qq(res.2$p.value, main="vitA", cex.main=1, cex.axis=0.8, cex.lab=1))
   ##4
   res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="vitD")
   print( qq(res.2$p.value, main="vitD", cex.main=1, cex.axis=0.8, cex.lab=1))
   ##5
   res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="vitE")
   print( qq(res.2$p.value, main="vitE", cex.main=1, cex.axis=0.8, cex.lab=1))
   ##6
   ##res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="water")
   ##print( qq(res.2$p.value, main="water", cex.main=1, cex.axis=0.8, cex.lab=1))
   ##7
   res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="zinc")
   print( qq(res.2$p.value, main="zinc", cex.main=1, cex.axis=0.8, cex.lab=1))   
   print( mtext(oneMCl, side=4, line=0.5, cex=1, col="blue"))  
}
dev.off()


###
### 3, canno plots


#------------------plot_loop--------------------------------
#nrow(res) #if qq plot has been called before this step, this already is with dropped NA values
res.1 <- res %>% drop_na(p.value)
max.y <- floor(max(-log10(res.1$p.value)))/2
#max.y
png(paste0("./1_DiffRNA_2.outs/Figure1.3_qq_clusters_together_treatandind_filt_", reso,"_", FDR, "_", pos_foldchange, "_cn_control.png"), width=2000, height=4000, pointsize=22, res=200)
par(mar=c(4,4,2,2),mgp=c(2,1,0))
x <- matrix(1:8, 4, 2, byrow=T)
layout(x)
# MCls <- c("CD4Naive", "TCM", "NKcell", "Bcell", "CD8Naive", "TEM", "Monocyte", "CTL", "dnT")
#MCls <- unique(cvt.2$MCls)
MCls <- unique(res$MCls)
for (oneMCl in MCls){
    ##1
    res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="caffeine")%>% drop_na(p.value)
    plot(-log10(1:length(res.2$p.value)/length(res.2$p.value)),
         -log10(sort(res.2$p.value)),
         main=oneMCl,
         pch = 19,
         ylim = c(0, max.y),
         ylab="observed -log10(p)",
         xlab="expected -log10(p)",
         col="red")
    ##2
    res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="nicotine")%>% drop_na(p.value)
    points(-log10(1:length(res.2$p.value)/length(res.2$p.value)),-log10(sort(res.2$p.value)),  col="tan", pch = 19)   
    ##3 
    res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="zinc")%>% drop_na(p.value)
    points(-log10(1:length(res.2$p.value)/length(res.2$p.value)),-log10(sort(res.2$p.value)),  col="maroon3", pch = 19)
    ##4
    res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="vitA")%>% drop_na(p.value)
    points(-log10(1:length(res.2$p.value)/length(res.2$p.value)),-log10(sort(res.2$p.value)),  col="tan4", pch = 19)
    ##5
    res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="vitD")%>% drop_na(p.value)
    points(-log10(1:length(res.2$p.value)/length(res.2$p.value)),-log10(sort(res.2$p.value)),  col="seagreen4", pch = 19)
    ##6
    ##res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="water")%>% drop_na(p.value)
    ##points(-log10(1:length(res.2$p.value)/length(res.2$p.value)),-log10(sort(res.2$p.value)),  col="grey", pch = 19)
    ##7
    res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="vitE")%>% drop_na(p.value)
    points(-log10(1:length(res.2$p.value)/length(res.2$p.value)),-log10(sort(res.2$p.value)),  col="salmon3", pch = 19)
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







#######################################################
# table of up and down regulated genes #2
#######################################################
#------------------------------------------------------
# Nicole script
#------------------------------------------------------
res <- res.DE

FDR <- 0.05
pos_foldchange <- 0.00
neg_foldchange <- 0.00

#------all res------
## ## res.DE.0.05 <- readRDS("./1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.05.rds")
## ## res.DE.0.1 <- readRDS("./1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1.rds")
## ## res.DE.0.15 <- readRDS("./1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.15.rds")

## ## res.DE.0.05.2 <- res.DE.0.05 %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=FDR)  %>% dplyr::filter(abs(estimate) >= pos_foldchange) %>% dplyr::select(contrast, MCls, gene) %>% unique()
## ## res.DE.0.1.2 <- res.DE.0.1 %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=FDR)  %>% dplyr::filter(abs(estimate) >= pos_foldchange) %>% dplyr::select(contrast, MCls,  gene) %>% unique()
## ## res.DE.0.15.2 <- res.DE.0.15 %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=FDR)  %>% dplyr::filter(abs(estimate) >= pos_foldchange) %>% dplyr::select(contrast, MCls,  gene) %>% unique()

## res.DE.0.05.2 <- res.DE.0.05 %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=FDR)  %>% dplyr::filter(abs(estimate) >= pos_foldchange) %>% dplyr::select(contrast, gene) %>% unique()
## res.DE.0.1.2 <- res.DE.0.1 %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=FDR)  %>% dplyr::filter(abs(estimate) >= pos_foldchange) %>% dplyr::select(contrast, gene) %>% unique()
## res.DE.0.15.2 <- res.DE.0.15 %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=FDR)  %>% dplyr::filter(abs(estimate) >= pos_foldchange) %>% dplyr::select(contrast, gene) %>% unique()
## res.DE.0.2.2 <- res.DE.0.2 %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=FDR)  %>% dplyr::filter(abs(estimate) >= pos_foldchange) %>% dplyr::select(contrast, gene) %>% unique()



## a <-as.data.frame(table(res.DE.0.05.2$contrast))
## colnames(a) <- c('treats','0.05')
## #a
## b <-as.data.frame(table(res.DE.0.1.2$contrast))
## colnames(b) <- c('treats','0.1')
## #b
## c <- as.data.frame(table(res.DE.0.15.2$contrast))
## colnames(c) <- c('treats','0.15')
## #c
## d <- as.data.frame(table(res.DE.0.2.2$contrast))
## colnames(d) <- c('treats','0.2')
## #c
## #d <- merge(x = a, y = c(b,c), by = "Var1")
## e <- left_join(a,b, by = "treats") %>% left_join(., c, by='treats') %>% left_join(., d, by='treats')
## e
## #colnames(res.DE.0.05.2)
#------all res------

## nrow(res %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=FDR)  %>% dplyr::filter(estimate >=pos_foldchange))
## nrow(res %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=FDR) %>% dplyr::filter(estimate <=neg_foldchange))

res[!is.na(res$p.adjusted) & res$p.adjusted<=FDR & res$estimate>=pos_foldchange, "Significance"] <- "Upregulated"
res[!is.na(res$p.adjusted) & res$p.adjusted<=FDR & res$estimate<=neg_foldchange, "Significance"] <- "Downregulated"
res[!is.na(res$p.adjusted) & res$p.adjusted>FDR,"Significance"] <- "Not Significant"
res[is.na(res$p.adjusted),"Significance"] <- "NA"
#res2 <- res %>% dplyr::filter(Significance %in% c("Upregulated", "Downregulated"))
#res

#---------------------------
res_up <- res %>% dplyr::filter(Significance %in% "Upregulated")
#nrow(res_up)

up <- table(res_up$MCls, res_up$contrast)
#up

write.csv(up, paste0("./1_DiffRNA_2.outs/up_treat&ind_filt_", reso, "_", FDR, "_", pos_foldchange, "_cn.csv"), quote = FALSE, row.names = TRUE)

#up <- read.csv(file = './1_DiffRNA_2.outs/up_treat&ind_filt.csv')
#up


## MCls <- unique(res_up$MCls)
## #MCls
## res_up.unique <- map_dfr(MCls, function(oneX){
##       time0 <- Sys.time()
##         ###
##         res_up.2 <- res_up%>%dplyr::filter(MCls==oneX)
##         ##
##         res_up.3 <- res_up%>%dplyr::filter(MCls!=oneX)
##             #mutate(count=ifelse(nrow(cvt.2)>=8,1,0)) 
##         ##
##         other.genes <- unique(res_up.3$gene)
##         res_up.2.2 <- res_up.2%>%dplyr::filter(!gene %in% other.genes)
##         time1 <- Sys.time()
##         elapsed <- difftime(time1, time0, units="mins")
##         cat(oneX, elapsed, "Done\n")
##         res_up.2.2
##       })
## #res_up.unique
## #nrow(res_up.unique)

## up.unique <- table(res_up.unique$MCls, res_up.unique$contrast)
## #up

## write.csv(up.unique, paste0("./1_DiffRNA_2.outs/up_treat&ind_filt_unique_", reso, "_", FDR, "_", pos_foldchange, ".csv"), quote = FALSE, row.names = TRUE)



#---------------------------
res_down <- res %>% dplyr::filter(Significance %in% "Downregulated")

down <- table(res_down$MCls, res_down$contrast)
#down

write.csv(down, paste0("./1_DiffRNA_2.outs/down_treat&ind_filt_", reso, "_", FDR, "_", pos_foldchange, "_cn.csv"), quote = FALSE, row.names = TRUE)

#down <- read.csv(file = './1_DiffRNA_2.outs/down_treat&ind_filt.csv')
#down


## MCls <- unique(res_down$MCls)
## #MCls
## res_down.unique <- map_dfr(MCls, function(oneX){
##       time0 <- Sys.time()
##         ###
##         res_down.2 <- res_down%>%dplyr::filter(MCls==oneX)
##         ##
##         res_down.3 <- res_down%>%dplyr::filter(MCls!=oneX)
##             #mutate(count=ifelse(nrow(cvt.2)>=8,1,0)) 
##         ##
##         other.genes <- unique(res_down.3$gene)
##         res_down.2.2 <- res_down.2%>%dplyr::filter(!gene %in% other.genes)
##         time1 <- Sys.time()
##         elapsed <- difftime(time1, time0, units="mins")
##         cat(oneX, elapsed, "Done\n")
##         res_down.2.2
##       })
## #res_down.unique
## #nrow(res_down.unique)


## down.unique <- table(res_down.unique$MCls, res_down.unique$contrast)
## #down
## write.csv(down.unique, paste0("./1_DiffRNA_2.outs/down_treat&ind_filt_unique_", reso, "_", FDR, "_", pos_foldchange, ".csv"), quote = FALSE, row.names = TRUE)


# -----------number of cells------------------------------------------------
#cvt.2

dd3 <- left_join(cvt.2, dd2)
#head(dd3)

dd3 <- dd3 %>% dplyr::filter(dd3$MCls %in% res$MCls, !dd3$treat %in% "etOH")
#head(dd3)
#nrow(dd3)
#sum(dd3$ncell)

dd4 <- dd3 %>% group_by(treat, MCls) %>% summarise(ncell_count = sum(ncell))
#dd4

ncell_table <- data.frame(pivot_wider(dd4, names_from=treat, values_from=ncell_count))

rownames(ncell_table) <- ncell_table$MCls
ncell_table$MCls <- NULL
#ncell_table

write.csv(ncell_table, paste0("./1_DiffRNA_2.outs/ncell_treat$ind_filt_", reso, "_", FDR, "_", pos_foldchange, "_cn.csv"), quote = FALSE, row.names = TRUE)

#ncell_table <- read.csv("./1_DiffRNA_2.outs/ncell_treat$ind_filt.csv")


#########################################################
# matrix_plots
#########################################################
library(reshape2)

up.matrix <-as.matrix(up) # create a numeric matrix object
up.melt <- melt(up.matrix)

#up.unique.matrix <-as.matrix(up.unique) # create a numeric matrix object
#up.unique.melt <- melt(up.unique.matrix)

down.matrix <- as.matrix(down)
down.melt <- melt(down.matrix)

#down.unique.matrix <- as.matrix(down.unique)
#down.unique.melt <- melt(down.unique.matrix)


#ncell
#rownames(ncell)<- ncell$X
#ncell$X <- NULL
ncell.matrix <- as.matrix(ncell_table)
ncell.melt <- melt(ncell.matrix)

#comb
#rownames(comb)<- comb$X
#comb$X <- NULL
#comb.2 <- subset(comb, select = -c(etOH))
#comb.2 <- comb.2 %>% dplyr::filter(rownames(comb.3) %in% rownames(up))
comb.3.matrix <- as.matrix(comb.2)
comb.3.matrix <- subset.matrix(comb.3.matrix, select = -c(etOH))
comb.3.melt <- melt(comb.3.matrix)

up.p <- ggplot(up.melt, aes(x = Var2, y = Var1)) +
    geom_tile(aes(fill=value), color="black")+
    scale_fill_gradient(low = "white", high = "white")+
    geom_text(aes(label = value), color = "black", size = 5) +
    ggtitle(paste0("Upregulated", "_", reso, "_", FDR, "-", pos_foldchange)) +
    theme(axis.text.x = element_text(angle = 90),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position="none")


down.p <- ggplot(down.melt, aes(x = Var2, y = Var1)) +
      geom_tile(aes(fill=value), color="black")+
      scale_fill_gradient(low = "white", high = "white")+
      geom_text(aes(label = value), color = "black", size = 5) +
      ggtitle(paste0("Downregulated", "_", reso, "_", FDR, "-", pos_foldchange)) +
      theme(axis.text.x = element_text(angle = 90),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position="none")

## up.unique.p <- ggplot(up.unique.melt, aes(x = Var2, y = Var1)) +
##     geom_tile(aes(fill=value))+
##     geom_text(aes(label = value), color = "white", size = 4) +
##     ggtitle("Upregulated.unique") +
##     theme(axis.text.x = element_text(angle = 90))

## down.unique.p <- ggplot(down.unique.melt, aes(x = Var2, y = Var1)) +
##       geom_tile(aes(fill=value))+
##       geom_text(aes(label = value), color = "white", size = 4) +
##       ggtitle("Downregulated.unique") +
##       theme(axis.text.x = element_text(angle = 90))


ncell.p <- ggplot(ncell.melt, aes(x = Var2, y = Var1)) +
      geom_tile(aes(fill=value), color="black")+
      scale_fill_gradient(low = "white", high = "white")+
      geom_text(aes(label = value), color = "black", size = 5) +
      ggtitle(paste0("ncells", "_", reso, "_", FDR, "-", pos_foldchange)) +
      theme(axis.text.x = element_text(angle = 90),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position="none")


comb.3.p <- ggplot(comb.3.melt, aes(x = Var2, y = Var1)) +
      geom_tile(aes(fill=value), color="black")+
      scale_fill_gradient(low = "white", high = "white")+
      geom_text(aes(label = value), color = "black", size = 5) +
      ggtitle(paste0("Combinations", "_", reso, "_", FDR, "-", pos_foldchange)) +
      theme(axis.text.x = element_text(angle = 90),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position="none")


## ncell.p <- ggplot(ncell.melt, aes(x = Var2, y = Var1)) +
##       geom_tile(aes(fill=value))+
##       geom_text(aes(label = value), color = "white", size = 4) +
##       ggtitle("ncells") +
##       theme(axis.text.x = element_text(angle = 90))

## comb.3.p <- ggplot(comb.3.melt, aes(x = Var2, y = Var1)) +
##       geom_tile(aes(fill=value))+
##       geom_text(aes(label = value), color = "white", size = 4) +
##       ggtitle("Combinations") +
##       theme(axis.text.x = element_text(angle = 90))

library("ggpubr")

figure <- ggarrange(comb.3.p, ncell.p, up.p, down.p +
                    font("x.text", size = 14),
                    ncol = 2, nrow = 2)

png(paste0("./1_DiffRNA_2.outs/plot_all_treatandind_filt_all_", reso, "_", FDR, "_", pos_foldchange, "_cn.png"), width=2000, height=2000, pointsize=16, res=175)
print(figure)
dev.off()



##############################################################
###
### 1. MA plots
##############################################################
figfn <- paste0("./1_DiffRNA_2.outs/Figure1.1_MA_clusters_treat&ind_filt_", reso, "_", FDR, "_", pos_foldchange, "_cn.png")
png(figfn, width=3500, height=4500, pointsize=18, res=225)
par(mar=c(4,4,2,2),mgp=c(2,1,0))
x <- matrix(1:63, 9, 7, byrow=T)
layout(x)
#fn <- "./1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered.rds"
#res <- read_rds(fn)%>%mutate(color=ifelse((p.adjusted<0.1)&(!is.na(p.adjusted)), T, F))
res.1 <- res%>%mutate(color=ifelse((p.adjusted<=FDR)&(!is.na(p.adjusted)), T, F))
# MCls <- c("CD4Naive", "TCM", "NKcell", "Bcell", "CD8Naive", "TEM", "Monocyte", "CTL", "dnT")
MCls <- unique(cvt.2$MCls)
for (oneMCl in MCls){
  ##1
  res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="caffeine")%>%
    dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
  print(plotMA(res.2[,1:3], colLine="NA", main="caffeine", cex.main=1, cex.axis=0.8, cex.lab=1))
  ###2    
  res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="nicotine")%>%
    dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
  print(plotMA(res.2[,1:3], colLine="NA", main="nicotine", cex.main=1, cex.axis=0.8, cex.lab=1))
  ###3
  res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="vitA")%>%
    dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
  print(plotMA(res.2[,1:3], colLine="NA", main="vitA", cex.main=1, cex.axis=0.8, cex.lab=1))
  ###4
  res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="vitD")%>%
    dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
  print(plotMA(res.2[,1:3], colLine="NA", main="vitD", cex.main=1, cex.axis=0.8, cex.lab=1))
  ###5
  res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="vitE")%>%
    dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
  print(plotMA(res.2[,1:3], colLine="NA", main="vitE", cex.main=1, cex.axis=0.8, cex.lab=1))
  ###6    
  res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="water")%>%
    dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
  print(plotMA(res.2[,1:3], colLine="NA", main="water", cex.main=1, cex.axis=0.8, cex.lab=1))
  ###7    
  res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="zinc")%>%
    dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
  print(plotMA(res.2[,1:3], colLine="NA", main="zinc", cex.main=1, cex.axis=0.8, cex.lab=1))
  print(mtext(oneMCl, side=4, line=0.5, cex=1, col="blue")) 
}
dev.off()











cat("#################################", "\n")
cat("RNA analysis done", "\n")
cat("#################################", "\n")



##########################################################################################################################################################
##########################################################################################################################################################
# ATAC
##########################################################################################################################################################
##########################################################################################################################################################

####################################
### 1. Generate pseudo-bulk data ###
####################################
# combined <- read_rds("../1_processing/3_Clustering.outs/1_seurat.cluster.combined.rds")
# colnames(combined@meta.data)

x <- combined@assays$ATAC@ranges
#x

xrange <- ranges(x)
#xrange

count <- combined@assays$ATAC@counts
#head(count)

anno <- data.frame(rn=rownames(count), rnz=rowSums(count),
                   chr=as.character(seqnames(x)), start=start(xrange), end=end(xrange))
# head(anno)

autosome <- paste("chr", as.character(1:22), sep="")
# autosome

annoSel <- anno%>%dplyr::filter(rnz>0, chr%in%autosome)
# head(annoSel)

Y <- count[annoSel$rn,]
# head(Y)

# meta <- atac@meta.data%>%mutate(treat=gsub(".*-ATAC-|_.*", "", NEW_BARCODE)) 
meta <- combined@meta.data
meta$treat <- meta$treats
# colnames(meta)

# cell <- c("0"="Tcell", "1"="Tcell", "2"="NKcell", "3"="Bcell",
#           "4"="Tcell", "5"="Tcell", "6"="Monocyte", "7"="Tcell", "8"="Tcell",
#           "9"="Tcell", "10"="Tcell", "11"="DC", "12"="DC")

## cell <- c("0"="CD4Naive", "1"="TCM", "2"="NKcell", "3"="Bcell",
##           "4"="CD8Naive", "5"="TEM", "6"="Monocyte", "7"="CTL", "8"="dnT",
##           "9"="Treg", "10"="Platelet", "11"="DC", "12"="HSPC", "13"="Bcell")

#cell <- c("0"="CD4Naive", "1"="TCM", "2"="NKcell", "3"="Bcell",
#          "4"="CD8Naive", "5"="TEM", "6"="Monocyte", "7"="CTL", "8"="dnT",
#          "9"="Treg", "10"="Platelet", "11"="DC")

cell <- c("0"="CD4Naive", "1"="TCM", "2"="NKcell", "3"="TEM",
          "4"="Bcell", "5"="CD8Naive", "6"="Monocyte", "7"="dnT",
          "8"="MAIT","9"="Platelet", "10"="DC")

#cell <- c("0"="Tcell", "1"="Tcell", "2"="Tcell", "3"="NKcell",
#         "4"="Bcell", "5"="Monocyte", "6"="Platelet", "7"="DC")

meta$MCls <- cell[as.character(meta$wsnn_res.0.1)]
#colnames(meta)

meta <- meta%>%mutate(bti=paste(wsnn_res.0.1, MCls, treat, SNG.BEST.GUESS, sep="_"))%>%dplyr::select(NEW_BARCODE, bti)
#colnames(meta)

dd <- meta%>%group_by(bti)%>%summarise(ncell=n(),.groups="drop")
#dd

#write_rds(dd, "./1_DiffPeak.outs/0_ncell.rds")



##pseudo-bulk peak data
bti <- factor(meta$bti)
#head(bti)

X <- model.matrix(~0+bti)
#head(X)

YtX <- Y %*% X
YtX <- as.matrix(YtX)
#head(YtX)

colnames(YtX) <- gsub("^bti", "", colnames(YtX))
#colnames(YtX)

###rnz>0,chr:1-22, 111,750*400
# opfn <- "./1_DiffPeak.outs/1_YtX.comb.rds" 
# write_rds(YtX, file=opfn)

###rnz>100, ncell>20, 111,746 * 331
#use when starting analysis again
#dd <- read_rds("./1_DiffPeak.outs/0_ncell.rds")%>%filter(ncell>20)
#YtX <- read_rds("./1_DiffPeak.outs/1_YtX.comb.rds")
#head(YtX)

dd <- dd %>% filter(ncell>20)

anno <- data.frame(rn=rownames(YtX), rnz=rowSums(YtX))
#head(anno)

annoSel <- anno%>%filter(rnz>100)
#head(annoSel)

YtX_sel <- YtX[annoSel$rn, dd$bti]
#head(YtX_sel)

nrow(YtX_sel)
ncol(YtX_sel)

opfn <- paste0("./1_DiffPeak.outs/1_YtX.sel_", reso, "_cn.rds")
write_rds(YtX_sel, file=opfn)

#YtX_sel <- read_rds(paste0("./1_DiffPeak.outs/1_YtX.sel_", reso, "_cn.rds"))

YtX_sel <- read_rds("./1_DiffPeak.outs/1_YtX.sel_0.1_cn.rds")


#####################################################
### 2. Differential analysis for CRP, High vs Low ###
#####################################################
# ### function, adjusted p value
# myqval <- function(pval){
#   qval <- pval
#   ii0 <- !is.na(pval)
#   qval[ii0] <- qvalue(pval[ii0],pi0=1)$qvalues
#   qval
# }

#### Read data
#YtX <- read_rds("./1_DiffPeak.outs/1_YtX.sel.rds")

YtX <- YtX_sel
#str(YtX)

bti2 <- colnames(YtX)
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

# filtered
#-------------------------------------------------------------------------
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

write.csv(comb.2, paste0("./1_DiffPeak.outs/combinations_2_", reso, "_cn.csv"), quote = FALSE, row.names = TRUE)

#comb.2 <- read.csv("./1_DiffPeak.outs/combinations_2.csv")
#comb.2

bti3 <- cvt.2$bti 
#bti3

#ncol(YtX)
YtX <- YtX[,bti3]
#ncol(YtX)

dd2 <- dd %>% dplyr::filter(bti %in% bti3)
#dd2
#sum(dd2$ncell)





#################################################################
### 2.1 call DESeq ###
#################################################################
contrast.list <- list("caffeine"=c("treat", "caffeine", "water"),
                      "nicotine"=c("treat", "nicotine", "water"),
                      "vitA"=c("treat", "vitA", "etOH"),
                      "vitD"=c("treat", "vitD", "etOH"),
                      "vitE"=c("treat", "vitE", "etOH"),
                      "water"=c("treat", "water", "etOH"),
                      "zinc"=c("treat", "zinc", "water"))


cat("2.1", "run DESeq batch separately", "\n")

#MCls <- c("Bcell", "Monocyte", "NK", "Tcell", "DC")
#MCls <- c("Bcell", "Monocyte")
#MCls  <- c("CD4Naive", "TCM", "NKcell", "Bcell",
#           "CD8Naive", "TEM", "Monocyte", "CTL", "dnT")
MCls <- unique(cvt.2$MCls)
MCls

reso = 0.1

treatment=c("caffeine"="red", "nicotine"="tan",
            "vitA"="tan4", "vitD"="seagreen4",
            "vitE"="salmon3", "water"="grey",
            "zinc"="maroon3", "etOH"="grey")

###
res.DE.ATAC <- map_dfr(MCls, function(oneX){
  time0 <- Sys.time()
  ###
  cvt0 <- cvt.2%>%dplyr::filter(MCls==oneX)
  YtX0 <- YtX[,cvt0$bti]
  ##
  dds <- DESeqDataSetFromMatrix(YtX0, cvt0, ~ sampleID + treat)
  vd <- vst(dds)
  #p <- plotPCA(vd, intgroup=c("treat"))
  #png(paste0("./1_DiffPeak.outs/plotpca_treatandind_filt_", reso, "_", oneX, "_cn.png"), width=2000, height=2000, pointsize=14, res=175)
  #print(p)
  #dev.off()
  data <- plotPCA(vd, intgroup=c("treat"), returnData=TRUE)
  data <- data %>% mutate(sampleID=str_split(data$name, "_", simplify=T)[,4])
head(data)
  #p <- qplot(PC1,PC2,color=treatment,data=data,size=I(10))
  p <- ggplot(data, aes(x=PC1, y=PC2, colour = factor(sampleID))) +
    geom_point() +
#    scale_colour_manual(values = treatment)+
    ggtitle(oneX)
  ggsave(paste0("./1_DiffPeak.outs/plotpca_treatandind_filt_", reso, "_", oneX, "_cn_2.png"),p)
  dds <- DESeq(dds)
  res <- contrast.list%>%map(~results(dds, contrast=.x))
  res2 <- tibble(contrast=names(res),MCls=oneX, data=map(res,tidy))%>%unnest(data)
  ##
  time1 <- Sys.time()
  elapsed <- difftime(time1, time0, units="mins")
  cat(oneX, elapsed, "Done\n")
  res2
})
#res.DE.ATAC
#table(res.DE.ATAC$MCls)
#table(res.DE.ATAC$contrast) 

opfn <- paste0("./1_DiffPeak.outs/2.0_DESeq.results_clusters_treatandind_filtered_", reso, "_cn.rds")
write_rds(res.DE.ATAC, opfn)

res.DE.ATAC <- readRDS("./1_DiffPeak.outs/2.0_DESeq.results_clusters_treatandind_filtered_0.1_cn.rds")

# nrow(res.DE.ATAC)
# nrow(res.DE.ATAC %>% drop_na(p.value))
# nrow(res.DE.ATAC %>% drop_na(p.adjusted))


## ######################################################
## # significant genes
## ######################################################
# res <- res.DE.ATAC
## padj_cutoff <- 0.05

## sig_res <- dplyr::filter(res, p.adjusted < padj_cutoff) %>% dplyr::arrange(p.adjusted)

## #write.csv(sig_res, paste0("./1_DiffPeak.outs/sig_genes.csv"), quote = FALSE, row.names = FALSE)

## #sig_res <- read.csv("./1_DiffPeak.outs/sig_genes.csv")

## head(sig_res)

## sig_res_split <- split(sig_res, f = list(sig_res$contrast,sig_res$MCls))

## head(sig_res_split)
## #sig_res_split[["vitD.TEM"]]

## length(sig_res_split)

## names(sig_res_split)

## sig_res_split2 = list()

## for (i in names(sig_res_split)){
##     print(i)
##     sig_res_split2[[i]] <- sig_res_split[[i]]  %>% dplyr::arrange(p.adjusted) %>% dplyr::pull(gene) %>% head(n=20)
## }

## sig_res_split2



#######################################################
# table of up and down regulated genes #1
#######################################################
#------------------------------------------------------
# Nicole script
#------------------------------------------------------
res <- res.DE.ATAC

FDR <- 0.1
pos_foldchange <- 0.25
neg_foldchange <- -0.25


#------all res------
## res.DE.ATAC.0.05 <- readRDS("./1_DiffPeak.outs/2.0_DESeq.results_clusters_treatandind_filtered_0.05_cn.rds")
## res.DE.ATAC.0.1 <- readRDS("./1_DiffPeak.outs/2.0_DESeq.results_clusters_treatandind_filtered_0.1_cn.rds")
## res.DE.ATAC.0.15 <- readRDS("./1_DiffPeak.outs/2.0_DESeq.results_clusters_treatandind_filtered_0.15_cn.rds")
## res.DE.ATAC.0.2 <- readRDS("./1_DiffPeak.outs/2.0_DESeq.results_clusters_treatandind_filtered_0.2_cn.rds")


## ## res.DE.ATAC.0.05.2 <- res.DE.ATAC.0.05 %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=FDR)  %>% dplyr::filter(abs(estimate) >= pos_foldchange) %>% dplyr::select(contrast, MCls, gene) %>% unique()
## ## res.DE.ATAC.0.1.2 <- res.DE.ATAC.0.1 %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=FDR)  %>% dplyr::filter(abs(estimate) >= pos_foldchange) %>% dplyr::select(contrast, MCls,  gene) %>% unique()
## ## res.DE.ATAC.0.15.2 <- res.DE.ATAC.0.15 %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=FDR)  %>% dplyr::filter(abs(estimate) >= pos_foldchange) %>% dplyr::select(contrast, MCls,  gene) %>% unique()

## res.DE.ATAC.0.05.2 <- res.DE.ATAC.0.05 %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=FDR)  %>% dplyr::filter(abs(estimate) >= pos_foldchange) %>% dplyr::select(contrast, gene) %>% unique()
## res.DE.ATAC.0.1.2 <- res.DE.ATAC.0.1 %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=FDR)  %>% dplyr::filter(abs(estimate) >= pos_foldchange) %>% dplyr::select(contrast, gene) %>% unique()
## res.DE.ATAC.0.15.2 <- res.DE.ATAC.0.15 %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=FDR)  %>% dplyr::filter(abs(estimate) >= pos_foldchange) %>% dplyr::select(contrast, gene) %>% unique()
## res.DE.ATAC.0.2.2 <- res.DE.ATAC.0.2 %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=FDR)  %>% dplyr::filter(abs(estimate) >= pos_foldchange) %>% dplyr::select(contrast, gene) %>% unique()



## a <-as.data.frame(table(res.DE.ATAC.0.05.2$contrast))
## colnames(a) <- c('treats','0.05')
## #a
## b <-as.data.frame(table(res.DE.ATAC.0.1.2$contrast))
## colnames(b) <- c('treats','0.1')
## #b
## c <- as.data.frame(table(res.DE.ATAC.0.15.2$contrast))
## colnames(c) <- c('treats','0.15')
## #c
## d <- as.data.frame(table(res.DE.ATAC.0.2.2$contrast))
## colnames(d) <- c('treats','0.2')
## #d
## #e <- merge(x = a, y = c(b,c), by = "Var1")
## e <- left_join(a,b, by = "treats") %>% left_join(., c, by='treats') %>% left_join(., d, by='treats')
## e
#------all res------






res[!is.na(res$p.adjusted) & res$p.adjusted<=FDR & res$estimate>=pos_foldchange, "Significance"] <- "Upregulated"
res[!is.na(res$p.adjusted) & res$p.adjusted<=FDR & res$estimate<=neg_foldchange, "Significance"] <- "Downregulated"
res[!is.na(res$p.adjusted) & res$p.adjusted>FDR,"Significance"] <- "Not Significant"
res[is.na(res$p.adjusted),"Significance"] <- "NA"
#res2 <- res %>% dplyr::filter(Significance %in% c("Upregulated", "Downregulated"))



#---------------------------
res_sig <- res %>% dplyr::filter(Significance %in% c("Upregulated", "Downregulated"))
nrow(res_sig)

sig <- table(res_sig$MCls, res_sig$contrast)
sig

#write.csv(up, paste0("./1_DiffRNA_2.outs/up_treat&ind_filt_", reso, "_", FDR, "_", pos_foldchange, "_cn.csv"), quote = FALSE, row.names = TRUE)

#up <- read.csv(file = './1_DiffRNA_2.outs/up_treat&ind_filt.csv')
#up


## MCls <- unique(res_up$MCls)
## #MCls
## res_up.unique <- map_dfr(MCls, function(oneX){
##       time0 <- Sys.time()
##         ###
##         res_up.2 <- res_up%>%dplyr::filter(MCls==oneX)
##         ##
##         res_up.3 <- res_up%>%dplyr::filter(MCls!=oneX)
##             #mutate(count=ifelse(nrow(cvt.2)>=8,1,0)) 
##         ##
##         other.genes <- unique(res_up.3$gene)
##         res_up.2.2 <- res_up.2%>%dplyr::filter(!gene %in% other.genes)
##         time1 <- Sys.time()
##         elapsed <- difftime(time1, time0, units="mins")
##         cat(oneX, elapsed, "Done\n")
##         res_up.2.2
##       })
## #res_up.unique
## #nrow(res_up.unique)

## up.unique <- table(res_up.unique$MCls, res_up.unique$contrast)
## #up

## write.csv(up.unique, paste0("./1_DiffRNA_2.outs/up_treat&ind_filt_unique_", reso, "_", FDR, "_", pos_foldchange, ".csv"), quote = FALSE, row.names = TRUE)

library(reshape2)
sig.matrix <- as.matrix(sig)
sig.melt <- melt(sig.matrix)


sig.p <- ggplot(sig.melt, aes(x = Var2, y = Var1)) +
      geom_tile(aes(fill=value), color="black")+
      scale_fill_gradient(low = "white", high = "white")+
      geom_text(aes(label = value), color = "black", size = 5) +
      ggtitle(paste0("Significant", "_", reso, "_", FDR, "-", pos_foldchange)) +
      theme(axis.text.x = element_text(angle = 90),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position="none")

## up.unique.p <- ggplot(up.unique.melt, aes(x = Var2, y = Var1)) +
##     geom_tile(aes(fill=value))+
##     geom_text(aes(label = value), color = "white", size = 4) +
##     ggtitle("Upregulated.unique") +
##     theme(axis.text.x = element_text(angle = 90))


library("ggpubr")

figure <- ggarrange(sig.p +
                    font("x.text", size = 14))
                   # ncol = 2, nrow = 2)

png(paste0("./1_DiffPeak.outs/plot_sig_treatandind_filt_all_", reso, "_", FDR, "_", pos_foldchange, "_cn.png"), width=500, height=500, pointsize=16, res=125)
print(figure)
dev.off()






#---------------------------
res_up <- res %>% dplyr::filter(Significance %in% "Upregulated")
#nrow(res_up)

up <- table(res_up$MCls, res_up$contrast)
#up

write.csv(up, paste0("./1_DiffPeak.outs/up_treat&ind_filt_", reso, "_", FDR, "_", pos_foldchange, "_cn.csv"), quote = FALSE, row.names = TRUE)

#up <- read.csv(file = './1_DiffPeak.outs/up_treat&ind_filt.csv')
#up


## MCls <- unique(res_up$MCls)
## #MCls
## res_up.unique <- map_dfr(MCls, function(oneX){
##       time0 <- Sys.time()
##         ###
##         res_up.2 <- res_up%>%dplyr::filter(MCls==oneX)
##         ##
##         res_up.3 <- res_up%>%dplyr::filter(MCls!=oneX)
##             #mutate(count=ifelse(nrow(cvt.2)>=8,1,0)) 
##         ##
##         other.genes <- unique(res_up.3$gene)
##         res_up.2.2 <- res_up.2%>%dplyr::filter(!gene %in% other.genes)
##         time1 <- Sys.time()
##         elapsed <- difftime(time1, time0, units="mins")
##         cat(oneX, elapsed, "Done\n")
##         res_up.2.2
##       })
## #res_up.unique
## #nrow(res_up.unique)

## up.unique <- table(res_up.unique$MCls, res_up.unique$contrast)
## #up

## write.csv(up.unique, paste0("./1_DiffPeak.outs/up_treat&ind_filt_unique_", reso, "_", FDR, "_", pos_foldchange, ".csv"), quote = FALSE, row.names = TRUE)



#---------------------------
res_down <- res %>% dplyr::filter(Significance %in% "Downregulated")

down <- table(res_down$MCls, res_down$contrast)
#down

write.csv(down, paste0("./1_DiffPeak.outs/down_treat&ind_filt_", reso, "_", FDR, "_", pos_foldchange, "_cn.csv"), quote = FALSE, row.names = TRUE)

#down <- read.csv(file = './1_DiffPeak.outs/down_treat&ind_filt.csv')
#down


## MCls <- unique(res_down$MCls)
## #MCls
## res_down.unique <- map_dfr(MCls, function(oneX){
##       time0 <- Sys.time()
##         ###
##         res_down.2 <- res_down%>%dplyr::filter(MCls==oneX)
##         ##
##         res_down.3 <- res_down%>%dplyr::filter(MCls!=oneX)
##             #mutate(count=ifelse(nrow(cvt.2)>=8,1,0)) 
##         ##
##         other.genes <- unique(res_down.3$gene)
##         res_down.2.2 <- res_down.2%>%dplyr::filter(!gene %in% other.genes)
##         time1 <- Sys.time()
##         elapsed <- difftime(time1, time0, units="mins")
##         cat(oneX, elapsed, "Done\n")
##         res_down.2.2
##       })
## #res_down.unique
## #nrow(res_down.unique)


## down.unique <- table(res_down.unique$MCls, res_down.unique$contrast)
## #down
## write.csv(down.unique, paste0("./1_DiffPeak.outs/down_treat&ind_filt_unique_", reso, "_", FDR, "_", pos_foldchange, ".csv"), quote = FALSE, row.names = TRUE)


# -----------number of cells------------------------------------------------
#cvt.2

dd3 <- left_join(cvt.2, dd2)
#head(dd3)

dd3 <- dd3 %>% dplyr::filter(dd3$MCls %in% res$MCls, !dd3$treat %in% "etOH")
#head(dd3)
#nrow(dd3)
#sum(dd3$ncell)

dd4 <- dd3 %>% group_by(treat, MCls) %>% summarise(ncell_count = sum(ncell))
#dd4

ncell_table <- data.frame(pivot_wider(dd4, names_from=treat, values_from=ncell_count))

rownames(ncell_table) <- ncell_table$MCls
ncell_table$MCls <- NULL
#ncell_table

write.csv(ncell_table, paste0("./1_DiffPeak.outs/ncell_treat$ind_filt_", reso, "_", FDR, "_", pos_foldchange, "_cn.csv"), quote = FALSE, row.names = TRUE)

#ncell_table <- read.csv("./1_DiffPeak.outs/ncell_treat$ind_filt.csv")


#########################################################
# matrix_plots
#########################################################
library(reshape2)

up.matrix <-as.matrix(up) # create a numeric matrix object
up.melt <- melt(up.matrix)

#up.unique.matrix <-as.matrix(up.unique) # create a numeric matrix object
#up.unique.melt <- melt(up.unique.matrix)

down.matrix <- as.matrix(down)
down.melt <- melt(down.matrix)

#down.unique.matrix <- as.matrix(down.unique)
#down.unique.melt <- melt(down.unique.matrix)


#ncell
#rownames(ncell)<- ncell$X
#ncell$X <- NULL
ncell.matrix <- as.matrix(ncell_table)
ncell.melt <- melt(ncell.matrix)

#comb
#rownames(comb)<- comb$X
#comb$X <- NULL
#comb.2 <- subset(comb, select = -c(etOH))
#comb.2 <- comb.2 %>% dplyr::filter(rownames(comb.3) %in% rownames(up))
comb.3.matrix <- as.matrix(comb.2)
comb.3.matrix <- subset.matrix(comb.3.matrix, select = -c(etOH))
comb.3.melt <- melt(comb.3.matrix)

up.p <- ggplot(up.melt, aes(x = Var2, y = Var1)) +
    geom_tile(aes(fill=value), color="black")+
    scale_fill_gradient(low = "white", high = "white")+
    geom_text(aes(label = value), color = "black", size = 5) +
    ggtitle(paste0("Upregulated", "_", reso, "_", FDR, "-", pos_foldchange)) +
    theme(axis.text.x = element_text(angle = 90),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position="none")



down.p <- ggplot(down.melt, aes(x = Var2, y = Var1)) +
      geom_tile(aes(fill=value), color="black")+
      scale_fill_gradient(low = "white", high = "white")+
      geom_text(aes(label = value), color = "black", size = 5) +
      ggtitle(paste0("Downregulated", "_", reso, "_", FDR, "-", pos_foldchange)) +
      theme(axis.text.x = element_text(angle = 90),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position="none")

## up.unique.p <- ggplot(up.unique.melt, aes(x = Var2, y = Var1)) +
##     geom_tile(aes(fill=value))+
##     geom_text(aes(label = value), color = "white", size = 4) +
##     ggtitle("Upregulated.unique") +
##     theme(axis.text.x = element_text(angle = 90))

## down.unique.p <- ggplot(down.unique.melt, aes(x = Var2, y = Var1)) +
##       geom_tile(aes(fill=value))+
##       geom_text(aes(label = value), color = "white", size = 4) +
##       ggtitle("Downregulated.unique") +
##       theme(axis.text.x = element_text(angle = 90))


ncell.p <- ggplot(ncell.melt, aes(x = Var2, y = Var1)) +
      geom_tile(aes(fill=value), color="black")+
      scale_fill_gradient(low = "white", high = "white")+
      geom_text(aes(label = value), color = "black", size = 5) +
      ggtitle(paste0("ncells", "_", reso, "_", FDR, "-", pos_foldchange)) +
      theme(axis.text.x = element_text(angle = 90),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position="none")


comb.3.p <- ggplot(comb.3.melt, aes(x = Var2, y = Var1)) +
      geom_tile(aes(fill=value), color="black")+
      scale_fill_gradient(low = "white", high = "white")+
      geom_text(aes(label = value), color = "black", size = 5) +
      ggtitle(paste0("Combinations", "_", reso, "_", FDR, "-", pos_foldchange)) +
      theme(axis.text.x = element_text(angle = 90),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position="none")


library("ggpubr")
figure <- ggarrange(comb.3.p, ncell.p, up.p, down.p +
                    font("x.text", size = 14),
                    ncol = 2, nrow = 2)

png(paste0("./1_DiffPeak.outs/plot_all_treatandind_filt_all_", reso, "_", FDR, "_", pos_foldchange, "_cn.png"), width=2000, height=2000, pointsize=16, res=175)
print(figure)
dev.off()


##############################################################
###
### 1. MA plots
##############################################################
figfn <- paste0("./1_DiffPeak.outs/Figure1.1_MA_clusters_treat&ind_filt_", reso, "_", FDR, "_", pos_foldchange, "_cn.png")
png(figfn, width=3500, height=4500, pointsize=18, res=225)
par(mar=c(4,4,2,2),mgp=c(2,1,0))
x <- matrix(1:63, 9, 7, byrow=T)
layout(x)
#fn <- "./1_DiffPeak.outs/2.0_DESeq.results_clusters_treat&ind_filtered.rds"
#res <- read_rds(fn)%>%mutate(color=ifelse((p.adjusted<0.1)&(!is.na(p.adjusted)), T, F))
res.1 <- res%>%mutate(color=ifelse((p.adjusted<=FDR)&(!is.na(p.adjusted)), T, F))
# MCls <- c("CD4Naive", "TCM", "NKcell", "Bcell", "CD8Naive", "TEM", "Monocyte", "CTL", "dnT")
MCls <- unique(cvt.2$MCls)
for (oneMCl in MCls){
  ##1
  res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="caffeine")%>%
    dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
  print(plotMA(res.2[,1:3], colLine="NA", main="caffeine", cex.main=1, cex.axis=0.8, cex.lab=1))
  ###2    
  res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="nicotine")%>%
    dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
  print(plotMA(res.2[,1:3], colLine="NA", main="nicotine", cex.main=1, cex.axis=0.8, cex.lab=1))
  ###3
  res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="vitA")%>%
    dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
  print(plotMA(res.2[,1:3], colLine="NA", main="vitA", cex.main=1, cex.axis=0.8, cex.lab=1))
  ###4
  res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="vitD")%>%
    dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
  print(plotMA(res.2[,1:3], colLine="NA", main="vitD", cex.main=1, cex.axis=0.8, cex.lab=1))
  ###5
  res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="vitE")%>%
    dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
  print(plotMA(res.2[,1:3], colLine="NA", main="vitE", cex.main=1, cex.axis=0.8, cex.lab=1))
  ###6    
  res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="water")%>%
    dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
  print(plotMA(res.2[,1:3], colLine="NA", main="water", cex.main=1, cex.axis=0.8, cex.lab=1))
  ###7    
  res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="zinc")%>%
    dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
  print(plotMA(res.2[,1:3], colLine="NA", main="zinc", cex.main=1, cex.axis=0.8, cex.lab=1))
  print(mtext(oneMCl, side=4, line=0.5, cex=1, col="blue")) 
}
dev.off()




### 2, qq plots
figfn <- paste0("./1_DiffPeak.outs/Figure1.2_qq_clusters_treat&ind_filt_", reso, "_", FDR, "_", pos_foldchange, "_cn.png")
png(figfn, width=3500, height=4500, pointsize=18, res=225)
par(mar=c(4,4,2,2),mgp=c(2,1,0))
x <- matrix(1:63, 9, 7, byrow=T)
layout(x)
#fn <- "./1_DiffPeak.outs/2.0_DESeq.results_clusters_treat&ind_filtered.rds"
#res <- read_rds(fn)%>%drop_na(p.value)
res.1 <- res%>%drop_na(p.value)
# MCls <- c("CD4Naive", "TCM", "NKcell", "Bcell", "CD8Naive", "TEM", "Monocyte", "CTL", "dnT")
MCls <- unique(cvt.2$MCls)
for (oneMCl in MCls){
  ##1
  res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="caffeine")
  print( qq(res.2$p.value, main="caffeine", cex.main=1, cex.axis=0.8, cex.lab=1))
  ##2
  res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="nicotine")
  print( qq(res.2$p.value, main="nicotine", cex.main=1, cex.axis=0.8, cex.lab=1))
  ##3 
  res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="vitA")
  print( qq(res.2$p.value, main="vitA", cex.main=1, cex.axis=0.8, cex.lab=1))
  ##4
  res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="vitD")
  print( qq(res.2$p.value, main="vitD", cex.main=1, cex.axis=0.8, cex.lab=1))
  ##5
  res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="vitE")
  print( qq(res.2$p.value, main="vitE", cex.main=1, cex.axis=0.8, cex.lab=1))
  ##6
  res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="water")
  print( qq(res.2$p.value, main="water", cex.main=1, cex.axis=0.8, cex.lab=1))
  ##7
  res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="zinc")
  print( qq(res.2$p.value, main="zinc", cex.main=1, cex.axis=0.8, cex.lab=1))   
  print( mtext(oneMCl, side=4, line=0.5, cex=1, col="blue"))  
}
dev.off()


###
### 3, canno plots


#------------------plot_loop--------------------------------
#nrow(res) #if qq plot has been called before this step, this already is with dropped NA values
res.1 <- res %>% drop_na(p.value)
max.y <- floor(max(-log10(res.1$p.value)))/2
#max.y
png(paste0("./1_DiffPeak.outs/Figure1.3_qq_clusters_together_treatandind_filt_", reso,"_", FDR, "_", pos_foldchange, "_cn.png"), width=2000, height=5000, pointsize=18, res=200)
par(mar=c(4,4,2,2),mgp=c(2,1,0))
x <- matrix(1:10, 5, 2, byrow=T)
layout(x)
# MCls <- c("CD4Naive", "TCM", "NKcell", "Bcell", "CD8Naive", "TEM", "Monocyte", "CTL", "dnT")
MCls <- unique(cvt.2$MCls)
for (oneMCl in MCls){
  ##1
  res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="caffeine")%>% drop_na(p.value)
  plot(-log10(1:length(res.2$p.value)/length(res.2$p.value)),
       -log10(sort(res.2$p.value)),
       main=oneMCl,
       pch = 19,
       ylim = c(0, max.y),
       ylab="observed -log10(p)",
       xlab="expected -log10(p)",
       col="red")
  ##2
  res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="nicotine")%>% drop_na(p.value)
  points(-log10(1:length(res.2$p.value)/length(res.2$p.value)),-log10(sort(res.2$p.value)),  col="tan", pch = 19)   
  ##3 
  res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="zinc")%>% drop_na(p.value)
  points(-log10(1:length(res.2$p.value)/length(res.2$p.value)),-log10(sort(res.2$p.value)),  col="maroon3", pch = 19)
  ##4
  res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="vitA")%>% drop_na(p.value)
  points(-log10(1:length(res.2$p.value)/length(res.2$p.value)),-log10(sort(res.2$p.value)),  col="tan4", pch = 19)
  ##5
  res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="vitD")%>% drop_na(p.value)
  points(-log10(1:length(res.2$p.value)/length(res.2$p.value)),-log10(sort(res.2$p.value)),  col="seagreen4", pch = 19)
  ##6
  res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="water")%>% drop_na(p.value)
  points(-log10(1:length(res.2$p.value)/length(res.2$p.value)),-log10(sort(res.2$p.value)),  col="grey", pch = 19)
  ##7
  res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="vitE")%>% drop_na(p.value)
  points(-log10(1:length(res.2$p.value)/length(res.2$p.value)),-log10(sort(res.2$p.value)),  col="salmon3", pch = 19)
  ##control
  abline(0, 1, col = "black")
  ##legend
  legend("topleft",
         #inset=c(-0.2,0),
         cex = 1.2,
         pch = 19,
         c("caff","nic", "zinc", "vitA", "vitD", "vitE", "water"),
         fill=c("red","tan", "maroon3", "tan4", "seagreen4", "salmon3", "grey"))
}
dev.off()







#######################################################
# table of up and down regulated genes #2
#######################################################
#------------------------------------------------------
# Nicole script
#------------------------------------------------------
res <- res.DE.ATAC

FDR <- 0.05
pos_foldchange <- 0.00
neg_foldchange <- 0.00


#------all res------
## ## res.DE.ATAC.0.05 <- readRDS("./1_DiffPeak.outs/2.0_DESeq.results_clusters_treatandind_filtered_0.05.rds")
## ## res.DE.ATAC.0.1 <- readRDS("./1_DiffPeak.outs/2.0_DESeq.results_clusters_treatandind_filtered_0.1.rds")
## ## res.DE.ATAC.0.15 <- readRDS("./1_DiffPeak.outs/2.0_DESeq.results_clusters_treatandind_filtered_0.15.rds")

## ## res.DE.ATAC.0.05.2 <- res.DE.ATAC.0.05 %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=FDR)  %>% dplyr::filter(abs(estimate) >= pos_foldchange) %>% dplyr::select(contrast, MCls, gene) %>% unique()
## ## res.DE.ATAC.0.1.2 <- res.DE.ATAC.0.1 %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=FDR)  %>% dplyr::filter(abs(estimate) >= pos_foldchange) %>% dplyr::select(contrast, MCls,  gene) %>% unique()
## ## res.DE.ATAC.0.15.2 <- res.DE.ATAC.0.15 %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=FDR)  %>% dplyr::filter(abs(estimate) >= pos_foldchange) %>% dplyr::select(contrast, MCls,  gene) %>% unique()

## res.DE.ATAC.0.05.2 <- res.DE.ATAC.0.05 %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=FDR)  %>% dplyr::filter(abs(estimate) >= pos_foldchange) %>% dplyr::select(contrast, gene) %>% unique()
## res.DE.ATAC.0.1.2 <- res.DE.ATAC.0.1 %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=FDR)  %>% dplyr::filter(abs(estimate) >= pos_foldchange) %>% dplyr::select(contrast, gene) %>% unique()
## res.DE.ATAC.0.15.2 <- res.DE.ATAC.0.15 %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=FDR)  %>% dplyr::filter(abs(estimate) >= pos_foldchange) %>% dplyr::select(contrast, gene) %>% unique()
## res.DE.ATAC.0.2.2 <- res.DE.ATAC.0.2 %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=FDR)  %>% dplyr::filter(abs(estimate) >= pos_foldchange) %>% dplyr::select(contrast, gene) %>% unique()

## a <-as.data.frame(table(res.DE.ATAC.0.05.2$contrast))
## colnames(a) <- c('treats','0.05')
## #a
## b <-as.data.frame(table(res.DE.ATAC.0.1.2$contrast))
## colnames(b) <- c('treats','0.1')
## #b
## c <- as.data.frame(table(res.DE.ATAC.0.15.2$contrast))
## colnames(c) <- c('treats','0.15')
## #d
## d <- as.data.frame(table(res.DE.ATAC.0.2.2$contrast))
## colnames(d) <- c('treats','0.2')
## #d
## #e <- merge(x = a, y = c(b,c), by = "Var1")
## e <- left_join(a,b, by = "treats") %>% left_join(., c, by='treats') %>% left_join(., d, by='treats')
## e
#------all res------






res[!is.na(res$p.adjusted) & res$p.adjusted<=FDR & res$estimate>=pos_foldchange, "Significance"] <- "Upregulated"
res[!is.na(res$p.adjusted) & res$p.adjusted<=FDR & res$estimate<=neg_foldchange, "Significance"] <- "Downregulated"
res[!is.na(res$p.adjusted) & res$p.adjusted>FDR,"Significance"] <- "Not Significant"
res[is.na(res$p.adjusted),"Significance"] <- "NA"
#res2 <- res %>% dplyr::filter(Significance %in% c("Upregulated", "Downregulated"))

#---------------------------
res_up <- res %>% dplyr::filter(Significance %in% "Upregulated")
#nrow(res_up)

up <- table(res_up$MCls, res_up$contrast)
#up

write.csv(up, paste0("./1_DiffPeak.outs/up_treat&ind_filt_", reso, "_", FDR, "_", pos_foldchange, "_cn.csv"), quote = FALSE, row.names = TRUE)

#up <- read.csv(file = './1_DiffPeak.outs/up_treat&ind_filt.csv')
#up


## MCls <- unique(res_up$MCls)
## #MCls
## res_up.unique <- map_dfr(MCls, function(oneX){
##       time0 <- Sys.time()
##         ###
##         res_up.2 <- res_up%>%dplyr::filter(MCls==oneX)
##         ##
##         res_up.3 <- res_up%>%dplyr::filter(MCls!=oneX)
##             #mutate(count=ifelse(nrow(cvt.2)>=8,1,0)) 
##         ##
##         other.genes <- unique(res_up.3$gene)
##         res_up.2.2 <- res_up.2%>%dplyr::filter(!gene %in% other.genes)
##         time1 <- Sys.time()
##         elapsed <- difftime(time1, time0, units="mins")
##         cat(oneX, elapsed, "Done\n")
##         res_up.2.2
##       })
## #res_up.unique
## #nrow(res_up.unique)

## up.unique <- table(res_up.unique$MCls, res_up.unique$contrast)
## #up

## write.csv(up.unique, paste0("./1_DiffPeak.outs/up_treat&ind_filt_unique_", reso, "_", FDR, "_", pos_foldchange, ".csv"), quote = FALSE, row.names = TRUE)



#---------------------------
res_down <- res %>% dplyr::filter(Significance %in% "Downregulated")

down <- table(res_down$MCls, res_down$contrast)
#down

write.csv(down, paste0("./1_DiffPeak.outs/down_treat&ind_filt_", reso, "_", FDR, "_", pos_foldchange, "_cn.csv"), quote = FALSE, row.names = TRUE)

#down <- read.csv(file = './1_DiffPeak.outs/down_treat&ind_filt.csv')
#down


## MCls <- unique(res_down$MCls)
## #MCls
## res_down.unique <- map_dfr(MCls, function(oneX){
##       time0 <- Sys.time()
##         ###
##         res_down.2 <- res_down%>%dplyr::filter(MCls==oneX)
##         ##
##         res_down.3 <- res_down%>%dplyr::filter(MCls!=oneX)
##             #mutate(count=ifelse(nrow(cvt.2)>=8,1,0)) 
##         ##
##         other.genes <- unique(res_down.3$gene)
##         res_down.2.2 <- res_down.2%>%dplyr::filter(!gene %in% other.genes)
##         time1 <- Sys.time()
##         elapsed <- difftime(time1, time0, units="mins")
##         cat(oneX, elapsed, "Done\n")
##         res_down.2.2
##       })
## #res_down.unique
## #nrow(res_down.unique)


## down.unique <- table(res_down.unique$MCls, res_down.unique$contrast)
## #down
## write.csv(down.unique, paste0("./1_DiffPeak.outs/down_treat&ind_filt_unique_", reso, "_", FDR, "_", pos_foldchange, ".csv"), quote = FALSE, row.names = TRUE)


# -----------number of cells------------------------------------------------
#cvt.2

dd3 <- left_join(cvt.2, dd2)
#head(dd3)

dd3 <- dd3 %>% dplyr::filter(dd3$MCls %in% res$MCls, !dd3$treat %in% "etOH")
#head(dd3)
#nrow(dd3)
#sum(dd3$ncell)

dd4 <- dd3 %>% group_by(treat, MCls) %>% summarise(ncell_count = sum(ncell))
#dd4

ncell_table <- data.frame(pivot_wider(dd4, names_from=treat, values_from=ncell_count))

rownames(ncell_table) <- ncell_table$MCls
ncell_table$MCls <- NULL
#ncell_table

write.csv(ncell_table, paste0("./1_DiffPeak.outs/ncell_treat$ind_filt_", reso, "_", FDR, "_", pos_foldchange, "_cn.csv"), quote = FALSE, row.names = TRUE)

#ncell_table <- read.csv("./1_DiffPeak.outs/ncell_treat$ind_filt.csv")


#########################################################
# matrix_plots
#########################################################
library(reshape2)

up.matrix <-as.matrix(up) # create a numeric matrix object
up.melt <- melt(up.matrix)

#up.unique.matrix <-as.matrix(up.unique) # create a numeric matrix object
#up.unique.melt <- melt(up.unique.matrix)

down.matrix <- as.matrix(down)
down.melt <- melt(down.matrix)

#down.unique.matrix <- as.matrix(down.unique)
#down.unique.melt <- melt(down.unique.matrix)


#ncell
#rownames(ncell)<- ncell$X
#ncell$X <- NULL
ncell.matrix <- as.matrix(ncell_table)
ncell.melt <- melt(ncell.matrix)

#comb
#rownames(comb)<- comb$X
#comb$X <- NULL
#comb.2 <- subset(comb, select = -c(etOH))
#comb.2 <- comb.2 %>% dplyr::filter(rownames(comb.3) %in% rownames(up))
comb.3.matrix <- as.matrix(comb.2)
comb.3.matrix <- subset.matrix(comb.3.matrix, select = -c(etOH))
comb.3.melt <- melt(comb.3.matrix)

up.p <- ggplot(up.melt, aes(x = Var2, y = Var1)) +
    geom_tile(aes(fill=value), color="black")+
    scale_fill_gradient(low = "white", high = "white")+
    geom_text(aes(label = value), color = "black", size = 5) +
    ggtitle(paste0("Upregulated", "_", reso, "_", FDR, "-", pos_foldchange)) +
    theme(axis.text.x = element_text(angle = 90),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position="none")



down.p <- ggplot(down.melt, aes(x = Var2, y = Var1)) +
      geom_tile(aes(fill=value), color="black")+
      scale_fill_gradient(low = "white", high = "white")+
      geom_text(aes(label = value), color = "black", size = 5) +
      ggtitle(paste0("Downregulated", "_", reso, "_", FDR, "-", pos_foldchange)) +
      theme(axis.text.x = element_text(angle = 90),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position="none")

## up.unique.p <- ggplot(up.unique.melt, aes(x = Var2, y = Var1)) +
##     geom_tile(aes(fill=value))+
##     geom_text(aes(label = value), color = "white", size = 4) +
##     ggtitle("Upregulated.unique") +
##     theme(axis.text.x = element_text(angle = 90))

## down.unique.p <- ggplot(down.unique.melt, aes(x = Var2, y = Var1)) +
##       geom_tile(aes(fill=value))+
##       geom_text(aes(label = value), color = "white", size = 4) +
##       ggtitle("Downregulated.unique") +
##       theme(axis.text.x = element_text(angle = 90))


ncell.p <- ggplot(ncell.melt, aes(x = Var2, y = Var1)) +
      geom_tile(aes(fill=value), color="black")+
      scale_fill_gradient(low = "white", high = "white")+
      geom_text(aes(label = value), color = "black", size = 5) +
      ggtitle(paste0("ncells", "_", reso, "_", FDR, "-", pos_foldchange)) +
      theme(axis.text.x = element_text(angle = 90),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position="none")


comb.3.p <- ggplot(comb.3.melt, aes(x = Var2, y = Var1)) +
      geom_tile(aes(fill=value), color="black")+
      scale_fill_gradient(low = "white", high = "white")+
      geom_text(aes(label = value), color = "black", size = 5) +
      ggtitle(paste0("Combinations", "_", reso, "_", FDR, "-", pos_foldchange)) +
      theme(axis.text.x = element_text(angle = 90),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position="none")


library("ggpubr")
figure <- ggarrange(comb.3.p, ncell.p, up.p, down.p +
                    font("x.text", size = 14),
                    ncol = 2, nrow = 2)

png(paste0("./1_DiffPeak.outs/plot_all_treatandind_filt_all_", reso, "_", FDR, "_", pos_foldchange, "_cn.png"), width=2000, height=2000, pointsize=16, res=175)
print(figure)
dev.off()


##############################################################
###
### 1. MA plots
##############################################################
figfn <- paste0("./1_DiffPeak.outs/Figure1.1_MA_clusters_treat&ind_filt_", reso, "_", FDR, "_", pos_foldchange, "_cn.png")
png(figfn, width=3500, height=4500, pointsize=18, res=225)
par(mar=c(4,4,2,2),mgp=c(2,1,0))
x <- matrix(1:63, 9, 7, byrow=T)
layout(x)
#fn <- "./1_DiffPeak.outs/2.0_DESeq.results_clusters_treat&ind_filtered.rds"
#res <- read_rds(fn)%>%mutate(color=ifelse((p.adjusted<0.1)&(!is.na(p.adjusted)), T, F))
res.1 <- res%>%mutate(color=ifelse((p.adjusted<=FDR)&(!is.na(p.adjusted)), T, F))
# MCls <- c("CD4Naive", "TCM", "NKcell", "Bcell", "CD8Naive", "TEM", "Monocyte", "CTL", "dnT")
MCls <- unique(cvt.2$MCls)
for (oneMCl in MCls){
  ##1
  res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="caffeine")%>%
    dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
  print(plotMA(res.2[,1:3], colLine="NA", main="caffeine", cex.main=1, cex.axis=0.8, cex.lab=1))
  ###2    
  res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="nicotine")%>%
    dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
  print(plotMA(res.2[,1:3], colLine="NA", main="nicotine", cex.main=1, cex.axis=0.8, cex.lab=1))
  ###3
  res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="vitA")%>%
    dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
  print(plotMA(res.2[,1:3], colLine="NA", main="vitA", cex.main=1, cex.axis=0.8, cex.lab=1))
  ###4
  res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="vitD")%>%
    dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
  print(plotMA(res.2[,1:3], colLine="NA", main="vitD", cex.main=1, cex.axis=0.8, cex.lab=1))
  ###5
  res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="vitE")%>%
    dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
  print(plotMA(res.2[,1:3], colLine="NA", main="vitE", cex.main=1, cex.axis=0.8, cex.lab=1))
  ###6    
  res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="water")%>%
    dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
  print(plotMA(res.2[,1:3], colLine="NA", main="water", cex.main=1, cex.axis=0.8, cex.lab=1))
  ###7    
  res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="zinc")%>%
    dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
  print(plotMA(res.2[,1:3], colLine="NA", main="zinc", cex.main=1, cex.axis=0.8, cex.lab=1))
  print(mtext(oneMCl, side=4, line=0.5, cex=1, col="blue")) 
}
dev.off()








#####################################################################################################################################################
### if enriched in DARs ###
#####################################################################################################################################################
library(EnsDb.Hsapiens.v86)

DefaultAssay(combined) <- "ATAC"
#granges(combined[["ATAC"]])

gene.ranges = genes(EnsDb.Hsapiens.v86)
#head(gene.ranges)

seqlevelsStyle(gene.ranges) <- 'UCSC'
#head(gene.ranges)

peakAnno <- ClosestFeature(combined,
                           regions=granges(combined),
                           annotation = gene.ranges,
                           sep = c(':', '-')
                           )

# head(peakAnno)
# colnames(peakAnno)

peakAnno <- peakAnno%>%dplyr::select(gene_name,query_region)%>%mutate(query_region=gsub(":", "-", query_region))
# head(peakAnno)



# res.DE object -> DEGS
#---------------------------------------------
#res.DE <- readRDS("./1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1.rds")
#res.DE

res.DE <- res.DE%>%mutate(comb=paste(MCls, contrast, sep="_"))
#head(resDE)

#res.DE.2 <- res.DE %>% dplyr::filter(p.adjusted<0.1,abs(estimate)>0.25)
#unique(resDE$gene)





# res.DE.ATAC object
#---------------------------------------------
#res.DE.ATAC <- readRDS("./1_DiffPeak.outs/2.0_DESeq.results_clusters_treatandind_filtered_0.1.rds")
#res.DE.ATAC

res.DE.ATAC <- res.DE.ATAC%>%dplyr::rename("peak_region"="gene")%>%drop_na(p.value)

res.DE.ATAC <- res.DE.ATAC%>%
      left_join(peakAnno,by=c("peak_region"="query_region"))%>%
      mutate(comb=paste(MCls, contrast, sep="_"))
#res
#head(res$gene_name)

## colnames(res.DE.ATAC)
## nrow(res.DE.ATAC)
## nrow(res.DE.ATAC %>% drop_na(gene_name))

# res.DE.ATAC.1 <- res.DE.ATAC %>% dplyr::filter(p.adjusted<=FDR, abs(estimate)>=pos_foldchange, abs(estimate)<=neg_foldchange)

res.DE.ATAC.1 <- res.DE.ATAC %>% dplyr::filter(p.adjusted<=0.1, abs(estimate)>=0.25)

res.DE.ATAC.2 <- res.DE.ATAC %>% dplyr::filter(p.adjusted<=0.05)




# function
###-------------------------------------------
comb <- sort(unique(res.DE.ATAC.1$comb))
#comb
dfNew.1 <- map_dfr(comb, function(ii){
      res2 <- res.DE%>%dplyr::filter(comb==ii)
        DEP <- res.DE.ATAC.1%>%dplyr::filter(comb==ii)%>%dplyr::pull(gene_name)
        res2 <- res2%>%mutate(is_DEP=ifelse(gene%in%DEP,1,0))
        ###
        dx <- map_dfr(c(0,1),function(i){
                di <- res2%>%dplyr::filter(is_DEP==i)%>%arrange(p.value)
                    ngene <- nrow(di)
                    di <- di%>%mutate(observed=-log10(p.value), expected=-log10(ppoints(ngene)))
                    di
                  })
        dx
      })
# head(dfNew)
# table(dfNew$contrast)

comb <- sort(unique(res.DE.ATAC.2$comb))
dfNew.2 <- map_dfr(comb, function(ii){
      res2 <- res.DE%>%dplyr::filter(comb==ii)
        DEP <- res.DE.ATAC.2%>%dplyr::filter(comb==ii)%>%dplyr::pull(gene_name)
        res2 <- res2%>%mutate(is_DEP=ifelse(gene%in%DEP,1,0))
        ###
        dx <- map_dfr(c(0,1),function(i){
                di <- res2%>%dplyr::filter(is_DEP==i)%>%arrange(p.value)
                    ngene <- nrow(di)
                    di <- di%>%mutate(observed=-log10(p.value), expected=-log10(ppoints(ngene)))
                    di
                  })
        dx
      })

### plot
###-----------------------------------------------
#lab1 <- c(unique(dfNew$contrast))
#labl
lab1 <- c("caffeine"="caffeine", "nicotine"="nicotine", "vitA"="vitA", "vitD"="vitD", "vitE"="vitE", "water"="water", "zinc"="zinc")

p1 <- ggplot(dfNew.1, aes(x=expected,y=observed, colour=factor(is_DEP)))+
      ggrastr::rasterise(geom_point(size=0.3),dpi=300)+
        geom_abline(colour="red")+
        scale_colour_manual(values=c("0"="grey40","1"="green"),
                                                  labels=c("0"="Not DAR", "1"="DAR"),
                                                  guide=guide_legend(override.aes=list(size=1)))+
        facet_grid(MCls~contrast, scales="free", labeller=labeller(contrast=lab1))+
        xlab(bquote("Expected"~ -log[10]~"("~italic(plain(P))~")"))+
        ylab(bquote("Observed"~-log[10]~"("~italic(plain(P))~")"))+
        theme_bw()+
        theme(legend.title=element_blank(),strip.text=element_text(size=12))
figfn <- paste0("./1_DiffPeak.outs/Figure1.3_qq.DEP_", reso,"_0.1_0.25_cn.png")
png(figfn, width=750, height=750, res=120)
print(p1)
dev.off()

p2 <- ggplot(dfNew.2, aes(x=expected,y=observed, colour=factor(is_DEP)))+
      ggrastr::rasterise(geom_point(size=0.3),dpi=300)+
        geom_abline(colour="red")+
        scale_colour_manual(values=c("0"="grey40","1"="green"),
                                                  labels=c("0"="Not DAR", "1"="DAR"),
                                                  guide=guide_legend(override.aes=list(size=1)))+
        facet_grid(MCls~contrast, scales="free", labeller=labeller(contrast=lab1))+
        xlab(bquote("Expected"~ -log[10]~"("~italic(plain(P))~")"))+
        ylab(bquote("Observed"~-log[10]~"("~italic(plain(P))~")"))+
        theme_bw()+
        theme(legend.title=element_blank(),strip.text=element_text(size=12))
figfn <- paste0("./1_DiffPeak.outs/Figure1.3_qq.DEP_", reso, "_0.05_0_cn.png")
png(figfn, width=750, height=750, res=120)
print(p2)
dev.off()


cat("#################################", "\n")
cat("ATAC analysis done", "\n")
cat("#################################", "\n")

v
#print("total number of genes <0.1 FDR")
cat("total number of genes <0.1 FDR", "\n")
nrow(res.DE %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=0.1))

#print("total number of genes <0.05 FDR")
cat("total number of genes <0.05 FDR", "\n")
nrow(res.DE %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=0.05))

#print("total number of ATAC features <0.1 FDR")
cat("total number of ATAC features <0.1 FDR", "\n")
nrow(res.DE.ATAC %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=0.1))

#print("total number of ATAC features <0.05 FDR")
cat("total number of ATAC features <0.05 FDR", "\n")
nrow(res.DE.ATAC %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=0.05))

#print("total number of genes with p.adj")
cat("total number of genes with p.adj", "\n")
nrow(res.DE %>% drop_na(p.adjusted))

#print("total number of ATAC features with p.adj")
cat("total number of ATAC features with p.adj", "\n")
nrow(res.DE.ATAC %>% drop_na(p.adjusted))

combined[["ATAC"]]
combined[["RNA"]]










#############################################################################################################################
#############################################################################################################################
#### Sarah's data <- bulk RNA seq ###  
#############################################################################################################################
#############################################################################################################################
####################################
### 1. Generate pseudo-bulk data ###
####################################
deg.counts <- read.table("/wsu/home/fq/fq93/fq9336/rprdata/mohammedhusain/DEG/after_correcting_script/batch_effect_ethnicity_only_protein_coding/dge_counts.txt")


combined <- read_rds("../1_processing/3_Clustering.outs/1_seurat.cluster.combined.mitofilt.rds")
#colnames(combined@meta.data)

#combined <- subset(x = combined, subset = percent.mt.RNA < 20)
#colnames(combined@meta.data)
#max(combined@meta.data$percent.mt.RNA)

#--------------------------------------------------------------------------
## meta <- atac@meta.data%>%mutate(treat=gsub(".*-ATAC-|_.*", "", NEW_BARCODE)) 
meta <- combined@meta.data

meta$treat <- meta$treats

table(meta$wsnn_res.0.1)

#colnames(meta)
#table(meta$wsnn_res.0.05)

#----cell labeling----
# cell <- c("0"="Tcell", "1"="Tcell", "2"="NKcell", "3"="Bcell",
# "4"="Tcell", "5"="Tcell", "6"="Monocyte", "7"="Tcell", "8"="Tcell",
# "9"="Tcell", "10"="Tcell", "11"="DC", "12"="DC")

## cell <- c("0"="CD4Naive", "1"="TCM", "2"="NKcell", "3"="Bcell",
##           "4"="CD8Naive", "5"="TEM", "6"="Monocyte", "7"="CTL", "8"="dnT",
##           "9"="Treg", "10"="Platelet", "11"="DC", "12"="HSPC", "13"="Bcell")


## cell <- c("0"="CD4Naive", "1"="TCM", "2"="NKcell", "3"="Bcell",
##           "4"="CD8Naive", "5"="TEM", "6"="Monocyte", "7"="CTL", "8"="dnT",
##           "9"="Treg", "10"="Platelet", "11"="DC")

cell <- c("0"="CD4Naive", "1"="TCM", "2"="NKcell", "3"="TEM",
          "4"="Bcell", "5"="CD8Naive", "6"="Monocyte", "7"="dnT", "8"="MAIT",
          "9"="Platelet", "10"="DC")

#cell <- c("0"="Tcell", "1"="Tcell", "2"="Tcell", "3"="NKcell",
#         "4"="Bcell", "5"="Monocyte", "6"="Platelet", "7"="DC")


#meta$MCls <-NULL
meta$MCls <- cell[as.character(meta$wsnn_res.0.1)]
colnames(meta)
#table(meta$MCls)
#----cell labeling----


meta <- meta%>%
  mutate(bti=paste(wsnn_res.0.1, MCls, treat, SNG.BEST.GUESS, sep="_"))%>%
  dplyr::select(NEW_BARCODE, bti)
#colnames(meta)
head(meta)

dd <- meta%>%group_by(bti)%>%summarise(ncell=n(),.groups="drop")
dd

#write_rds(dd, "./1_DiffRNA_2.outs/0_ncell_clusters.rds")

#dd <- read_rds("./1_DiffRNA_2.outs/0_ncell_clusters.rds")
#dd


#----------------------------------------------------------------------------
count <- combined@assays$RNA@counts
#head(count)
#nrow(count)

count_genes <- rownames(count)
#head(count_genes)
#length(count_genes)





# filtering X and Y
#-------------------------------------------------------------------
library(biomaRt)

mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
#mart

genes.table <- try(biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description", "gene_biotype", "chromosome_name", "start_position"), mart = mart, useCache = F))
#colnames(genes.table)
#table(genes.table$chromosome_name)
#nrow(genes.table)

genes.table.2 <- genes.table%>%dplyr::filter(chromosome_name == c("X", "Y"), external_gene_name != "", external_gene_name %in% count_genes)
#table(genes.table.2$chromosome_name)
#colnames(genes.table.2)
#nrow(genes.table.2)

no_genes <- genes.table.2$external_gene_name
length(no_genes)

#check <- no_genes %in% count_genes
#check

#--------------------------------------------------------------
yes_genes <- count_genes[!(count_genes %in% no_genes)]

#length(yes_genes)

count.2 <- count[yes_genes,]
#head(count.2)
#head(rownames(count.2))
nrow(count.2)
ncol(count.2)

#---------------------------------------------------------------
anno <- data.frame(rn=rownames(count.2), rnz=rowSums(count.2))
#str(anno)
head(anno)

#autosome <- paste("chr", as.character(1:22), sep="")
#autosome

annoSel <- anno%>%dplyr::filter(rnz>0)
head(annoSel)

Y <- count[annoSel$rn,]
head(Y)


#------------------------------------------------------------------
##pseudo-bulk peak data
bti <- factor(meta$bti)
sort(unique(bti))

X <- model.matrix(~0+bti)
#head(X)

YtX <- Y %*% X
YtX <- as.matrix(YtX)
head(YtX)

nrow(YtX)
ncol(YtX)

colnames(YtX) <- gsub("^bti", "", colnames(YtX))
#colnames(YtX)

###rnz>0,chr:1-22, 111,750*400
#opfn <- "./1_DiffRNA_2.outs/1_YtX.comb_clusters.rds" 
#write_rds(YtX, file=opfn)


#-------------------------------------------------------------------
###rnz>100, ncell>20, 111,746 * 331
#dd <- read_rds("./1_DiffRNA_2.outs/0_ncell_clusters.rds")
#dd <- read_rds("./1_DiffRNA_2.outs/0_ncell_clusters.rds")%>%filter(ncell>20)

dd <- dd %>%filter(ncell>20)
#dd
#table(dd$bti)

## YtX <- read_rds("./1_DiffRNA_2.outs/1_YtX.comb_clusters.rds")
## head(YtX)
## ncol(YtX)
## nrow(YtX)

anno <- data.frame(rn=rownames(YtX), rnz=rowSums(YtX))
head(anno)

annoSel <- anno%>%filter(rnz>100)
#head(annoSel)

YtX_sel <- YtX[annoSel$rn, dd$bti]
head(YtX_sel)

nrow(YtX_sel)
ncol(YtX_sel)

opfn <- paste0("./1_DiffRNA_2.outs/1_YtX.sel_clusters_", reso, "_cn.rds")
write_rds(YtX_sel, file=opfn)

YtX_sel <- read_rds("./1_DiffRNA_2.outs/1_YtX.sel_clusters_0.1_cn.rds")

#####################################################
### 2. Differential analysis for CRP, High vs Low ###
#####################################################
### function, adjusted p value
## myqval <- function(pval){
##   qval <- pval
##   ii0 <- !is.na(pval)
##   qval[ii0] <- qvalue(pval[ii0],pi0=1)$qvalues
##   qval
## }


#### Read data
#YtX <- read_rds("./1_DiffRNA_2.outs/1_YtX.sel_clusters.rds")

YtX  <- YtX_sel
#ncol(YtX)
#nrow(YtX)
#str(YtX)
#colnames(YtX)

bti2 <- colnames(YtX)
#length(bti2)

cvt0 <- str_split(bti2, "_", simplify=T)
#head(cvt0)

cvt <- data.frame(bti=bti2, MCls=paste(cvt0[,1], cvt0[,2], sep="_"), treat=cvt0[,3], sampleID=cvt0[,4])
#cvt

#comb <- table(cvt$MCls, cvt$treat)
#write.csv(comb, paste0("./1_DiffRNA_2.outs/combinations.csv"), quote = FALSE, row.names = TRUE)
#comb <- read.csv("./1_DiffRNA_2.outs/combinations.csv")
#comb
#dd2 <- dd %>% dplyr::filter(bti %in% bti2)
#sum(dd2$ncell)




# filtered
#-------------------------------------------------------------------------
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
#comb.2

write.csv(comb.2, paste0("./1_DiffRNA_2.outs/combinations_2_", reso, "_cn.csv"), quote = FALSE, row.names = TRUE)

#comb.2 <- read.csv("./1_DiffRNA_2.outs/combinations_2.csv")
#comb.2

bti3 <- cvt.2$bti 
#bti3

#ncol(YtX)
YtX <- YtX[,bti3]
#ncol(YtX)

dd2 <- dd %>% dplyr::filter(bti %in% bti3)
#dd2
#sum(dd2$ncell)





#################################################################
### 2.1 call DESeq ###
#################################################################
contrast.list <- list("caffeine"=c("treat", "caffeine", "water"),
                      "nicotine"=c("treat", "nicotine", "water"),
                      "vitA"=c("treat", "vitA", "etOH"),
                      "vitD"=c("treat", "vitD", "etOH"),
                      "vitE"=c("treat", "vitE", "etOH"),
                      "water"=c("treat", "water", "etOH"),
                      "zinc"=c("treat", "zinc", "water"))

cat("2.1", "run DESeq batch separately", "\n")


#MCls <- c("CD4Naive", "TCM", "NKcell", "Bcell", "CD8Naive", "TEM", "Monocyte", "CTL", "dnT")

MCls <- unique(cvt.2$MCls)
MCls


#----INDIVIDUAL STEP----
MCls <- "0_CD4Naive"

treatment=c("caffeine"="red", "nicotine"="tan",
            "vitA"="tan4", "vitD"="seagreen4",
            "vitE"="salmon3", "water"="grey", "zinc"="maroon3")


cvt0 <- cvt.2%>%dplyr::filter(MCls=="0_CD4Naive")
YtX0 <- YtX[,cvt0$bti]
##

dds <- DESeqDataSetFromMatrix(YtX0, cvt0, ~ sampleID + treat) #also controlling individual effect

vd <- vst(dds)

treatment=c("caffeine"="red", "nicotine"="tan",
            "vitA"="tan4", "vitD"="seagreen4",
            "vitE"="salmon3", "water"="grey",
            "zinc"="maroon3", "etOH"="grey")

cols <- c("red", "grey", "tan", "tan4", "seagreen4", "salmon3", "grey", "maroon3")

color_key <- c("red", "tan", "tan4", "seagreen4", "salmon3", "grey", "maroon3")
names(color_key) <- c("caffeine", "nicotine", "vitA", "vitD", "vitE", "water", "zinc")

data <- plotPCA(vd, intgroup=c("treat"), returnData=TRUE)
head(data)

data <- data %>% mutate(sampleID=str_split(data$name, "_", simplify=T)[,4])
head(data)

#p <- qplot(PC1,PC2,color=treat,data=data,size=I(10))

#p <- ggplot(data, aes(x=PC1, y=PC2)) + geom_point(aes(colour = factor(treat)))
#p + scale_colour_manual(values = treatment)
#p + scale_fill_manual(values = color_key[names(color_key) %in% data$treat])

p <- ggplot(data, aes(x=PC1, y=PC2, colour = factor(treat))) +
  geom_point() +
  scale_colour_manual(values = treatment)+
  ggtitle("")
ggsave("./1_DiffRNA_2.outs/PCA_condition.png",p)

png(paste0("./1_DiffRNA_2.outs/plotpca_treatandind_filt_", reso, "_", oneX, "_cn.png"), width=2000, height=2000, pointsize=14, res=175)
print(p)
dev.off()

rlang::last_error()


#----FUNCTION----
res.DE <- map_dfr(MCls, function(oneX){
  time0 <- Sys.time()
  ###
  cvt0 <- cvt.2%>%dplyr::filter(MCls==oneX)
  YtX0 <- YtX[,cvt0$bti]
  ##
  dds <- DESeqDataSetFromMatrix(YtX0, cvt0, ~ sampleID + treat) #also controlling individual effect
  vd <- vst(dds)
  #p <- plotPCA(vd, intgroup=c("treat"))
  #png(paste0("./1_DiffRNA_2.outs/plotpca_treatandind_filt_", reso, "_", oneX, "_cn.png"), width=2000, height=2000, pointsize=14, res=175)
  #print(p)
  #dev.off()
  data <- plotPCA(vd, intgroup=c("treat"), returnData=TRUE)
  #p <- qplot(PC1,PC2,color=treatment,data=data,size=I(10))
  p <- ggplot(data, aes(x=PC1, y=PC2, colour = factor(treat))) +
    geom_point() +
    scale_colour_manual(values = treatment)+
    ggtitle(oneX)
  ggsave(paste0("./1_DiffRNA_2.outs/plotpca_treatandind_filt_", reso, "_", oneX, "_cn.png"),p)
  dds <- DESeq(dds)
  res <- contrast.list%>%map(~results(dds, contrast=.x))
  res2 <- tibble(contrast=names(res),MCls=oneX, data=map(res,tidy))%>%unnest(data)
  ##
  time1 <- Sys.time()
  elapsed <- difftime(time1, time0, units="mins")
  cat(oneX, elapsed, "Done\n")
  res2
})
## head(res.DE)
## table(res.DE$contrast)
## table(res.DE$MCls)



opfn <- paste0("./1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_", reso, "_cn.rds")
write_rds(res.DE, opfn)

res.DE <- readRDS("./1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn.rds")

head(res.DE)

#nrow(res.DE)
#nrow(res.DE %>% drop_na(p.adjusted))
#nrow(res.DE %>% drop_na(p.value))

