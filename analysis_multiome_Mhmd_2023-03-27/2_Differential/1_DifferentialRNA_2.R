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
reso = 0.05 

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
#colnames(meta)

# cell <- c("0"="Tcell", "1"="Tcell", "2"="NKcell", "3"="Bcell",
# "4"="Tcell", "5"="Tcell", "6"="Monocyte", "7"="Tcell", "8"="Tcell",
# "9"="Tcell", "10"="Tcell", "11"="DC", "12"="DC")

## cell <- c("0"="CD4Naive", "1"="TCM", "2"="NKcell", "3"="Bcell",
##           "4"="CD8Naive", "5"="TEM", "6"="Monocyte", "7"="CTL", "8"="dnT",
##           "9"="Treg", "10"="Platelet", "11"="DC", "12"="HSPC")

## cell <- c("0"="CD4Naive", "1"="TCM", "2"="NKcell", "3"="TEM",
##           "4"="Bcell", "5"="CD8Naive", "6"="Monocyte", "7"="dnT", "8"="MAIT",
##           "9"="Platelet", "10"="DC")

cell <- c("0"="Tcell", "1"="Tcell", "2"="Tcell", "3"="NKcell",
          "4"="Bcell", "5"="Monocyte", "6"="Platelet", "7"="DC")

#meta$MCls <- NULL
#meta$MCls <- cell[as.character(meta$seurat_cluster)]
meta$MCls <- cell[as.character(meta$wsnn_res.0.05)]
#meta$MCls <- cell[as.character(paste0("meta$wsnn_res.", "0.1"))]
#colnames(meta)
#table(meta$MCls)

meta <- meta%>%
   mutate(bti=paste(MCls, treat, SNG.BEST.GUESS, sep="_"))%>%
   dplyr::select(NEW_BARCODE, bti)

colnames(meta)

dd <- meta%>%group_by(bti)%>%summarise(ncell=n(),.groups="drop")

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
#length(no_genes)

#check <- no_genes %in% count_genes
#check

#--------------------------------------------------------------
yes_genes <- count_genes[!(count_genes %in% no_genes)]

# length(yes_genes)

count.2 <- count[yes_genes,]
#head(count.2)
#head(rownames(count.2))
#nrow(count.2)


#---------------------------------------------------------------
anno <- data.frame(rn=rownames(count.2), rnz=rowSums(count.2))
#str(anno)
#head(anno)

#autosome <- paste("chr", as.character(1:22), sep="")
#autosome

annoSel <- anno%>%dplyr::filter(rnz>0)
#head(annoSel)

Y <- count[annoSel$rn,]
#head(Y)


#------------------------------------------------------------------
##pseudo-bulk peak data
bti <- factor(meta$bti)
#bti

X <- model.matrix(~0+bti)
#head(X)

YtX <- Y %*% X

YtX <- as.matrix(YtX)
#head(YtX)

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
#head(anno)

annoSel <- anno%>%filter(rnz>100)
#head(annoSel)

YtX_sel <- YtX[annoSel$rn, dd$bti]
#head(YtX_sel)
#ncol(YtX_sel)
#nrow(YtX_sel)

opfn <- paste0("./1_DiffRNA_2.outs/1_YtX.sel_clusters_", reso, ".rds")
write_rds(YtX_sel, file=opfn)



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

cvt <- data.frame(bti=bti2, MCls=cvt0[,1], treat=cvt0[,2], sampleID=cvt0[,3])

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

comb.2 <- table(cvt.2$MCls, cvt.2$treat)
#comb.2

write.csv(comb.2, paste0("./1_DiffRNA_2.outs/combinations_2_", reso, ".csv"), quote = FALSE, row.names = TRUE)

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
#MCls

#MCls = "Bcell"

## cvt0 <- cvt.2%>%dplyr::filter(MCls=="Bcell")
## YtX0 <- YtX[,cvt0$bti]
## ##
## dds <- DESeqDataSetFromMatrix(YtX0, cvt0, ~ sampleID + treat) #also controlling individual effect
## #dds
## #rld <- rlog(dds)
## vd <- vst(dds)
## p <- plotPCA(vd,  intgroup=c("treat"))
## png(paste0("./1_DiffRNA_2.outs/plotpca_treatandind_filt_", reso, ".png"), width=2000, height=2000, pointsize=14, res=175)
## print(p)
## dev.off()
## dds <- DESeq(dds)
## res <- contrast.list%>%map(~results(dds, contrast=.x))
## res2 <- tibble(contrast=names(res),MCls=oneX, data=map(res,tidy))%>%unnest(data)  


res.DE <- map_dfr(MCls, function(oneX){
  time0 <- Sys.time()
###
  cvt0 <- cvt.2%>%dplyr::filter(MCls==oneX)
  YtX0 <- YtX[,cvt0$bti]
##
  dds <- DESeqDataSetFromMatrix(YtX0, cvt0, ~ sampleID + treat) #also controlling individual effect
  vd <- vst(dds)
  p <- plotPCA(vd, intgroup=c("treat"))
  png(paste0("./1_DiffRNA_2.outs/plotpca_treatandind_filt_", reso, "_", oneX, ".png"), width=2000, height=2000, pointsize=14, res=175)
  print(p)
  dev.off()
  dds <- DESeq(dds)
  res <- contrast.list%>%map(~results(dds, contrast=.x))
  res2 <- tibble(contrast=names(res),MCls=oneX, data=map(res,tidy))%>%unnest(data)  
##
  time1 <- Sys.time()
  elapsed <- difftime(time1, time0, units="mins")
  cat(oneX, elapsed, "Done\n")
  res2
})
#head(res.DE)
#table(res$contrast)
#table(res$MCls)

opfn <- paste0("./1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_", reso, ".rds")
write_rds(res.DE, opfn)

#res.DE <- readRDS("./1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1.rds")

#res.DE
#nrow(res.DE)
#nrow(res.DE %>% drop_na(p.adjusted))
#nrow(res.DE %>% drop_na(p.value))



## ######################################################
## # significant genes
## ######################################################
# res <- res.DE
## padj_cutoff <- 0.05

## sig_res <- dplyr::filter(res, p.adjusted < padj_cutoff) %>% dplyr::arrange(p.adjusted)

## #write.csv(sig_res, paste0("./1_DiffRNA_2.outs/sig_genes.csv"), quote = FALSE, row.names = FALSE)

## #sig_res <- read.csv("./1_DiffRNA_2.outs/sig_genes.csv")

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
res <- res.DE
FDR <- 0.1
pos_foldchange <- 0.25
neg_foldchange <- -0.25

nrow(res %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=FDR)  %>% dplyr::filter(estimate >=pos_foldchange))
nrow(res %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=FDR) %>% dplyr::filter(estimate <=neg_foldchange))

res[!is.na(res$p.adjusted) & res$p.adjusted<=FDR & res$estimate>=pos_foldchange, "Significance"] <- "Upregulated"
res[!is.na(res$p.adjusted) & res$p.adjusted<=FDR & res$estimate<=neg_foldchange, "Significance"] <- "Downregulated"
res[!is.na(res$p.adjusted) & res$p.adjusted>FDR,"Significance"] <- "Not Significant"
res[is.na(res$p.adjusted),"Significance"] <- "NA"
#res2 <- res %>% dplyr::filter(Significance %in% c("Upregulated", "Downregulated"))
#res

#---------------------------
res_up <- res %>% dplyr::filter(Significance %in% "Upregulated")
#nrow(res_up)

## res_up.uniquegenes <-unique(res_up$gene)
## length(res_up.uniquegenes)
## #res_up.uniquegenes

## res_up.2 <- res_up %>% dplyr::filter(gene %in% res_up.uniquegenes)
## nrow(res_up.2)

MCls <- unique(res_up$MCls)
#MCls
res_up.unique <- map_dfr(MCls, function(oneX){
      time0 <- Sys.time()
        ###
        res_up.2 <- res_up%>%dplyr::filter(MCls==oneX)
        ##
        res_up.3 <- res_up%>%dplyr::filter(MCls!=oneX)
            #mutate(count=ifelse(nrow(cvt.2)>=8,1,0)) 
        ##
        other.genes <- unique(res_up.3$gene)
        res_up.2.2 <- res_up.2%>%dplyr::filter(!gene %in% other.genes)
        time1 <- Sys.time()
        elapsed <- difftime(time1, time0, units="mins")
        cat(oneX, elapsed, "Done\n")
        res_up.2.2
      })
#res_up.unique
#nrow(res_up.unique)

up <- table(res_up.unique$MCls, res_up.unique$contrast)
#up
write.csv(up, paste0("./1_DiffRNA_2.outs/up_treat&ind_filt_unique_", reso, "_", FDR, "_", pos_foldchange, ".csv"), quote = FALSE, row.names = TRUE)



#---------------------------
res_down <- res %>% dplyr::filter(Significance %in% "Downregulated")

MCls <- unique(res_down$MCls)
#MCls
res_down.unique <- map_dfr(MCls, function(oneX){
      time0 <- Sys.time()
        ###
        res_down.2 <- res_down%>%dplyr::filter(MCls==oneX)
        ##
        res_down.3 <- res_down%>%dplyr::filter(MCls!=oneX)
            #mutate(count=ifelse(nrow(cvt.2)>=8,1,0)) 
        ##
        other.genes <- unique(res_down.3$gene)
        res_down.2.2 <- res_down.2%>%dplyr::filter(!gene %in% other.genes)
        time1 <- Sys.time()
        elapsed <- difftime(time1, time0, units="mins")
        cat(oneX, elapsed, "Done\n")
        res_down.2.2
      })
#res_down.unique
#nrow(res_down.unique)


down <- table(res_down.unique$MCls, res_down.unique$contrast)
#down
write.csv(down, paste0("./1_DiffRNA_2.outs/down_treat&ind_filt_unique_", reso, "_", FDR, "_", pos_foldchange, ".csv"), quote = FALSE, row.names = TRUE)

#up <- read.csv(file = './1_DiffRNA_2.outs/up_treat&ind_filt.csv')
#up
#down <- read.csv(file = './1_DiffRNA_2.outs/down_treat&ind_filt.csv')
#down

# #-----------number of cells------------------------------------------------
# #cvt.2
# 
# dd3 <- left_join(cvt.2, dd2)
# #head(dd3)
# 
# dd3 <- dd3 %>% dplyr::filter(dd3$MCls %in% res$MCls, !dd3$treat %in% "etOH")
# #head(dd3)
# #nrow(dd3)
# #sum(dd3$ncell)
# 
# dd4 <- dd3 %>% group_by(treat, MCls) %>% summarise(ncell_count = sum(ncell))
# #dd4
# 
# ncell_table <- data.frame(pivot_wider(dd4, names_from=treat, values_from=ncell_count))
# 
# rownames(ncell_table) <- ncell_table$MCls
# ncell_table$MCls <- NULL
# #ncell_table
# 
# write.csv(ncell_table, paste0("./1_DiffRNA_2.outs/ncell_treat$ind_filt_", reso, "_", FDR, "_", pos_foldchange, ".csv"), quote = FALSE, row.names = TRUE)
# 
# #ncell_table <- read.csv("./1_DiffRNA_2.outs/ncell_treat$ind_filt.csv")


#########################################################
# matrix_plots
#########################################################
library(reshape2)

up.matrix <-as.matrix(up) # create a numeric matrix object
up.melt <- melt(up.matrix)

#down
#rownames(down)<- down$X
#down$X <- NULL
down.matrix <- as.matrix(down)
down.melt <- melt(down.matrix)

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
    geom_tile(aes(fill=value))+
    geom_text(aes(label = value), color = "white", size = 4) +
    ggtitle("Upregulated") +
    theme(axis.text.x = element_text(angle = 90))

down.p <- ggplot(down.melt, aes(x = Var2, y = Var1)) +
      geom_tile(aes(fill=value))+
      geom_text(aes(label = value), color = "white", size = 4) +
      ggtitle("Downregulated") +
      theme(axis.text.x = element_text(angle = 90))

ncell.p <- ggplot(ncell.melt, aes(x = Var2, y = Var1)) +
      geom_tile(aes(fill=value))+
      geom_text(aes(label = value), color = "white", size = 4) +
      ggtitle("ncells") +
      theme(axis.text.x = element_text(angle = 90))

comb.3.p <- ggplot(comb.3.melt, aes(x = Var2, y = Var1)) +
      geom_tile(aes(fill=value))+
      geom_text(aes(label = value), color = "white", size = 4) +
      ggtitle("Combinations") +
      theme(axis.text.x = element_text(angle = 90))

library("ggpubr")
figure <- ggarrange(comb.3.p, ncell.p, up.p, down.p +
                    font("x.text", size = 12),
                    ncol = 2, nrow = 2)

png(paste0("./1_DiffRNA_2.outs/plot_all_treatandind_filt_unique_", reso, "_", FDR, "_", pos_foldchange, ".png"), width=2000, height=2000, pointsize=14, res=175)
print(figure)
dev.off()


##############################################################
###
### 1. MA plots
##############################################################
figfn <- paste0("./1_DiffRNA_2.outs/Figure1.1_MA_clusters_treat&ind_filt_", reso, "_", FDR, "_", pos_foldchange, ".png")
png(figfn, width=3500, height=4500, pointsize=12, res=225)
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




### 2, qq plots
figfn <- paste0("./1_DiffRNA_2.outs/Figure1.2_qq_clusters_treat&ind_filt_", reso, "_", FDR, "_", pos_foldchange, ".png")
png(figfn, width=3500, height=4500, pointsize=16, res=225)
par(mar=c(4,4,2,2),mgp=c(2,1,0))
x <- matrix(1:63, 9, 7, byrow=T)
layout(x)
#fn <- "./1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered.rds"
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
png(paste0("./1_DiffRNA_2.outs/Figure1.3_qq_clusters_together_treatandind_filt_", reso,"_", FDR, "_", pos_foldchange, ".png"), width=2000, height=5000, pointsize=16, res=200)
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
           cex = 1,
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
res <- res.DE
FDR <- 0.05
pos_foldchange <- 0.00
neg_foldchange <- 0.00

nrow(res %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=FDR)  %>% dplyr::filter(estimate >=pos_foldchange))
nrow(res %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=FDR) %>% dplyr::filter(estimate <neg_foldchange))

res[!is.na(res$p.adjusted) & res$p.adjusted<=FDR & res$estimate>pos_foldchange, "Significance"] <- "Upregulated"
res[!is.na(res$p.adjusted) & res$p.adjusted<=FDR & res$estimate<=neg_foldchange, "Significance"] <- "Downregulated"
res[!is.na(res$p.adjusted) & res$p.adjusted>FDR,"Significance"] <- "Not Significant"
res[is.na(res$p.adjusted),"Significance"] <- "NA"
#res2 <- res %>% dplyr::filter(Significance %in% c("Upregulated", "Downregulated"))

#---------------------------
res_up <- res %>% dplyr::filter(Significance %in% "Upregulated")
#nrow(res_up)

## res_up.uniquegenes <-unique(res_up$gene)
## length(res_up.uniquegenes)
## #res_up.uniquegenes

## res_up.2 <- res_up %>% dplyr::filter(gene %in% res_up.uniquegenes)
## nrow(res_up.2)

MCls <- unique(res_up$MCls)
#MCls
res_up.unique <- map_dfr(MCls, function(oneX){
  time0 <- Sys.time()
  ###
  res_up.2 <- res_up%>%dplyr::filter(MCls==oneX)
  ##
  res_up.3 <- res_up%>%dplyr::filter(MCls!=oneX)
  #mutate(count=ifelse(nrow(cvt.2)>=8,1,0)) 
  ##
  other.genes <- unique(res_up.3$gene)
  res_up.2.2 <- res_up.2%>%dplyr::filter(!gene %in% other.genes)
  time1 <- Sys.time()
  elapsed <- difftime(time1, time0, units="mins")
  cat(oneX, elapsed, "Done\n")
  res_up.2.2
})
#res_up.unique
#nrow(res_up.unique)

up <- table(res_up.unique$MCls, res_up.unique$contrast)
#up
write.csv(up, paste0("./1_DiffRNA_2.outs/up_treat&ind_filt_unique_", reso, "_", FDR, "_", pos_foldchange, ".csv"), quote = FALSE, row.names = TRUE)



#---------------------------
res_down <- res %>% dplyr::filter(Significance %in% "Downregulated")

MCls <- unique(res_down$MCls)
#MCls
res_down.unique <- map_dfr(MCls, function(oneX){
  time0 <- Sys.time()
  ###
  res_down.2 <- res_down%>%dplyr::filter(MCls==oneX)
  ##
  res_down.3 <- res_down%>%dplyr::filter(MCls!=oneX)
  #mutate(count=ifelse(nrow(cvt.2)>=8,1,0)) 
  ##
  other.genes <- unique(res_down.3$gene)
  res_down.2.2 <- res_down.2%>%dplyr::filter(!gene %in% other.genes)
  time1 <- Sys.time()
  elapsed <- difftime(time1, time0, units="mins")
  cat(oneX, elapsed, "Done\n")
  res_down.2.2
})
#res_down.unique
#nrow(res_down.unique)


down <- table(res_down.unique$MCls, res_down.unique$contrast)
#down
write.csv(down, paste0("./1_DiffRNA_2.outs/down_treat&ind_filt_unique_", reso, "_", FDR, "_", pos_foldchange, ".csv"), quote = FALSE, row.names = TRUE)


# res_up <- res %>% dplyr::filter(Significance %in% "Upregulated")
# res_down <- res %>% dplyr::filter(Significance %in% "Downregulated")
# 
# up <- table(res_up$MCls, res_up$contrast)
# #up
# write.csv(up, paste0("./1_DiffRNA_2.outs/up_treat&ind_filt_", reso, "_", FDR, "_", pos_foldchange, ".csv"), quote = FALSE, row.names = TRUE)
# 
# down <- table(res_down$MCls, res_down$contrast)
# #down
# write.csv(down, paste0("./1_DiffRNA_2.outs/down_treat&ind_filt_", reso, "_", FDR, "_", pos_foldchange, ".csv"), quote = FALSE, row.names = TRUE)
# 
# #up <- read.csv(file = './1_DiffRNA_2.outs/up_treat&ind_filt.csv')
# #up
# #down <- read.csv(file = './1_DiffRNA_2.outs/down_treat&ind_filt.csv')
# #down

#-----------number of cells------------------------------------------------
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

write.csv(ncell_table, paste0("./1_DiffRNA_2.outs/ncell_treat$ind_filt_", reso, "_", FDR, "_", pos_foldchange, ".csv"), quote = FALSE, row.names = TRUE)

#ncell_table <- read.csv("./1_DiffRNA_2.outs/ncell_treat$ind_filt.csv")


#########################################################
# matrix_plots
#########################################################
#library(reshape2)

up.matrix <-as.matrix(up) # create a numeric matrix object
up.melt <- melt(up.matrix)

#down
#rownames(down)<- down$X
#down$X <- NULL
down.matrix <- as.matrix(down)
down.melt <- melt(down.matrix)

#ncell
#rownames(ncell)<- ncell$X
#ncell$X <- NULL
ncell.matrix <- as.matrix(ncell_table)
ncell.melt <- melt(ncell.matrix)

#comb
#rownames(comb)<- comb$X
#comb$X <- NULL
#comb <- comb.2
#comb.2 <- subset(comb, select = -c(etOH))
#comb.2 <- comb.2 %>% dplyr::filter(rownames(comb.2) %in% rownames(up))
comb.3.matrix <- as.matrix(comb.2)
comb.3.matrix <- subset.matrix(comb.3.matrix, select = -c(etOH))
comb.3.melt <- melt(comb.3.matrix)

up.p <- ggplot(up.melt, aes(x = Var2, y = Var1)) +
    geom_tile(aes(fill=value))+
    geom_text(aes(label = value), color = "white", size = 4) +
    ggtitle("Upregulated") +
    theme(axis.text.x = element_text(angle = 90))

down.p <- ggplot(down.melt, aes(x = Var2, y = Var1)) +
      geom_tile(aes(fill=value))+
      geom_text(aes(label = value), color = "white", size = 4) +
      ggtitle("Downregulated") +
      theme(axis.text.x = element_text(angle = 90))

ncell.p <- ggplot(ncell.melt, aes(x = Var2, y = Var1)) +
      geom_tile(aes(fill=value))+
      geom_text(aes(label = value), color = "white", size = 4) +
      ggtitle("ncells") +
      theme(axis.text.x = element_text(angle = 90))

comb.3.p <- ggplot(comb.3.melt, aes(x = Var2, y = Var1)) +
      geom_tile(aes(fill=value))+
      geom_text(aes(label = value), color = "white", size = 4) +
      ggtitle("Combinations") +
      theme(axis.text.x = element_text(angle = 90))

library("ggpubr")
figure <- ggarrange(comb.3.p, ncell.p, up.p, down.p +
                    font("x.text", size = 12),
                    ncol = 2, nrow = 2)

png(paste0("./1_DiffRNA_2.outs/plot_all_treatandind_filt_unique_", reso, "_", FDR, "_", pos_foldchange, ".png"), width=2000, height=2000, pointsize=14, res=175)
print(figure)
dev.off()


##############################################################
###
### 1. MA plots
##############################################################
figfn <- paste0("./1_DiffRNA_2.outs/Figure1.1_MA_clusters_treat&ind_filt_", reso, "_", FDR, "_", pos_foldchange, ".png")
png(figfn, width=3500, height=4500, pointsize=12, res=225)
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




## ### 2, qq plots
## figfn <- paste0("./1_DiffRNA_2.outs/Figure1.2_qq_clusters_treat&ind_filt_", reso, "_", FDR, "_", pos_foldchange, ".png")
## png(figfn, width=3500, height=4500, pointsize=16, res=225)
## par(mar=c(4,4,2,2),mgp=c(2,1,0))
## x <- matrix(1:63, 9, 7, byrow=T)
## layout(x)
## #fn <- "./1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered.rds"
## #res <- read_rds(fn)%>%drop_na(p.value)
## res.1 <- res%>%drop_na(p.value)
## #MCls <- c("CD4Naive", "TCM", "NKcell", "Bcell", "CD8Naive", "TEM", "Monocyte", "CTL", "dnT")
## MCls <- unique(cvt.2$MCls)
## for (oneMCl in MCls){
##   ##1
##   res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="caffeine")
##   print( qq(res.2$p.value, main="caffeine", cex.main=1, cex.axis=0.8, cex.lab=1))
##   ##2
##   res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="nicotine")
##   print( qq(res.2$p.value, main="nicotine", cex.main=1, cex.axis=0.8, cex.lab=1))
##   ##3 
##   res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="vitA")
##   print( qq(res.2$p.value, main="vitA", cex.main=1, cex.axis=0.8, cex.lab=1))
##   ##4
##   res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="vitD")
##   print( qq(res.2$p.value, main="vitD", cex.main=1, cex.axis=0.8, cex.lab=1))
##   ##5
##   res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="vitE")
##   print( qq(res.2$p.value, main="vitE", cex.main=1, cex.axis=0.8, cex.lab=1))
##   ##6
##   res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="water")
##   print( qq(res.2$p.value, main="water", cex.main=1, cex.axis=0.8, cex.lab=1))
##   ##7
##   res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="zinc")
##   print( qq(res.2$p.value, main="zinc", cex.main=1, cex.axis=0.8, cex.lab=1))   
##   print( mtext(oneMCl, side=4, line=0.5, cex=1, col="blue"))  
## }
## dev.off()


## ###
## ### 3, canno plots


## #------------------plot_loop--------------------------------
## #nrow(res) #if qq plot has been called before this step, this already is with dropped NA values
## res.1 <- res %>% drop_na(p.value)
## max.y <- floor(max(-log10(res.1$p.value)))/2
## #max.y
## png(paste0("./1_DiffRNA_2.outs/Figure1.3_qq_clusters_together_treatandind_filt_", reso,"_", FDR, "_", pos_foldchange, ".png"), width=2000, height=5000, pointsize=16, res=200)
## par(mar=c(4,4,2,2),mgp=c(2,1,0))
## x <- matrix(1:10, 5, 2, byrow=T)
## layout(x)
## # MCls <- c("CD4Naive", "TCM", "NKcell", "Bcell", "CD8Naive", "TEM", "Monocyte", "CTL", "dnT")
## MCls <- unique(cvt.2$MCls)
## for (oneMCl in MCls){
##   ##1
##   res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="caffeine")%>% drop_na(p.value)
##   plot(-log10(1:length(res.2$p.value)/length(res.2$p.value)),
##        -log10(sort(res.2$p.value)),
##        main=oneMCl,
##        pch = 19,
##        ylim = c(0, max.y),
##        ylab="observed -log10(p)",
##        xlab="expected -log10(p)",
##        col="red")
##   ##2
##   res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="nicotine")%>% drop_na(p.value)
##   points(-log10(1:length(res.2$p.value)/length(res.2$p.value)),-log10(sort(res.2$p.value)),  col="tan", pch = 19)   
##   ##3 
##   res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="zinc")%>% drop_na(p.value)
##   points(-log10(1:length(res.2$p.value)/length(res.2$p.value)),-log10(sort(res.2$p.value)),  col="maroon3", pch = 19)
##   ##4
##   res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="vitA")%>% drop_na(p.value)
##   points(-log10(1:length(res.2$p.value)/length(res.2$p.value)),-log10(sort(res.2$p.value)),  col="tan4", pch = 19)
##   ##5
##   res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="vitD")%>% drop_na(p.value)
##   points(-log10(1:length(res.2$p.value)/length(res.2$p.value)),-log10(sort(res.2$p.value)),  col="seagreen4", pch = 19)
##   ##6
##   res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="water")%>% drop_na(p.value)
##   points(-log10(1:length(res.2$p.value)/length(res.2$p.value)),-log10(sort(res.2$p.value)),  col="grey", pch = 19)
##   ##7
##   res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="vitE")%>% drop_na(p.value)
##   points(-log10(1:length(res.2$p.value)/length(res.2$p.value)),-log10(sort(res.2$p.value)),  col="salmon3", pch = 19)
##   ##control
##   abline(0, 1, col = "black")
##   ##legend
##   legend("topleft",
##          #inset=c(-0.2,0),
##          cex = 1,
##          pch = 19,
##          c("caff","nic", "zinc", "vitA", "vitD", "vitE", "water"),
##          fill=c("red","tan", "maroon3", "tan4", "seagreen4", "salmon3", "grey"))
## }
## dev.off()













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
##           "9"="Treg", "10"="Platelet", "11"="DC", "12"="HSPC")

## cell <- c("0"="CD4Naive", "1"="TCM", "2"="NKcell", "3"="TEM",
##           "4"="Bcell", "5"="CD8Naive", "6"="Monocyte", "7"="dnT",
##           "8"="MAIT","9"="Platelet", "10"="DC")

cell <- c("0"="Tcell", "1"="Tcell", "2"="Tcell", "3"="NKcell",
          "4"="Bcell", "5"="Monocyte", "6"="Platelet", "7"="DC")

meta$MCls <- cell[as.character(meta$wsnn_res.0.05)]

meta <- meta%>%mutate(bti=paste(MCls, treat, SNG.BEST.GUESS, sep="_"))%>%dplyr::select(NEW_BARCODE, bti)
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

opfn <- paste0("./1_DiffPeak.outs/1_YtX.sel_", reso, ".rds")
write_rds(YtX_sel, file=opfn)



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

cvt <- data.frame(bti=bti2, MCls=cvt0[,1], treat=cvt0[,2], sampleID=cvt0[,3])
# cvt
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
#comb.2

write.csv(comb.2, paste0("./1_DiffPeak.outs/combinations_2_", reso, ".csv"), quote = FALSE, row.names = TRUE)

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
#MCls

###

res.DE.ATAC <- map_dfr(MCls, function(oneX){
  time0 <- Sys.time()
  ###
  cvt0 <- cvt.2%>%dplyr::filter(MCls==oneX)
  YtX0 <- YtX[,cvt0$bti]
  ##
  dds <- DESeqDataSetFromMatrix(YtX0, cvt0, ~ sampleID + treat)
  dds <- DESeq(dds)
  #PCA plot color dot by treatment
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

opfn <- paste0("./1_DiffPeak.outs/2.0_DESeq.results_clusters_treatandind_filtered_", reso, ".rds")
write_rds(res.DE.ATAC, opfn)

res.DE.ATAC <- readRDS("./1_DiffPeak.outs/2.0_DESeq.results_clusters_treatandind_filtered_0.1.rds")

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


res[!is.na(res$p.adjusted) & res$p.adjusted<=FDR & res$estimate>=pos_foldchange, "Significance"] <- "Upregulated"
res[!is.na(res$p.adjusted) & res$p.adjusted<=FDR & res$estimate<=neg_foldchange, "Significance"] <- "Downregulated"
res[!is.na(res$p.adjusted) & res$p.adjusted>FDR,"Significance"] <- "Not Significant"
res[is.na(res$p.adjusted),"Significance"] <- "NA"
#res2 <- res %>% dplyr::filter(Significance %in% c("Upregulated", "Downregulated"))

#---------------------------
res_up <- res %>% dplyr::filter(Significance %in% "Upregulated")
#nrow(res_up)

## res_up.uniquegenes <-unique(res_up$gene)
## length(res_up.uniquegenes)
## #res_up.uniquegenes

## res_up.2 <- res_up %>% dplyr::filter(gene %in% res_up.uniquegenes)
## nrow(res_up.2)

MCls <- unique(res_up$MCls)
#MCls
res_up.unique <- map_dfr(MCls, function(oneX){
  time0 <- Sys.time()
  ###
  res_up.2 <- res_up%>%dplyr::filter(MCls==oneX)
  ##
  res_up.3 <- res_up%>%dplyr::filter(MCls!=oneX)
  #mutate(count=ifelse(nrow(cvt.2)>=8,1,0)) 
  ##
  other.genes <- unique(res_up.3$gene)
  res_up.2.2 <- res_up.2%>%dplyr::filter(!gene %in% other.genes)
  time1 <- Sys.time()
  elapsed <- difftime(time1, time0, units="mins")
  cat(oneX, elapsed, "Done\n")
  res_up.2.2
})
#res_up.unique
#nrow(res_up.unique)

up <- table(res_up.unique$MCls, res_up.unique$contrast)
#up
write.csv(up, paste0("./1_DiffRNA_2.outs/up_treat&ind_filt_unique_", reso, "_", FDR, "_", pos_foldchange, ".csv"), quote = FALSE, row.names = TRUE)



#---------------------------
res_down <- res %>% dplyr::filter(Significance %in% "Downregulated")

MCls <- unique(res_down$MCls)
#MCls
res_down.unique <- map_dfr(MCls, function(oneX){
  time0 <- Sys.time()
  ###
  res_down.2 <- res_down%>%dplyr::filter(MCls==oneX)
  ##
  res_down.3 <- res_down%>%dplyr::filter(MCls!=oneX)
  #mutate(count=ifelse(nrow(cvt.2)>=8,1,0)) 
  ##
  other.genes <- unique(res_down.3$gene)
  res_down.2.2 <- res_down.2%>%dplyr::filter(!gene %in% other.genes)
  time1 <- Sys.time()
  elapsed <- difftime(time1, time0, units="mins")
  cat(oneX, elapsed, "Done\n")
  res_down.2.2
})
#res_down.unique
#nrow(res_down.unique)


down <- table(res_down.unique$MCls, res_down.unique$contrast)
#down
write.csv(down, paste0("./1_DiffRNA_2.outs/down_treat&ind_filt_unique_", reso, "_", FDR, "_", pos_foldchange, ".csv"), quote = FALSE, row.names = TRUE)


# res_up <- res %>% dplyr::filter(Significance %in% "Upregulated")
# res_down <- res %>% dplyr::filter(Significance %in% "Downregulated")
# 
# up <- table(res_up$MCls, res_up$contrast)
# #up
# write.csv(up, paste0("./1_DiffPeak.outs/up_treat&ind_filt_", reso, "_", FDR, "_", pos_foldchange, ".csv"), quote = FALSE, row.names = TRUE)
# 
# down <- table(res_down$MCls, res_down$contrast)
# #down
# write.csv(down, paste0("./1_DiffPeak.outs/down_treat&ind_filt_", reso, "_", FDR, "_", pos_foldchange, ".csv"), quote = FALSE, row.names = TRUE)
# 
# #up <- read.csv(file = './1_DiffPeak.outs/up_treat&ind_filt.csv')
# #up
# #down <- read.csv(file = './1_DiffPeak.outs/down_treat&ind_filt.csv')
# #down

#-----------number of cells------------------------------------------------
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

write.csv(ncell_table, paste0("./1_DiffPeak.outs/ncell_treat$ind_filt_", reso, "_", FDR, "_", pos_foldchange, ".csv"), quote = FALSE, row.names = TRUE)

#ncell <- read.csv("./1_DiffPeak.outs/ncell_treat$ind_filt.csv")


#########################################################
# matrix_plots
#########################################################
library(reshape2)

up.matrix <-as.matrix(up) # create a numeric matrix object
up.melt <- melt(up.matrix)

#down
#rownames(down)<- down$X
#down$X <- NULL
down.matrix <- as.matrix(down)
down.melt <- melt(down.matrix)

#ncell
#rownames(ncell)<- ncell$X
#ncell$X <- NULL
ncell.matrix <- as.matrix(ncell_table)
ncell.melt <- melt(ncell.matrix)

#comb
#rownames(comb)<- comb$X
#comb$X <- NULL
#comb <- comb.2
#comb.2 <- comb.2 %>% dplyr::filter(rownames(comb.2) %in% rownames(up))
comb.3.matrix <- as.matrix(comb.2)
comb.3.matrix <- subset.matrix(comb.3.matrix, select = -c(etOH))
comb.3.melt <- melt(comb.3.matrix)

up.p <- ggplot(up.melt, aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill=value))+
  geom_text(aes(label = value), color = "white", size = 4) +
  ggtitle("Upregulated") +
  theme(axis.text.x = element_text(angle = 90))

down.p <- ggplot(down.melt, aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill=value))+
  geom_text(aes(label = value), color = "white", size = 4) +
  ggtitle("Downregulated") +
  theme(axis.text.x = element_text(angle = 90))

ncell.p <- ggplot(ncell.melt, aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill=value))+
  geom_text(aes(label = value), color = "white", size = 4) +
  ggtitle("ncells") +
  theme(axis.text.x = element_text(angle = 90))

comb.3.p <- ggplot(comb.3.melt, aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill=value))+
  geom_text(aes(label = value), color = "white", size = 4) +
  ggtitle("Combinations") +
  theme(axis.text.x = element_text(angle = 90))

library("ggpubr")
figure <- ggarrange(comb.3.p, ncell.p, up.p, down.p +
                      font("x.text", size = 12),
                    ncol = 2, nrow = 2)

png(paste0("./1_DiffPeak.outs/plot_all_treatandind_filt_unique_", reso, "_", FDR, "_", pos_foldchange, ".png"), width=2000, height=2000, pointsize=14, res=175)
print(figure)
dev.off()


##############################################################
###
### 1. MA plots
##############################################################
figfn <- paste0("./1_DiffPeak.outs/Figure1.1_MA_clusters_treat&ind_filt_", reso, "_", FDR, "_", pos_foldchange, ".png")
png(figfn, width=3500, height=4500, pointsize=12, res=225)
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
figfn <- paste0("./1_DiffPeak.outs/Figure1.2_qq_clusters_treat&ind_filt_", reso, "_", FDR, "_", pos_foldchange, ".png")
png(figfn, width=3500, height=4500, pointsize=16, res=225)
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
png(paste0("./1_DiffPeak.outs/Figure1.3_qq_clusters_together_treatandind_filt_", reso,"_", FDR, "_", pos_foldchange, ".png"), width=2000, height=5000, pointsize=16, res=200)
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
         cex = 1,
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


res[!is.na(res$p.adjusted) & res$p.adjusted<=FDR & res$estimate>pos_foldchange, "Significance"] <- "Upregulated"
res[!is.na(res$p.adjusted) & res$p.adjusted<=FDR & res$estimate<=neg_foldchange, "Significance"] <- "Downregulated"
res[!is.na(res$p.adjusted) & res$p.adjusted>FDR,"Significance"] <- "Not Significant"
res[is.na(res$p.adjusted),"Significance"] <- "NA"
#res2 <- res %>% dplyr::filter(Significance %in% c("Upregulated", "Downregulated"))

#---------------------------
res_up <- res %>% dplyr::filter(Significance %in% "Upregulated")
#nrow(res_up)

## res_up.uniquegenes <-unique(res_up$gene)
## length(res_up.uniquegenes)
## #res_up.uniquegenes

## res_up.2 <- res_up %>% dplyr::filter(gene %in% res_up.uniquegenes)
## nrow(res_up.2)

MCls <- unique(res_up$MCls)
#MCls
res_up.unique <- map_dfr(MCls, function(oneX){
  time0 <- Sys.time()
  ###
  res_up.2 <- res_up%>%dplyr::filter(MCls==oneX)
  ##
  res_up.3 <- res_up%>%dplyr::filter(MCls!=oneX)
  #mutate(count=ifelse(nrow(cvt.2)>=8,1,0)) 
  ##
  other.genes <- unique(res_up.3$gene)
  res_up.2.2 <- res_up.2%>%dplyr::filter(!gene %in% other.genes)
  time1 <- Sys.time()
  elapsed <- difftime(time1, time0, units="mins")
  cat(oneX, elapsed, "Done\n")
  res_up.2.2
})
#res_up.unique
#nrow(res_up.unique)

up <- table(res_up.unique$MCls, res_up.unique$contrast)
#up
write.csv(up, paste0("./1_DiffRNA_2.outs/up_treat&ind_filt_unique_", reso, "_", FDR, "_", pos_foldchange, ".csv"), quote = FALSE, row.names = TRUE)



#---------------------------
res_down <- res %>% dplyr::filter(Significance %in% "Downregulated")

MCls <- unique(res_down$MCls)
#MCls
res_down.unique <- map_dfr(MCls, function(oneX){
  time0 <- Sys.time()
  ###
  res_down.2 <- res_down%>%dplyr::filter(MCls==oneX)
  ##
  res_down.3 <- res_down%>%dplyr::filter(MCls!=oneX)
  #mutate(count=ifelse(nrow(cvt.2)>=8,1,0)) 
  ##
  other.genes <- unique(res_down.3$gene)
  res_down.2.2 <- res_down.2%>%dplyr::filter(!gene %in% other.genes)
  time1 <- Sys.time()
  elapsed <- difftime(time1, time0, units="mins")
  cat(oneX, elapsed, "Done\n")
  res_down.2.2
})
#res_down.unique
#nrow(res_down.unique)


down <- table(res_down.unique$MCls, res_down.unique$contrast)
#down
write.csv(down, paste0("./1_DiffRNA_2.outs/down_treat&ind_filt_unique_", reso, "_", FDR, "_", pos_foldchange, ".csv"), quote = FALSE, row.names = TRUE)


# res_up <- res %>% dplyr::filter(Significance %in% "Upregulated")
# res_down <- res %>% dplyr::filter(Significance %in% "Downregulated")
# 
# up <- table(res_up$MCls, res_up$contrast)
# #up
# write.csv(up, paste0("./1_DiffPeak.outs/up_treat&ind_filt_", reso, "_", FDR, "_", pos_foldchange, ".csv"), quote = FALSE, row.names = TRUE)
# 
# down <- table(res_down$MCls, res_down$contrast)
# #down
# write.csv(down, paste0("./1_DiffPeak.outs/down_treat&ind_filt_", reso, "_", FDR, "_", pos_foldchange, ".csv"), quote = FALSE, row.names = TRUE)
# 
# #up <- read.csv(file = './1_DiffPeak.outs/up_treat&ind_filt.csv')
# #up
# #down <- read.csv(file = './1_DiffPeak.outs/down_treat&ind_filt.csv')
# #down

#-----------number of cells------------------------------------------------
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

write.csv(ncell_table, paste0("./1_DiffPeak.outs/ncell_treat$ind_filt_", reso, "_", FDR, "_", pos_foldchange, ".csv"), quote = FALSE, row.names = TRUE)

#ncell <- read.csv("./1_DiffPeak.outs/ncell_treat$ind_filt.csv")


#########################################################
# matrix_plots
#########################################################
#library(reshape2)

up.matrix <-as.matrix(up) # create a numeric matrix object
up.melt <- melt(up.matrix)

#down
#rownames(down)<- down$X
#down$X <- NULL
down.matrix <- as.matrix(down)
down.melt <- melt(down.matrix)

#ncell
#rownames(ncell)<- ncell$X
#ncell$X <- NULL
ncell.matrix <- as.matrix(ncell_table)
ncell.melt <- melt(ncell.matrix)

#comb
#rownames(comb)<- comb$X
#comb$X <- NULL
#comb <- comb.2
#comb.2 <- comb.2 %>% dplyr::filter(rownames(comb.2) %in% rownames(up))
comb.3.matrix <- as.matrix(comb.2)
comb.3.matrix <- subset.matrix(comb.3.matrix, select = -c(etOH))
comb.3.melt <- melt(comb.3.matrix)

up.p <- ggplot(up.melt, aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill=value))+
  geom_text(aes(label = value), color = "white", size = 4) +
  ggtitle("Upregulated") +
  theme(axis.text.x = element_text(angle = 90))

down.p <- ggplot(down.melt, aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill=value))+
  geom_text(aes(label = value), color = "white", size = 4) +
  ggtitle("Downregulated") +
  theme(axis.text.x = element_text(angle = 90))

ncell.p <- ggplot(ncell.melt, aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill=value))+
  geom_text(aes(label = value), color = "white", size = 4) +
  ggtitle("ncells") +
  theme(axis.text.x = element_text(angle = 90))

comb.3.p <- ggplot(comb.3.melt, aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill=value))+
  geom_text(aes(label = value), color = "white", size = 4) +
  ggtitle("Combinations") +
  theme(axis.text.x = element_text(angle = 90))

library("ggpubr")
figure <- ggarrange(comb.3.p, ncell.p, up.p, down.p +
                      font("x.text", size = 12),
                    ncol = 2, nrow = 2)

png(paste0("./1_DiffPeak.outs/plot_all_treatandind_filt_unique_", reso, "_", FDR, "_", pos_foldchange, ".png"), width=2000, height=2000, pointsize=14, res=175)
print(figure)
dev.off()


##############################################################
###
### 1. MA plots
##############################################################
figfn <- paste0("./1_DiffPeak.outs/Figure1.1_MA_clusters_treat&ind_filt_", reso, "_", FDR, "_", pos_foldchange, ".png")
png(figfn, width=3500, height=4500, pointsize=12, res=225)
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




## ### 2, qq plots
## figfn <- paste0("./1_DiffPeak.outs/Figure1.2_qq_clusters_treat&ind_filt_", reso, "_", FDR, "_", pos_foldchange, ".png")
## png(figfn, width=3500, height=4500, pointsize=16, res=225)
## par(mar=c(4,4,2,2),mgp=c(2,1,0))
## x <- matrix(1:63, 9, 7, byrow=T)
## layout(x)
## #fn <- "./1_DiffPeak.outs/2.0_DESeq.results_clusters_treat&ind_filtered.rds"
## #res <- read_rds(fn)%>%drop_na(p.value)
## res.1 <- res%>%drop_na(p.value)
## # MCls <- c("CD4Naive", "TCM", "NKcell", "Bcell", "CD8Naive", "TEM", "Monocyte", "CTL", "dnT")
## MCls <- unique(cvt.2$MCls)
## for (oneMCl in MCls){
##   ##1
##   res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="caffeine")
##   print( qq(res.2$p.value, main="caffeine", cex.main=1, cex.axis=0.8, cex.lab=1))
##   ##2
##   res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="nicotine")
##   print( qq(res.2$p.value, main="nicotine", cex.main=1, cex.axis=0.8, cex.lab=1))
##   ##3 
##   res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="vitA")
##   print( qq(res.2$p.value, main="vitA", cex.main=1, cex.axis=0.8, cex.lab=1))
##   ##4
##   res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="vitD")
##   print( qq(res.2$p.value, main="vitD", cex.main=1, cex.axis=0.8, cex.lab=1))
##   ##5
##   res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="vitE")
##   print( qq(res.2$p.value, main="vitE", cex.main=1, cex.axis=0.8, cex.lab=1))
##   ##6
##   res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="water")
##   print( qq(res.2$p.value, main="water", cex.main=1, cex.axis=0.8, cex.lab=1))
##   ##7
##   res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="zinc")
##   print( qq(res.2$p.value, main="zinc", cex.main=1, cex.axis=0.8, cex.lab=1))   
##   print( mtext(oneMCl, side=4, line=0.5, cex=1, col="blue"))  
## }
## dev.off()


## ###
## ### 3, canno plots


## #------------------plot_loop--------------------------------
## #nrow(res) #if qq plot has been called before this step, this already is with dropped NA values
## res.1 <- res %>% drop_na(p.value)
## max.y <- floor(max(-log10(res.1$p.value)))/2
## #max.y
## png(paste0("./1_DiffPeak.outs/Figure1.3_qq_clusters_together_treatandind_filt_", reso,"_", FDR, "_", pos_foldchange, ".png"), width=2000, height=5000, pointsize=16, res=200)
## par(mar=c(4,4,2,2),mgp=c(2,1,0))
## x <- matrix(1:10, 5, 2, byrow=T)
## layout(x)
## # MCls <- c("CD4Naive", "TCM", "NKcell", "Bcell", "CD8Naive", "TEM", "Monocyte", "CTL", "dnT")
## MCls <- unique(cvt.2$MCls)
## for (oneMCl in MCls){
##   ##1
##   res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="caffeine")%>% drop_na(p.value)
##   plot(-log10(1:length(res.2$p.value)/length(res.2$p.value)),
##        -log10(sort(res.2$p.value)),
##        main=oneMCl,
##        pch = 19,
##        ylim = c(0, max.y),
##        ylab="observed -log10(p)",
##        xlab="expected -log10(p)",
##        col="red")
##   ##2
##   res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="nicotine")%>% drop_na(p.value)
##   points(-log10(1:length(res.2$p.value)/length(res.2$p.value)),-log10(sort(res.2$p.value)),  col="tan", pch = 19)   
##   ##3 
##   res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="zinc")%>% drop_na(p.value)
##   points(-log10(1:length(res.2$p.value)/length(res.2$p.value)),-log10(sort(res.2$p.value)),  col="maroon3", pch = 19)
##   ##4
##   res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="vitA")%>% drop_na(p.value)
##   points(-log10(1:length(res.2$p.value)/length(res.2$p.value)),-log10(sort(res.2$p.value)),  col="tan4", pch = 19)
##   ##5
##   res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="vitD")%>% drop_na(p.value)
##   points(-log10(1:length(res.2$p.value)/length(res.2$p.value)),-log10(sort(res.2$p.value)),  col="seagreen4", pch = 19)
##   ##6
##   res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="water")%>% drop_na(p.value)
##   points(-log10(1:length(res.2$p.value)/length(res.2$p.value)),-log10(sort(res.2$p.value)),  col="grey", pch = 19)
##   ##7
##   res.2 <- res.1%>%dplyr::filter(MCls==oneMCl, contrast=="vitE")%>% drop_na(p.value)
##   points(-log10(1:length(res.2$p.value)/length(res.2$p.value)),-log10(sort(res.2$p.value)),  col="salmon3", pch = 19)
##   ##control
##   abline(0, 1, col = "black")
##   ##legend
##   legend("topleft",
##          #inset=c(-0.2,0),
##          cex = 1,
##          pch = 19,
##          c("caff","nic", "zinc", "vitA", "vitD", "vitE", "water"),
##          fill=c("red","tan", "maroon3", "tan4", "seagreen4", "salmon3", "grey"))
## }
## dev.off()












## #####################################################################################################################################################
## ### if enriched in DEGs ###
## #####################################################################################################################################################
## DefaultAssay(combined) <- "ATAC"
## #granges(combined[["ATAC"]])

## gene.ranges = genes(EnsDb.Hsapiens.v86)
## #head(gene.ranges)

## seqlevelsStyle(gene.ranges) <- 'UCSC'
## #head(gene.ranges)

## peakAnno <- ClosestFeature(combined,
##                            regions=granges(combined),
##                            annotation = gene.ranges,
##                            sep = c(':', '-')
## )

## # head(peakAnno)
## # colnames(peakAnno)

## peakAnno <- peakAnno%>%dplyr::select(gene_name,query_region)%>%mutate(query_region=gsub(":", "-", query_region))
## # head(peakAnno)








## # res.DE.ATAC object
## #---------------------------------------------
## res <- res.DE.ATAC 
## #nrow(res)
## #table(res$MCls)
## #res <- res%>%mutate(MCls=gsub("NK", "NKcell", MCls))

## res <- res%>%dplyr::rename("peak_region"="gene")%>%drop_na(p.value)
## #res

## res <- res%>%
##   left_join(peakAnno,by=c("peak_region"="query_region"))%>%
##   mutate(comb=paste(MCls, contrast, sep="_"))
## #res
## #head(res$gene_name)


## # res.DE object -> DEGS
## #---------------------------------------------
## resDE <- res.DE
## #resDE
## #resDE <- resDE %>% filter(qval<0.1,abs(beta)>0.5)

## resDE <- resDE %>% dplyr::filter(p.adjusted<0.1,abs(estimate)>0.25)
## #unique(resDE$gene)



## # ENSEMBL to symbol
## #-------------------------------------------
## library(clusterProfiler)
## library(org.Hs.eg.db)

## #x <- bitr(unique(resDE$gene),fromType="ENSEMBL",toType=c("SYMBOL"),OrgDb=org.Hs.eg.db)
## #x <- bitr(unique(resDE$gene),fromType="GENENAME",toType=c("SYMBOL"),OrgDb=org.Hs.eg.db)
## #x <- bitr(unique(resDE$gene),fromType="GENENAME",toType=c("SYMBOL"),EnsDb=EnsDb.Hsapiens.v86)
## #head(x)

## #showMethods(keys)
## #help("SYMBOL")

## #resDE <- resDE%>%
## #   inner_join(x,by=c("gene"="ENSEMBL"))%>%
## #   mutate(comb=paste(MCls, contrast, sep="_"))


## resDE <- resDE%>%mutate(comb=paste(MCls, contrast, sep="_"))
## #head(resDE)


## # trial
## #---------------------------------------------
## table(res$gene_name)

## resDE.genes <- resDE$gene

## head(resDE.genes)

## new_df <- res %>% dplyr::filter(gene_name %in% resDE.genes)

## new_df


## # function
## ###-------------------------------------------
## comb <- sort(unique(resDE$comb))
## #comb

## res2 <- res%>%dplyr::filter(comb=="TEM_nicotine")
## #res2

## DEG <- resDE%>%dplyr::filter(comb=="TEM_nicotine")%>%dplyr::pull(gene)
## #DEG

## res2 <- res2%>%mutate(is_DEG=ifelse(gene_name%in%DEG,1,0))
## #table(res2$is_DEG)

## di <- res2%>%dplyr::filter(is_DEG==1)%>%arrange(p.value)
## #head(di$is_DEG)
## #ngene <- nrow(di)
## #ngene

## di <- di%>%mutate(observed=-log10(p.value), expected=-log10(ppoints(ngene)))
## #di


## ## dfNew <- map_dfr(comb, function(ii){
## ##   res2 <- res%>%filter(comb==ii)
## ##   DEG <- resDE%>%filter(comb==ii)%>%dplyr::pull(SYMBOL)
## ##   res2 <- res2%>%mutate(is_DEG=ifelse(gene_name%in%DEG,1,0))
## ##   ###
## ##   dx <- map_dfr(c(0,1),function(i){
## ##      di <- res2%>%filter(is_DEG==i)%>%arrange(p.value)
## ##      ngene <- nrow(di)
## ##      di <- di%>%mutate(observed=-log10(p.value), expected=-log10(ppoints(ngene)))
## ##      di
## ##   })
## ##   dx
## ## })

## dfNew <- map_dfr(comb, function(ii){
##   res2 <- res%>%dplyr::filter(comb==ii)
##   DEG <- resDE%>%dplyr::filter(comb==ii)%>%dplyr::pull(gene)
##   res2 <- res2%>%mutate(is_DEG=ifelse(gene_name%in%DEG,1,0))
##   ###
##   dx <- map_dfr(c(0,1),function(i){
##     di <- res2%>%dplyr::filter(is_DEG==i)%>%arrange(p.value)
##     ngene <- nrow(di)
##     di <- di%>%mutate(observed=-log10(p.value), expected=-log10(ppoints(ngene)))
##     di
##   })
##   dx
## })

## head(dfNew)

## table(dfNew$contrast)



## ### plot
## ###-----------------------------------------------
## #lab1 <- c(unique(dfNew$contrast))
## #labl
## lab1 <- c("caffeine"="caffeine", "nicotine"="nicotine", "vitA"="vitA", "vitD"="vitD", "vitE"="vitE", "water"="water", "zinc"="zinc")
## p1 <- ggplot(dfNew, aes(x=expected,y=observed, colour=factor(is_DEG)))+
##   ggrastr::rasterise(geom_point(size=0.3),dpi=300)+
##   geom_abline(colour="red")+
##   scale_colour_manual(values=c("0"="grey40","1"="green"),
##                       labels=c("0"="Not DEG", "1"="DEG"),
##                       guide=guide_legend(override.aes=list(size=1)))+
##   facet_grid(MCls~contrast, scales="free", labeller=labeller(contrast=lab1))+
##   xlab(bquote("Expected"~ -log[10]~"("~italic(plain(P))~")"))+
##   ylab(bquote("Observed"~-log[10]~"("~italic(plain(P))~")"))+
##   theme_bw()+
##   theme(legend.title=element_blank(),strip.text=element_text(size=12))

## figfn <- "./1_DiffPeak.outs/Figure1.3_qq.DEG.png"
## png(figfn, width=750, height=750, res=120)
## print(p1)
## dev.off()



#####################################################################################################################################################
### if enriched in DEGs ###
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
reso = 0.1

res.DE <- readRDS("./1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1.rds")

#res.DE

res.DE <- res.DE%>%mutate(comb=paste(MCls, contrast, sep="_"))
#head(resDE)

#res.DE.2 <- res.DE %>% dplyr::filter(p.adjusted<0.1,abs(estimate)>0.25)
#unique(resDE$gene)





# res.DE.ATAC object
#---------------------------------------------
res.DE.ATAC <- readRDS("./1_DiffPeak.outs/2.0_DESeq.results_clusters_treatandind_filtered_0.1.rds")

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
                                                  labels=c("0"="Not DEP", "1"="DEP"),
                                                  guide=guide_legend(override.aes=list(size=1)))+
        facet_grid(MCls~contrast, scales="free", labeller=labeller(contrast=lab1))+
        xlab(bquote("Expected"~ -log[10]~"("~italic(plain(P))~")"))+
        ylab(bquote("Observed"~-log[10]~"("~italic(plain(P))~")"))+
        theme_bw()+
        theme(legend.title=element_blank(),strip.text=element_text(size=12))
figfn <- paste0("./1_DiffPeak.outs/Figure1.3_qq.DEP_", reso,"_0.1_0.25.png")
png(figfn, width=750, height=750, res=120)
print(p1)
dev.off()

p2 <- ggplot(dfNew.2, aes(x=expected,y=observed, colour=factor(is_DEP)))+
      ggrastr::rasterise(geom_point(size=0.3),dpi=300)+
        geom_abline(colour="red")+
        scale_colour_manual(values=c("0"="grey40","1"="green"),
                                                  labels=c("0"="Not DEP", "1"="DEP"),
                                                  guide=guide_legend(override.aes=list(size=1)))+
        facet_grid(MCls~contrast, scales="free", labeller=labeller(contrast=lab1))+
        xlab(bquote("Expected"~ -log[10]~"("~italic(plain(P))~")"))+
        ylab(bquote("Observed"~-log[10]~"("~italic(plain(P))~")"))+
        theme_bw()+
        theme(legend.title=element_blank(),strip.text=element_text(size=12))
figfn <- paste0("./1_DiffPeak.outs/Figure1.3_qq.DEP_", reso, "_0.05_0.png")
png(figfn, width=750, height=750, res=120)
print(p2)
dev.off()


nrow(res.DE %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=0.1))
nrow(res.DE %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=0.05))
nrow(res.DE.ATAC %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=0.1))
nrow(res.DE.ATAC %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=0.05))

nrow(res.DE %>% drop_na(p.adjusted))
nrow(res.DE.ATAC %>% drop_na(p.adjusted))


combined[["ATAC"]]
combined[["RNA"]]

