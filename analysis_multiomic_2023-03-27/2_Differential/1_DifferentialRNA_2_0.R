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
library(cowplot,lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(grid)
library(gridExtra)
library(ggExtra)
library(RColorBrewer)
library(ggsci)
library(RColorBrewer)
library(viridis)
theme_set(theme_grey())

outdir <- "./1_DiffRNA_2.outs/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)

###
### Differential peaks analysis for raw features, that were not re-call peaks grouped for each cell-types
### Last motified by Julong Wei



####################################
### 1. Generate pseudo-bulk data ###
####################################
combined <- read_rds("../1_processing/3_Clustering.outs/1_seurat.cluster.combined.rds")
colnames(combined@meta.data)


combined <- subset(x = combined, subset = percent.mt.RNA < 20)
combined
max(combined@meta.data$percent.mt.RNA)


## meta <- atac@meta.data%>%
##    mutate(treat=gsub(".*-ATAC-|_.*", "", NEW_BARCODE)) 

meta <- combined@meta.data
head(meta)

nrow(meta)
colnames(meta)
table(meta$treat)

number_of_reads <- aggregate(nCount_RNA~treats,meta,sum)
number_of_reads

number_of_ATAC_reads <- aggregate(nCount_ATAC~treats,meta,sum)
number_of_ATAC_reads


meta$treat <- meta$treats
colnames(meta)

# cell <- c("0"="Tcell", "1"="Tcell", "2"="NKcell", "3"="Bcell",
#           "4"="Tcell", "5"="Tcell", "6"="Monocyte", "7"="Tcell", "8"="Tcell",
#           "9"="Tcell", "10"="Tcell", "11"="DC", "12"="DC")

# cell <- c("0"="CD4Naive", "1"="TCM", "2"="NKcell", "3"="Bcell",
#           "4"="CD8Naive", "5"="TEM", "6"="Monocyte", "7"="CTL", "8"="dnT",
#           "9"="Treg", "10"="Platelet", "11"="DC", "12"="HSPC")


cell <- c("0"="CD4Naive", "1"="TCM", "2"="NKcell", "3"="TEM",
          "4"="Bcell", "5"="CD8Naive", "6"="Monocyte", "7"="dnT", "8"="MAIT",
          "9"="Platelet", "10"="DC")


meta$MCls <- cell[as.character(meta$wsnn_res.0.1)]

colnames(meta)

meta <- meta%>%
   mutate(bti=paste(MCls, treat, SNG.BEST.GUESS, sep="_"))%>%
   dplyr::select(NEW_BARCODE, bti)
colnames(meta)

dd <- meta%>%group_by(bti)%>%summarise(ncell=n(),.groups="drop")
dd

write_rds(dd, "./1_DiffRNA_2.outs/0_ncell_clusters.rds")

dd <- read_rds("./1_DiffRNA_2.outs/0_ncell_clusters.rds")

dd


#----------------------------------------------------------------------------
#x <- atac@assays$ATAC@ranges
#x
#xrange <- ranges(x)
#xrange

#combined

count <- combined@assays$RNA@counts
head(count)

#count.1.col <- colnames(count)[colnames(count) %in% "scGxE_1-8*"]
#count.1.col
#count.1 <- select_counts(count,   )
                         
anno <- data.frame(rn=rownames(count), rnz=rowSums(count))
#str(anno)
head(anno)

#autosome <- paste("chr", as.character(1:22), sep="")
#autosome

annoSel <- anno%>%dplyr::filter(rnz>0)
head(annoSel)

Y <- count[annoSel$rn,]
#head(Y)


#------------------------------------------------------------------
##pseudo-bulk peak data
bti <- factor(meta$bti)
bti

X <- model.matrix(~0+bti)
head(X)


Y <- count
YtX <- Y %*% X

YtX <- as.matrix(YtX)
head(YtX)

colnames(YtX) <- gsub("^bti", "", colnames(YtX))
colnames(YtX)

###rnz>0,chr:1-22, 111,750*400
opfn <- "./1_DiffRNA_2.outs/1_YtX.comb_clusters.rds" 
write_rds(YtX, file=opfn)


#-------------------------------------------------------------------
###rnz>100, ncell>20, 111,746 * 331
dd <- read_rds("./1_DiffRNA_2.outs/0_ncell_clusters.rds")

dd <- read_rds("./1_DiffRNA_2.outs/0_ncell_clusters.rds")%>%filter(ncell>20)

table(dd$bti)

dd <- dd %>%filter(ncell>20)

dd

YtX <- read_rds("./1_DiffRNA_2.outs/1_YtX.comb_clusters.rds")
head(YtX)
ncol(YtX)
nrow(YtX)

anno <- data.frame(rn=rownames(YtX), rnz=rowSums(YtX))
head(anno)

annoSel <- anno%>%filter(rnz>100)
head(annoSel)

YtX_sel <- YtX[annoSel$rn, dd$bti]
head(YtX_sel)
ncol(YtX_sel)
nrow(YtX_sel)

opfn <- "./1_DiffRNA_2.outs/1_YtX.sel_clusters.rds"
write_rds(YtX_sel, file=opfn)





## genes >0 for each treatment
bti2 <- colnames(YtX)
length(bti2)

cvt0 <- str_split(bti2, "_", simplify=T)
head(cvt0)

cvt <- data.frame(bti=bti2, MCls=cvt0[,1], treat=cvt0[,2], sampleID=cvt0[,3])
head(cvt)
#table(cvt$MCls, cvt$sampleID) 
#table(cvt$treat, cvt$sampleID)
#table(cvt$treat)
#table(cvt$MCls)
#table(cvt$MCls, cvt$treat)
#comb <- table(cvt$MCls, cvt$treat)


treatments <- unique(cvt$treat)
treatments

for(i in treatments){
    cvt.1 <- cvt%>%dplyr::filter(treat==i)
    YtX.1 <- YtX[,cvt.1$bti]
    anno <- data.frame(rn=rownames(YtX.1), rnz=rowSums(YtX.1))
    ##head(anno)
    annoSel <- anno%>%filter(rnz>9)
    ##head(annoSel)
    print(i)
    print(nrow(annoSel))
    }











#####################################################
### 2. Differential analysis for CRP, High vs Low ###
#####################################################

### function, adjusted p value
myqval <- function(pval){
  qval <- pval
  ii0 <- !is.na(pval)
  qval[ii0] <- qvalue(pval[ii0],pi0=1)$qvalues
  qval
}


#### Read data
YtX <- read_rds("./1_DiffRNA_2.outs/1_YtX.sel_clusters.rds")

ncol(YtX)
nrow(YtX)

str(YtX)

colnames(YtX)

bti2 <- colnames(YtX)

length(bti2)

cvt0 <- str_split(bti2, "_", simplify=T)

head(cvt0)

cvt <- data.frame(bti=bti2, MCls=cvt0[,1], treat=cvt0[,2], sampleID=cvt0[,3])
#table(cvt$MCls, cvt$sampleID) 
#table(cvt$treat, cvt$sampleID)
#table(cvt$treat)
#table(cvt$MCls)
table(cvt$MCls, cvt$treat)

comb <- table(cvt$MCls, cvt$treat)

write.csv(comb, paste0("./1_DiffRNA_2.outs/combinations.csv"), quote = FALSE, row.names = TRUE)

comb <- read.csv("./1_DiffRNA_2.outs/combinations.csv")
comb

dd2 <- dd %>% dplyr::filter(bti %in% bti2)
sum(dd2$ncell)




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

cvt.filtered

cvt.2 <- cvt.filtered %>% dplyr::filter(count==1)
#cvt.2
#colnames(cvt.2)
#table(cvt.2$treat)
#table(cvt.2$MCls)

comb.2 <- table(cvt.2$MCls, cvt.2$treat)
comb.2

write.csv(comb.2, paste0("./1_DiffRNA_2.outs/combinations_2.csv"), quote = FALSE, row.names = TRUE)

comb.2 <- read.csv("./1_DiffRNA_2.outs/combinations_2.csv")
comb.2

bti3 <- cvt.2$bti 
#bti3

ncol(YtX)

YtX <- YtX[,bti3]
ncol(YtX)

dd2 <- dd %>% dplyr::filter(bti %in% bti3)
dd2
sum(dd2$ncell)





#################################################################
### 2.1 call DESeq ###
#################################################################
contrast.list <- list("caffeine"=c("treat", "caffeine", "water"),
                      "nicotine"=c("treat", "nicotine", "water"),
                      "vitA"=c("treat", "vitA", "etOH"),
                      "vitD"=c("treat", "vitD", "etOH"),
                      "vitE"=c("treat", "vitE", "etOH"),
                      "zinc"=c("treat", "zinc", "water"),
                      "water"=c("treat", "water", "etOH"))

cat("2.1", "run DESeq batch separately", "\n")

#MCls <- c("Bcell", "Monocyte", "NK", "Tcell", "DC")
###

#MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")

#MCls <- c("0", "1", "2", "3", "4", "5", "6", "7", "8")

MCls <- c("CD4Naive", "TCM", "NKcell", "Bcell",
          "CD8Naive", "TEM", "Monocyte", "CTL", "dnT")


## res <- map_dfr(MCls, function(oneX){
##   time0 <- Sys.time()
## ###
##   cvt0 <- cvt%>%dplyr::filter(MCls==oneX)
##   YtX0 <- YtX[,cvt0$bti]
## ##
##   dds <- DESeqDataSetFromMatrix(YtX0, cvt0, ~ sampleID + treat) #also controlling individual effect
##   dds <- DESeq(dds)
##   res <- contrast.list%>%map(~results(dds, contrast=.x))
##   res2 <- tibble(contrast=names(res),MCls=oneX, data=map(res,tidy))%>%unnest(data)  
## ##
##   time1 <- Sys.time()
##   elapsed <- difftime(time1, time0, units="mins")
##   cat(oneX, elapsed, "Done\n")
##   res2
## })

## head(res)

## table(res$contrast)

## table(res$MCls)

## opfn <- "./1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind.rds"
## write_rds(res, opfn)

## res <- readRDS("./1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind.rds")




res <- map_dfr(MCls, function(oneX){
  time0 <- Sys.time()
###
  cvt0 <- cvt.2%>%dplyr::filter(MCls==oneX)
  YtX0 <- YtX[,cvt0$bti]
##
  dds <- DESeqDataSetFromMatrix(YtX0, cvt0, ~ sampleID + treat) #also controlling individual effect
  dds <- DESeq(dds)
  res <- contrast.list%>%map(~results(dds, contrast=.x))
  res2 <- tibble(contrast=names(res),MCls=oneX, data=map(res,tidy))%>%unnest(data)  
##
  time1 <- Sys.time()
  elapsed <- difftime(time1, time0, units="mins")
  cat(oneX, elapsed, "Done\n")
  res2
})

head(res)

table(res$contrast)

table(res$MCls)

opfn <- "./1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered.rds"
write_rds(res, opfn)

res <- readRDS("./1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered.rds")








#colnames(combined@meta.data)
#table(combined@meta.data$seurat_clusters)

res

nrow(res)
nrow(res %>% drop_na(p.adjusted))
nrow(res %>% drop_na(p.value))

######################################################
# significant genes
######################################################
padj_cutoff <- 0.05

sig_res <- dplyr::filter(res, p.adjusted < padj_cutoff) %>% dplyr::arrange(p.adjusted)

#write.csv(sig_res, paste0("./1_DiffRNA_2.outs/sig_genes.csv"), quote = FALSE, row.names = FALSE)

#sig_res <- read.csv("./1_DiffRNA_2.outs/sig_genes.csv")

head(sig_res)

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


## #-----------------------------------------------------
## treats <- c(levels(factor(sig_res$contrast)))
## treats

## #MCls <- c("Bcell",  "Monocyte", "NKcell", "Tcell")
## MCls

## head(sig_res)

## table(sig_res$MCls, sig_res$contrast)

## #dnT_caffeine <- sig_res %>% dplyr::filter(MCls=="dnT", contrast=="caffeine") %>% dplyr::arrange(p.adjusted) %>% dplyr::pull(gene) %>% head(n=20)
## #dnT_caffeine

## #dnT_zinc <- sig_res %>% dplyr::filter(MCls=="dnT", contrast=="zinc") %>% dplyr::arrange(p.adjusted) %>% dplyr::pull(gene) %>% head(n=20)
## #dnT_zinc

## #-------------list----------------------------------
## sig_genes <- list()

## sig_genes

## for(i in length(MCls)){
##     new_list <- list()
##     sig_genes <- append(sig_genes[[1]], new_list)
##     }

## sig_genes

## length(sig_genes)


## #-------------loop----------------------------------
## for(MCl in MCls){
##     #print(MCl)
##     for (treat in treats){
##         #print(treat)
##         dfname <- paste0(MCl, "_", treat, sep ="")
##         print(dfname)
##         dfname <- sig_res %>% dplyr::filter(MCls==MCl) %>%  dplyr::filter(contrast==treat)  %>% dplyr::arrange(p.adjusted) %>% dplyr::pull(gene) %>% head(n=20)
##         print(dfname)
##     }}


## for(MCl in MCls){
##     #print(MCl)
##     for (treat in treats){
##         #print(treat)
##         dfname <- paste0(MCl, "_", treat, sep ="")
##         print(dfname)
##         dfname <- sig_res %>% dplyr::filter(MCls==MCl) %>%  dplyr::filter(contrast==treat)  %>% dplyr::arrange(p.adjusted) %>% dplyr::pull(gene) %>% head(n=20)
##         print(dfname)
##     }}

## Bcell_caffeine

## ## for(MCl in MCls){
## ##     MCl <- sig_res %>% dplyr::filter(MCls==MCl)  %>% dplyr::arrange(p.adjusted) %>% dplyr::pull(gene) %>% head(n=20)
## ## }
## ## MCl



#######################################################
# table of up and down regulated genes
#######################################################
#------------------------------------------------------
# Nicole script
#------------------------------------------------------
FDR <- 0.1
pos_foldchange <- 0.25
neg_foldchange <- -0.25

## min((res %>% drop_na(estimate))$estimate)

## nrow(res %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=0.1) %>% dplyr::filter(estimate <0.25, estimate >-0.25))

## nrow(res %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=0.1)  %>% dplyr::filter(estimate >=0.25))

## nrow(res %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=0.1) %>% dplyr::filter(estimate <=-0.25))

## nrow(res %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<=0.1) %>% dplyr::filter(estimate <=0.25))

#res$Significance <- NULL

# get DESeq results:
#res <- results(ddsFull,parallel=TRUE)

res[!is.na(res$p.adjusted) & res$p.adjusted<=FDR & res$estimate>=pos_foldchange, "Significance"] <- "Upregulated"

res[!is.na(res$p.adjusted) & res$p.adjusted<=FDR & res$estimate<=neg_foldchange, "Significance"] <- "Downregulated"

res[!is.na(res$p.adjusted) & res$p.adjusted>FDR,"Significance"] <- "Not Significant"

res[is.na(res$p.adjusted),"Significance"] <- "NA"

#res2 <- res %>% dplyr::filter(Significance %in% c("Upregulated", "Downregulated"))

res_up <- res %>% dplyr::filter(Significance %in% "Upregulated")

res_down <- res %>% dplyr::filter(Significance %in% "Downregulated")

up <- table(res_up$MCls, res_up$contrast)
up

write.csv(up, paste0("./1_DiffRNA_2.outs/up_treat&ind_filt.csv"), quote = FALSE, row.names = TRUE)

down <- table(res_down$MCls, res_down$contrast)
down

write.csv(down, paste0("./1_DiffRNA_2.outs/down_treat&ind_filt.csv"), quote = FALSE, row.names = TRUE)

## up_down <- table(res2$MCls, res2$contrast, res2$Significance)
## up_down
## write.csv(up_down, paste0("./1_DiffRNA_2.outs/up_down.csv"), quote = FALSE, row.names = FALSE)


up <- read.csv(file = './1_DiffRNA_2.outs/up_treat&ind_filt.csv')
#up

down <- read.csv(file = './1_DiffRNA_2.outs/down_treat&ind_filt.csv')
#down

## #### this is repeat from above, incase starting from here #############
## YtX <- read_rds("./1_DiffRNA_2.outs/1_YtX.sel_clusters.rds")

## colnames(YtX)
## length(YtX)

## bti2 <- colnames(YtX)
## bti2
## length(bti2)

## cvt0 <- str_split(bti2, "_", simplify=T)
## head(cvt0)

## cvt <- data.frame(bti=bti2, MCls=cvt0[,1], treat=cvt0[,2], sampleID=cvt0[,3])
## head(cvt)
## nrow(cvt)

## table(cvt$treat)
## table(cvt$MCls)

## dd2 <- dd %>% dplyr::filter(bti %in% bti2)
## head(dd2)
## nrow(dd2)
## sum(dd2$ncell)



#-----------number of cells------------------------------------------------
cvt.2

dd3 <- left_join(cvt.2, dd2)

head(dd3)

res

dd3 <- dd3 %>% dplyr::filter(dd3$MCls %in% res$MCls, !dd3$treat %in% "etOH")
head(dd3)
nrow(dd3)
sum(dd3$ncell)

dd4 <- dd3 %>% group_by(treat, MCls) %>% summarise(ncell_count = sum(ncell))

dd4

ncell_table <- data.frame(pivot_wider(dd4, names_from=treat, values_from=ncell_count))

rownames(ncell_table) <- ncell_table$MCls
ncell_table$MCls <- NULL

ncell_table

write.csv(ncell_table, paste0("./1_DiffRNA_2.outs/ncell_treat$ind_filt.csv"), quote = FALSE, row.names = TRUE)

ncell <- read.csv("./1_DiffRNA_2.outs/ncell_treat$ind_filt.csv")

#------------------------not using------------------------------------------
## res.df <- data.frame(res)
## head(res.df)

## table(res.df$MCls, res.df$contrast, res.df$Significance)


## # how many DEGs:
## cat("##10%FDR:\t",save.name,"\t",sum(res$padj<0.1,na.rm=TRUE),"\t",dim(cv)[1],"\n")

## # Save results.
## write.table(res,file=gzfile(paste("./all-out_data/", save.name,"_stats.txt.gz",sep="")),quote=F,sep="\t")

## # drop untested:
## res <- res[!is.na(res$padj),]

## # save DEGs:
## DEensg <- rownames(res[res$padj<0.1,])
## DEensgup <- rownames(res[res$padj<0.1 & res$log2FoldChange>0,])
## DEensgdown <- rownames(res[res$padj<0.1 & res$log2FoldChange<0,])
## write.table(DEensg,file=paste0("./DEGs/",save.name,"_DEGs_10FDR.txt"),quote=F,sep="\t", row.names=F,col.names=F)
## write.table(DEensgup,file=paste0("./DEGs/",save.name,"_up_DEGs_10FDR.txt"),quote=F,sep="\t", row.names=F,col.names=F)
## write.table(DEensgdown,file=paste0("./DEGs/",save.name,"_down_DEGs_10FDR.txt"),quote=F,sep="\t", row.names=F,col.names=F)

## # save all genes to use as background:
## bckgrnd <- rownames(res[!is.na(res$padj),])
## write.table(bckgrnd,file=paste0("./background_ensgs/background_ensgs_",save.name,".txt"),quote=F,sep="\t", row.names=F,col.names=F)


#########################################################
# matrix_plots
#########################################################
up
class(up)
dim(up)

rownames(up)<- up$X
up$X <- NULL
up
colnames(up)
rownames(up)

#up.df <- as.data.frame(up)

#library(plot.matrix)
## y<-matrix(up)
## y
## class(y)
## par(mar=c(5.1, 4.1, 4.1, 4.1)) # adapt margins
## plot(x)

library(reshape2)

up.matrix <-as.matrix(up) # create a numeric matrix object
#up.matrix
#class(up.matrix)
up.melt <- melt(up.matrix)
#up.melt
#colnames(up.melt)

#down
rownames(down)<- down$X
down$X <- NULL
down.matrix <- as.matrix(down)
down.melt <- melt(down.matrix)

#ncell
rownames(ncell)<- ncell$X
ncell$X <- NULL
ncell.matrix <- as.matrix(ncell)
ncell.melt <- melt(ncell.matrix)

#comb
comb <- comb.2
comb

#rownames(comb)<- comb$X
#comb$X <- NULL
comb.2 <- subset(comb, select = -c(etOH))
comb.2 <- comb.2 %>% dplyr::filter(rownames(comb.2) %in% rownames(up))
#comb.2 <- comb[-c("Platelet", "Treg")]
#comb.2
comb.2.matrix <- as.matrix(comb.2)
comb.2.melt <- melt(comb.2.matrix)
#comb.2.melt


#figfn <- "./1_DiffPeak.outs/Figure1.1_MA_clusters_treatandind.png"

#png("./1_DiffRNA_2.outs/plot_up_treatandind.png", width=500, height=500, pointsize=14, res=150)
up.p <- ggplot(up.melt, aes(x = Var2, y = Var1)) +
    geom_tile(aes(fill=value))+
    geom_text(aes(label = value), color = "white", size = 4) +
    ggtitle("Upregulated") +
    theme(axis.text.x = element_text(angle = 90))

#print(p)
#dev.off()

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

comb.2.p <- ggplot(comb.2.melt, aes(x = Var2, y = Var1)) +
      geom_tile(aes(fill=value))+
      geom_text(aes(label = value), color = "white", size = 4) +
      ggtitle("Combinations") +
      theme(axis.text.x = element_text(angle = 90))

library("ggpubr")

figure <- ggarrange(comb.2.p, ncell.p, up.p, down.p +
                    font("x.text", size = 12),
                    ncol = 2, nrow = 2)


png("./1_DiffRNA_2.outs/plot_all_treatandind_filt.png", width=2000, height=2000, pointsize=14, res=175)
print(figure)
dev.off()


##############################################################
###
### 1. MA plots
##############################################################
figfn <- "./1_DiffRNA_2.outs/Figure1.1_MA_clusters_treat&ind_filt.png"
#png(figfn, width=900, height=1000, pointsize=12, res=150)
png(figfn, width=3500, height=4500, pointsize=12, res=225)
par(mar=c(4,4,2,2),mgp=c(2,1,0))
x <- matrix(1:63, 9, 7, byrow=T)
layout(x)
fn <- "./1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered.rds"
res <- read_rds(fn)%>%
       mutate(color=ifelse((p.adjusted<0.1)&(!is.na(p.adjusted)), T, F))
#MCls <- c("Bcell",  "Monocyte", "NKcell", "Tcell")
#MCls <- c("0", "1", "2", "3", "4", "5", "6", "7", "8")
MCls <- c("CD4Naive", "TCM", "NKcell", "Bcell",
          "CD8Naive", "TEM", "Monocyte", "CTL", "dnT")
for (oneMCl in MCls){
##1
   res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="caffeine")%>%
           dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
   print(plotMA(res2[,1:3], colLine="NA", main="caffeine", cex.main=1, cex.axis=0.8, cex.lab=1))
###2    
   res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="nicotine")%>%
          dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
   print(plotMA(res2[,1:3], colLine="NA", main="nicotine", cex.main=1, cex.axis=0.8, cex.lab=1))
###3
   res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="vitA")%>%
           dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
   print(plotMA(res2[,1:3], colLine="NA", main="vitA", cex.main=1, cex.axis=0.8, cex.lab=1))
###4
   res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="vitD")%>%
           dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
   print(plotMA(res2[,1:3], colLine="NA", main="vitD", cex.main=1, cex.axis=0.8, cex.lab=1))
###5
    res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="vitE")%>%
           dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
   print(plotMA(res2[,1:3], colLine="NA", main="vitE", cex.main=1, cex.axis=0.8, cex.lab=1))
###6    
    res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="zinc")%>%
           dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
   print(plotMA(res2[,1:3], colLine="NA", main="zinc", cex.main=1, cex.axis=0.8, cex.lab=1))
###7
    res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="water")%>%
           dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
   print(plotMA(res2[,1:3], colLine="NA", main="water", cex.main=1, cex.axis=0.8, cex.lab=1))    
   print(mtext(oneMCl, side=4, line=0.5, cex=1, col="blue")) 
}
dev.off()



## figfn <- "./1_DiffRNA_2.outs/Figure1.1_MA_clusters_2.png"
## #png(figfn, width=900, height=1000, pointsize=12, res=150)
## png(figfn, width=900, height=1000, pointsize=12, res=300)
## par(mar=c(4,4,2,2),mgp=c(2,1,0))
## x <- matrix(1:63, 9, 7, byrow=T)
## layout(x)
## fn <- "./1_DiffRNA_2.outs/2.0_DESeq.results_clusters.rds"
## res <- read_rds(fn)%>%
##        mutate(color=ifelse((p.adjusted<0.1)&(!is.na(p.adjusted)), T, F))
## #MCls <- c("Bcell",  "Monocyte", "NKcell", "Tcell")
## MCls <- c("4", "5", "6", "7", "8")
## for (oneMCl in MCls){
## ##1
##    res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="caffeine")%>%
##            dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
##    print(plotMA(res2[,1:3], colLine="NA", main="caffeine", cex.main=1, cex.axis=0.8, cex.lab=1))
## ###2    
##    res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="nicotine")%>%
##           dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
##    print(plotMA(res2[,1:3], colLine="NA", main="nicotine", cex.main=1, cex.axis=0.8, cex.lab=1))
## ###3
##    res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="vitA")%>%
##            dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
##    print(plotMA(res2[,1:3], colLine="NA", main="vitA", cex.main=1, cex.axis=0.8, cex.lab=1))
## ###4
##    res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="vitD")%>%
##            dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
##    print(plotMA(res2[,1:3], colLine="NA", main="vitD", cex.main=1, cex.axis=0.8, cex.lab=1))
## ###5
##     res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="vitE")%>%
##            dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
##    print(plotMA(res2[,1:3], colLine="NA", main="vitE", cex.main=1, cex.axis=0.8, cex.lab=1))
## ###6    
##     res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="zinc")%>%
##            dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
##    print(plotMA(res2[,1:3], colLine="NA", main="zinc", cex.main=1, cex.axis=0.8, cex.lab=1))
## ###7
##     res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="water")%>%
##            dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
##    print(plotMA(res2[,1:3], colLine="NA", main="water", cex.main=1, cex.axis=0.8, cex.lab=1))    
##    print(mtext(oneMCl, side=4, line=0.5, cex=1, col="blue")) 
## }
## dev.off()





##Example 
## fn <- "./1_DP.outs/2.0_DESeq.results.rds"
## res <- read_rds(fn)%>%drop_na(p.value)%>%
##        mutate(color=ifelse(p.adjusted<0.1, T, F))
## oneMCl <- "NK"
## res2 <- res%>%dplyr::filter(MCls==oneMCl)%>%dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
## figfn <- "./1_DP.outs/Figure1.1_MA2.png"
## png(figfn, width=800, height=1000, pointsize=12, res=150)
## print(plotMA(res2[,1:3], colLine="NA", main=oneMCl, cex.main=1, cex.axis=0.8, cex.lab=1))
## dev.off()


### 2, qq plots
figfn <- "./1_DiffRNA_2.outs/Figure1.2_qq_clusters_treat&ind_filt.png"
#png(figfn, width=900, height=1000, pointsize=12, res=150)
png(figfn, width=3500, height=4500, pointsize=16, res=225)
par(mar=c(4,4,2,2),mgp=c(2,1,0))
x <- matrix(1:63, 9, 7, byrow=T)
layout(x)
fn <- "./1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered.rds"
res <- read_rds(fn)%>%drop_na(p.value)
#MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
#MCls <- c("0", "1", "2", "3", "4", "5", "6", "7", "8")
MCls <- c("CD4Naive", "TCM", "NKcell", "Bcell",
          "CD8Naive", "TEM", "Monocyte", "CTL", "dnT")
for (oneMCl in MCls){
   ##1
   res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="caffeine")
   print( qq(res2$p.value, main="caffeine", cex.main=1, cex.axis=0.8, cex.lab=1))
   ##2
   res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="nicotine")
   print( qq(res2$p.value, main="nicotine", cex.main=1, cex.axis=0.8, cex.lab=1))
   ##3 
   res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="vitA")
   print( qq(res2$p.value, main="vitA", cex.main=1, cex.axis=0.8, cex.lab=1))
   ##4
   res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="vitD")
   print( qq(res2$p.value, main="vitD", cex.main=1, cex.axis=0.8, cex.lab=1))
   ##5
   res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="vitE")
   print( qq(res2$p.value, main="vitE", cex.main=1, cex.axis=0.8, cex.lab=1))
   ##6
   res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="zinc")
   print( qq(res2$p.value, main="zinc", cex.main=1, cex.axis=0.8, cex.lab=1))
   ##7
   res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="water")
   print( qq(res2$p.value, main="water", cex.main=1, cex.axis=0.8, cex.lab=1))
   print( mtext(oneMCl, side=4, line=0.5, cex=1, col="blue"))  
}
dev.off()


###
### 3, canno plots


res

table(res$contrast)

res.B.caff <- res%>%dplyr::filter(MCls=="Bcell", contrast=="caffeine")
res.B.zinc <- res%>%dplyr::filter(MCls=="Bcell", contrast=="zinc")
res.B <- res%>%dplyr::filter(MCls=="Bcell")

res.B


#---------------stat_qq-------------------------------------
figfn <- "./1_DiffRNA_2.outs/Figure1.2_qq_clusters_together.png"
#png(figfn, width=900, height=1000, pointsize=12, res=150)
png(figfn, width=960, height=960, pointsize=16, res=200)
ggplot() +
      stat_qq(aes(sample = res.B.caff$p.value), colour = "green") +
      stat_qq(aes(sample = res.B.zinc$p.value), colour = "red") +
      geom_abline(aes(slope = 1, intercept = 0), linetype = 2)
dev.off()

#---------------q_plot-------------------------------------
png("./1_DiffRNA_2.outs/Figure1.2_qq_clusters_together.png", width=960, height=960, pointsize=16, res=200)
p<-qplot(sample = p.value, data = res.B, color=contrast)+theme_bw()
print(p)
dev.off()

png("./1_DiffRNA_2.outs/Figure1.2_qq_clusters_together.png", width=960, height=960, pointsize=16, res=200)
p<-qplot(sample = -log10(p.value), data = res.B, color=contrast)+theme_bw()
print(p)
dev.off()

#---------------qqnorm & qqline-------------------------------------
png("./1_DiffRNA_2.outs/Figure1.2_qq_clusters_together.png", width=960, height=960, pointsize=16, res=200)
qqnorm(-log10(res.B$p.value))
dev.off()

#---------------plot-------------------------------------
#https://stats.stackexchange.com/questions/92124/how-to-interpret-a-qq-plot-of-p-values
observedPValues_B.caff <- na.omit(res.B.caff$p.value)
observedPValues_B.zinc <- na.omit(res.B.zinc$p.value)
res
png("./1_DiffRNA_2.outs/Figure1.2_qq_clusters_together.png", width=960, height=960, pointsize=16, res=200)
plot(-log10(1:length(observedPValues_B.caff)/length(observedPValues_B.caff)),
     -log10(sort(observedPValues_B.caff)),
     main="Bcell",
     ylab="observed -log10(p)",
     xlab="expected -log10(p)",
     col="red")
points(-log10(1:length(observedPValues_B.zinc)/length(observedPValues_B.zinc)),
       -log10(sort(observedPValues_B.zinc)),
       col="blue")
abline(0, 1, col = "black")
dev.off()


#------------------plot_loop--------------------------------
colors()

nrow(res) #if qq plot has been called before this step, this already is with dropped NA values

res.drop.na <- res %>% drop_na(p.value)

nrow(res.drop.na)

max.y <- floor(max(-log10(res.drop.na$p.value)))/2

max.y

png("./1_DiffRNA_2.outs/Figure1.2_qq_clusters_together_treatandind.png", width=2000, height=5000, pointsize=16, res=200)
par(mar=c(4,4,2,2),mgp=c(2,1,0))
x <- matrix(1:10, 5, 2, byrow=T)
layout(x)
MCls <- c("CD4Naive", "TCM", "NKcell", "Bcell",
          "CD8Naive", "TEM", "Monocyte", "CTL", "dnT")
#MCls <- c("Bcell")
for (oneMCl in MCls){
    ##1
    res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="caffeine")%>% drop_na(p.value)
    plot(-log10(1:length(res2$p.value)/length(res2$p.value)),
         -log10(sort(res2$p.value)),
         main=oneMCl,
         pch = 19,
         ylim = c(0, max.y),
         ylab="observed -log10(p)",
         xlab="expected -log10(p)",
         col="red")
    ##2
    res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="nicotine")%>% drop_na(p.value)
    points(-log10(1:length(res2$p.value)/length(res2$p.value)),-log10(sort(res2$p.value)),  col="tan", pch = 19)   
    ##3 
    res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="zinc")%>% drop_na(p.value)
    points(-log10(1:length(res2$p.value)/length(res2$p.value)),-log10(sort(res2$p.value)),  col="maroon3", pch = 19)
    ##4
    res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="vitA")%>% drop_na(p.value)
    points(-log10(1:length(res2$p.value)/length(res2$p.value)),-log10(sort(res2$p.value)),  col="tan4", pch = 19)
    ##5
    res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="vitD")%>% drop_na(p.value)
    points(-log10(1:length(res2$p.value)/length(res2$p.value)),-log10(sort(res2$p.value)),  col="seagreen4", pch = 19)
    ##6
    res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="vitE")%>% drop_na(p.value)
    points(-log10(1:length(res2$p.value)/length(res2$p.value)),-log10(sort(res2$p.value)),  col="salmon3", pch = 19)
    ##7
    res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="water")%>% drop_na(p.value)
    points(-log10(1:length(res2$p.value)/length(res2$p.value)),-log10(sort(res2$p.value)),  col="orange", pch = 19)
    ##control
    abline(0, 1, col = "black")
    ##legend
    legend("topleft",
           #inset=c(-0.2,0),
           cex = 1,
           pch = 19,
           c("caff","nic", "zinc", "vitA", "vitD", "vitE", "water"),
           fill=c("red","tan", "maroon3", "tan4", "seagreen4", "salmon3", "orange"))
}
dev.off()



#----------plot function-individual plots----------------------------------------
MCls <- c("CD4Naive", "TCM", "NKcell", "Bcell",
          "CD8Naive", "TEM", "Monocyte", "CTL", "dnT")
#MCls <- c("Bcell")

for (oneMCl in MCls){
    png(paste("./1_DiffRNA_2.outs/Figure1.2_qq_clusters_", oneMCl, "_together_y60.png", sep=""), pointsize=12, width=960, height=960, res=200)
    ##1
    res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="caffeine")%>% drop_na(p.value)
    plot(-log10(1:length(res2$p.value)/length(res2$p.value)),
         -log10(sort(res2$p.value)),
         main=oneMCl,
         pch = 19,
         ylim = c(0, 60),
         ylab="observed -log10(p)",
         xlab="expected -log10(p)",
         col="red")
    ##2
    res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="nicotine")%>% drop_na(p.value)
    points(-log10(1:length(res2$p.value)/length(res2$p.value)),-log10(sort(res2$p.value)),  col="grey", pch = 19, ylim = c(0, 60))   
    ##3 
    res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="zinc")%>% drop_na(p.value)
    points(-log10(1:length(res2$p.value)/length(res2$p.value)),-log10(sort(res2$p.value)),  col="green", pch = 19, ylim = c(0, 60))
    ##4
    res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="vitA")%>% drop_na(p.value)
    points(-log10(1:length(res2$p.value)/length(res2$p.value)),-log10(sort(res2$p.value)),  col="blue", pch = 19, ylim = c(0, 60))
    ##5
    res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="vitD")%>% drop_na(p.value)
    points(-log10(1:length(res2$p.value)/length(res2$p.value)),-log10(sort(res2$p.value)),  col="pink", pch = 19, ylim = c(0, 60))
    ##6
    res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="vitE")%>% drop_na(p.value)
    points(-log10(1:length(res2$p.value)/length(res2$p.value)),-log10(sort(res2$p.value)),  col="yellow", pch = 19, ylim = c(0, 60))
    ##7
    res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="water")%>% drop_na(p.value)
    points(-log10(1:length(res2$p.value)/length(res2$p.value)),-log10(sort(res2$p.value)),  col="orange", pch = 19, ylim = c(0, 60))
    ##control
    abline(0, 1, col = "black")
    ##legend
    legend("topleft",
           #inset=c(-0.2,0)
           cex=1,
           pch=19,
           c("caff","nic", "zinc", "vitA", "vitD", "vitE", "water"),
           fill=c("red","grey", "green", "blue", "pink", "yellow", "orange"))
    dev.off()
}



#-----------trial_data------------------------
res3 <- (res%>%dplyr::filter(MCls=="Bcell", contrast=="zinc")%>% drop_na(p.value))$p.value

res4 <- -log10(res3)

max(res4)

res.filter <- res%>% drop_na(p.value)

table(res.filter$MCls, res.filter$contrast)

nrow(res%>%dplyr::filter(MCls=="CD4Naive", contrast=="caffeine")%>% drop_na(p.value))

nrow(res%>%dplyr::filter(MCls=="CD4Naive", contrast=="zinc")%>% drop_na(p.value))



#-----------ggplot2------------------------
res2 <- res%>%dplyr::filter(MCls=="Bcell", contrast=="caffeine")%>% drop_na(p.value)
p1 <- ggplot() +
      geom_point(data=res2, aes(-log10(1:length(p.value)/length(p.value)),
                                                             -log10(sort(p.value))),
                              color="red") +
      ggtitle("Bcell") +
      ylab("observed -log10(p)") +
      xlab("expected -log10(p)")

res2 <- res%>%dplyr::filter(MCls=="Bcell", contrast=="zinc")%>% drop_na(p.value)
p1 <- p1 + geom_point(data=res2, aes(-log10(1:length(p.value)/length(p.value)),
                                                               -log10(sort(p.value))),
                                            color="blue")
print(p1)




#-----------ggplot2-groupby------------------------
res.sort <- res %>% dplyr::arrange(MCls, contrast, p.value)
res.sort



MCls = unique(res.sort$MCls)
MCls

treats = unique(res.sort$contrast)
treats

res.split <- split(res, f = list(res$contrast,res$MCls))

MCls.treats<-names(res.split)
MCls.treats[1]

cvt0 <- str_split(MCls.treats[1], "\\.", simplify=T)
cvt0[1]

# MCls.treats<-nameslist()
# for (i in length(MCls)){
#   for (j in length(treats)){
#     new_element <- paste(MCls[i], treats[j])
#     MCls.treats[[length(MCls.treats) + 1]] <- new_element
#   }
# }
# MCls.treats



res.exp.p <- map_dfr(MCls.treats, function(oneX){
      time0 <- Sys.time()
        ###
        res.sort.2 <- res.sort%>%dplyr::filter(MCls==str_split(oneX, "\\.", simplify=T)[2], contrast==str_split(oneX, "\\.", simplify=T)[1])
        ##
        res.sort.2$exp.p.value <- 1:length(res.sort.2$p.value)/length(res.sort.2$p.value)
        ##
        time1 <- Sys.time()
        elapsed <- difftime(time1, time0, units="mins")
        cat(oneX, elapsed, "Done\n")
        res.sort.2
      })

res.exp.p

colnames(res.exp.p)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
res.exp.p.B <- res.exp.p %>% dplyr::filter(MCls=="Bcell")

p1 <- ggplot() +
      geom_point(data=res.exp.p.B, aes(-log10(exp.p.value),
                                                                   -log10(p.value), color=contrast)) +
      ggtitle("Bcell")+
      ylab("observed -log10(p)") +
      xlab("expected -log10(p)")
print(p1)



#------------------ggplot_groupby_loop--------------------------------
library(ggpubr)

#png("./1_DiffRNA_2.outs/Figure1.2_qq_clusters_ggplot_all_together.png", width=960, height=2250, pointsize=16, res=200)

## par(mar=c(4,4,2,2),mgp=c(2,1,0))
## x <- matrix(1:10, 5, 2, byrow=T)
## layout(x)
## x

MCls.types <- c("CD4Naive", "TCM", "NKcell", "Bcell",
                   "CD8Naive", "TEM", "Monocyte", "CTL", "dnT")
# MCls <- c("Bcell")

plot_list = list()

for (i in length(MCls.types)){
      ##1
      res.loop <- res.exp.p %>% dplyr::filter(MCls==MCls.types[i])
      p1 <- ggplot() +
              geom_point(data=res.loop, aes(-log10(exp.p.value),
                                            -log10(p.value), color=contrast)) +
              ggtitle(MCls.types[i])+
              ylab("observed -log10(p)") +
              xlab("expected -log10(p)")
    plot_list[[i]] = p1
    }

    

plot_list

png("./1_DiffRNA_2.outs/Figure1.2_qq_clusters_ggplot_all_together.png", width=960, height=2250, pointsize=16, res=200)
p <- ggarrange(plot_list[[1]], plot_list[[2]],
          plot_list[[3]], plot_list[[4]],
          plot_list[[5]], plot_list[[6]],
          plot_list[[7]], plot_list[[8]],
          plot_list[[9]],
          ncol = 2, nrow = 5)
print(p)
dev.off()






#-----------ggplot2 loop------------------------
# for pch -> http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
#MCls <- c("Bcell")
#colors <- c("caff"="red","nic"="grey", "zinc"="green", "vitA"="blue", "vitD"="pink", "vitE"="yellow", "water"="orange")
MCls <- c("Bcell")
# par(mar=c(4,4,2,2),mgp=c(2,1,0))
# x <- matrix(1:10, 5, 2, byrow=T)
# layout(x)
for (oneMCl in MCls){
      ##1
      res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="caffeine")%>% drop_na(p.value)
      p1 <- ggplot() +
              geom_point(data=res2, aes(-log10(1:length(p.value)/length(p.value)), -log10(sort(p.value))), colour="red") +
              ggtitle(oneMCl) +
              ylab("observed -log10(p)") +
              xlab("expected -log10(p)")+
              #scale_colour_manual(name = "Treatments", values = c("red","grey", "green", "blue", "pink", "yellow", "orange"), labels = c("caff","nic", "zinc", "vitA", "vitD", "vitE", "water"))
            ##2
            res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="nicotine")%>% drop_na(p.value)
      p1 <- p1 + geom_point(data=res2, aes(-log10(1:length(p.value)/length(p.value)),-log10(sort(p.value))), colour="grey")
      ##3
      res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="zinc")%>% drop_na(p.value)
      p1 <- p1 + geom_point(data=res2, aes(-log10(1:length(p.value)/length(p.value)),-log10(sort(p.value))), colour="green")
      ##4
      res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="vitA")%>% drop_na(p.value)
      p1 <- p1 + geom_point(data=res2, aes(-log10(1:length(p.value)/length(p.value)),-log10(sort(p.value))), colour="blue")
      ##5
      res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="vitD")%>% drop_na(p.value)
      p1 <- p1 + geom_point(data=res2, aes(-log10(1:length(p.value)/length(p.value)),-log10(sort(p.value))), colour="pink")
      ##6
      res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="vitE")%>% drop_na(p.value)
      p1 <- p1 + geom_point(data=res2, aes(-log10(1:length(p.value)/length(p.value)),-log10(sort(p.value))), colour="yellow")
      ##7
      res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="water")%>% drop_na(p.value)
      p1 <- p1 + geom_point(data=res2, aes(-log10(1:length(p.value)/length(p.value)),-log10(sort(p.value))), colour="orange")
      ##control
      p1 <- p1 + geom_abline(intercept=0, slope=1, colour = "black")
      ##legend
      #labs(colour = "Treatments")
      #p1 +
      print(p1)
    }


#"red"="red","grey"="grey", "green"="green", "blue"="blue", "pink"="pink", "yellow"="yellow", "orange"="orange"






## ###################################
## ### Table of Differential peaks ###
## ###################################
## fn <- "./1_DiffRNA_2.outs/2.0_DESeq.results.rds"
## res <- read_rds(fn)%>%drop_na(p.value)

## res%>%filter(p.adjusted<0.1,abs(estimate)>0.5)%>%
##    group_by(MCls, contrast)%>%
##    summarise(ngene=n(),.groups="drop")

## x <- res%>%filter(p.adjusted<0.1)%>%group_by(MCls)%>%nest()%>%
##     mutate(ngene=map(data, ~(.x)$gene))



## ###########################
## ### if enriched in DEGs ###
## ###########################

## atac <- read_rds("../1_processing/4.2_Integrate.outs/3_scATAC.annot.rds")
## DefaultAssay(atac) <- "ATAC"
## peakAnno <- ClosestFeature(atac, regions=granges(atac))
## peakAnno <- peakAnno%>%dplyr::select(gene_name,query_region)

## ### differential results
## res <- read_rds("./1_DiffRNA_2.outs/2.0_DESeq.results.rds")%>%as.data.frame()
## res <- res%>%mutate(MCls=gsub("NK", "NKcell", MCls))
## res <- res%>%dplyr::rename("peak_region"="gene")%>%drop_na(p.value)
## res <- res%>%
##    left_join(peakAnno,by=c("peak_region"="query_region"))%>%
##    mutate(comb=paste(MCls, contrast, sep="_"))

## ### previous identified DEGs
## fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
## resDE <- read_rds(fn)%>%filter(qval<0.1,abs(beta)>0.5)
## x <- bitr(unique(resDE$gene),fromType="ENSEMBL",toType=c("SYMBOL"),OrgDb=org.Hs.eg.db)
## resDE <- resDE%>%
##    inner_join(x,by=c("gene"="ENSEMBL"))%>%
##    mutate(comb=paste(MCls, contrast, sep="_"))


## ###
## comb <- sort(unique(resDE$comb))
## dfNew <- map_dfr(comb, function(ii){
##   res2 <- res%>%filter(comb==ii)
##   DEG <- resDE%>%filter(comb==ii)%>%dplyr::pull(SYMBOL)
##   res2 <- res2%>%mutate(is_DEG=ifelse(gene_name%in%DEG,1,0))
##   ###
##   dx <- map_dfr(c(0,1),function(i){
##      di <- res2%>%filter(is_DEG==i)%>%arrange(p.value)
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

## figfn <- "./1_DiffRNA_2.outs/Figure1.3_qq.DEG.png"
## png(figfn, width=750, height=750, res=120)
## print(p1)
## dev.off()
