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

outdir <- "./1_DiffPeak.outs/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)

###
### Differential peaks analysis for raw features, that were not re-call peaks grouped for each cell-types
### Last motified by Julong Wei



####################################
### 1. Generate pseudo-bulk data ###
####################################
combined <- read_rds("../1_processing/3_Clustering.outs/1_seurat.cluster.combined.rds")

colnames(combined@meta.data)

x <- combined@assays$ATAC@ranges

x

xrange <- ranges(x)

xrange

count <- combined@assays$ATAC@counts

head(count)

anno <- data.frame(rn=rownames(count), rnz=rowSums(count),
   chr=as.character(seqnames(x)), start=start(xrange), end=end(xrange))


## anno2 <- data.frame(rn=rownames(count), rnz=rowSums(count), chr=as.character(seqnames(x)))
## annoSel2 <- anno2%>%dplyr::filter(rnz>0, chr%in%autosome)
## Y2 <- count[annoSel2$rn,]
## head(Y2)

## anno3 <- data.frame(rn=rownames(count), rnz=rowSums(count))
## annoSel3 <- anno3%>%dplyr::filter(rnz>0)
## Y3 <- count[annoSel3$rn,]
## head(Y3)

head(anno)

autosome <- paste("chr", as.character(1:22), sep="")

autosome

annoSel <- anno%>%dplyr::filter(rnz>0, chr%in%autosome)

head(annoSel)

Y <- count[annoSel$rn,]

head(Y)

colnames(combined@meta.data)

## meta <- atac@meta.data%>%
##    mutate(treat=gsub(".*-ATAC-|_.*", "", NEW_BARCODE)) 

meta <- combined@meta.data

meta$treat <- meta$treats

colnames(meta)

# cell <- c("0"="Tcell", "1"="Tcell", "2"="NKcell", "3"="Bcell",
#           "4"="Tcell", "5"="Tcell", "6"="Monocyte", "7"="Tcell", "8"="Tcell",
#           "9"="Tcell", "10"="Tcell", "11"="DC", "12"="DC")


cell <- c("0"="CD4Naive", "1"="TCM", "2"="NKcell", "3"="Bcell",
          "4"="CD8Naive", "5"="TEM", "6"="Monocyte", "7"="CTL", "8"="dnT",
          "9"="Treg", "10"="Platelet", "11"="DC", "12"="HSPC")


meta$MCls <- cell[as.character(meta$seurat_cluster)]

meta <- meta%>%
   mutate(bti=paste(MCls, treat, SNG.BEST.GUESS, sep="_"))%>%
   dplyr::select(NEW_BARCODE, bti)


colnames(meta)

dd <- meta%>%group_by(bti)%>%summarise(ncell=n(),.groups="drop")

dd

write_rds(dd, "./1_DiffPeak.outs/0_ncell.rds")



##pseudo-bulk peak data
bti <- factor(meta$bti)

head(bti)

X <- model.matrix(~0+bti)

head(X)

YtX <- Y %*% X

YtX <- as.matrix(YtX)

head(YtX)

colnames(YtX) <- gsub("^bti", "", colnames(YtX))

colnames(YtX)

###rnz>0,chr:1-22, 111,750*400
opfn <- "./1_DiffPeak.outs/1_YtX.comb.rds" 
write_rds(YtX, file=opfn)





###rnz>100, ncell>20, 111,746 * 331
#use when starting analysis again
#dd <- read_rds("./1_DiffPeak.outs/0_ncell.rds")%>%filter(ncell>20)

#YtX <- read_rds("./1_DiffPeak.outs/1_YtX.comb.rds")
#head(YtX)

dd <- dd %>% filter(ncell>20)

anno <- data.frame(rn=rownames(YtX), rnz=rowSums(YtX))
head(anno)

annoSel <- anno%>%filter(rnz>100)
head(annoSel)

YtX_sel <- YtX[annoSel$rn, dd$bti]

head(YtX_sel)

opfn <- "./1_DiffPeak.outs/1_YtX.sel.rds"
write_rds(YtX_sel, file=opfn)



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
#YtX <- read_rds("./1_DiffPeak.outs/1_YtX.sel.rds")

YtX <- YtX_sel

str(YtX)

bti2 <- colnames(YtX)
bti2

cvt0 <- str_split(bti2, "_", simplify=T)
cvt0

cvt <- data.frame(bti=bti2, MCls=cvt0[,1], treat=cvt0[,2], sampleID=cvt0[,3])
cvt

table(cvt$treat)
table(cvt$MCls)
table(cvt$MCls, cvt$treat)

comb <- table(cvt$MCls, cvt$treat)

write.csv(comb, paste0("./1_DiffPeak.outs/combinations.csv"), quote = FALSE, row.names = TRUE)

comb <- read.csv("./1_DiffPeak.outs/combinations.csv")

comb

dd

dd2 <- dd %>% dplyr::filter(bti %in% bti2)
dd2
sum(dd2$ncell)


######################
### 2.1 call DESeq ###
######################
contrast.list <- list("caffeine"=c("treat", "caffeine", "water"),
                      "nicotine"=c("treat", "nicotine", "water"),
                      "vitA"=c("treat", "vitA", "etOH"),
                      "vitD"=c("treat", "vitD", "etOH"),
                      "vitE"=c("treat", "vitE", "etOH"),
                      "zinc"=c("treat", "zinc", "water"),
                      "water"=c("treat", "water", "etOH"))

cat("2.1", "run DESeq batch separately", "\n")

#MCls <- c("Bcell", "Monocyte", "NK", "Tcell", "DC")
#MCls <- c("Bcell", "Monocyte")

MCls  <- c("CD4Naive", "TCM", "NKcell", "Bcell",
          "CD8Naive", "TEM", "Monocyte", "CTL", "dnT")


###

res <- map_dfr(MCls, function(oneX){
  time0 <- Sys.time()
###
  cvt0 <- cvt%>%dplyr::filter(MCls==oneX)
  YtX0 <- YtX[,cvt0$bti]
##
  dds <- DESeqDataSetFromMatrix(YtX0, cvt0, ~ sampleID + treat)
  dds <- DESeq(dds)
  res <- contrast.list%>%map(~results(dds, contrast=.x))
  res2 <- tibble(contrast=names(res),MCls=oneX, data=map(res,tidy))%>%unnest(data)  
##
  time1 <- Sys.time()
  elapsed <- difftime(time1, time0, units="mins")
  cat(oneX, elapsed, "Done\n")
  res2
})

res
table(res$MCls)
table(res$contrast) 

opfn <- "./1_DiffPeak.outs/2.0_DESeq.results_clusters_treatandind.rds"
write_rds(res, opfn)

res <- readRDS("./1_DiffPeak.outs/2.0_DESeq.results_clusters_treatandind.rds")

nrow(res)
nrow(res %>% drop_na(p.value))
nrow(res %>% drop_na(p.adjusted))

######################################################
# significant genes
#####################################################
padj_cutoff <- 0.05

sig_res <- dplyr::filter(res, p.adjusted < padj_cutoff) %>% dplyr::arrange(p.adjusted)

sig_res_split <- split(sig_res, f = list(sig_res$contrast,sig_res$MCls))

sig_res_split2 = list()

for (i in names(sig_res_split)){
    print(i)
    sig_res_split2[[i]] <- sig_res_split[[i]]  %>% dplyr::arrange(p.adjusted) %>% dplyr::pull(gene) %>% head(n=20)
}

sig_res_split2

#######################################################
# table of up and down regulated peaks
#######################################################

#------------------------------------------------------
# Nicole script
#------------------------------------------------------
FDR <- 0.1
pos_foldchange <- 0.25
neg_foldchange <- -0.25


res[!is.na(res$p.adjusted) & res$p.adjusted<=FDR & res$estimate>=pos_foldchange, "Significance"] <- "Upregulated"

res[!is.na(res$p.adjusted) & res$p.adjusted<=FDR & res$estimate<=neg_foldchange, "Significance"] <- "Downregulated"

res[!is.na(res$p.adjusted) & res$p.adjusted>FDR,"Significance"] <- "Not Significant"

res[is.na(res$p.adjusted),"Significance"] <- "NA"

#res2 <- res %>% dplyr::filter(Significance %in% c("Upregulated", "Downregulated"))

res_up <- res %>% dplyr::filter(Significance %in% "Upregulated")

res_down <- res %>% dplyr::filter(Significance %in% "Downregulated")

up <- table(res_up$MCls, res_up$contrast)
up
write.csv(up, paste0("./1_DiffPeak.outs/up_treatandind.csv"), quote = FALSE, row.names = TRUE)

down <- table(res_down$MCls, res_down$contrast)
down
write.csv(down, paste0("./1_DiffPeak.outs/down_treatandind.csv"), quote = FALSE, row.names = TRUE)

# opening files
up <- read.csv(file="./1_DiffPeak.outs/up_treatandind.csv")
down <- read.csv(file="./1_DiffPeak.outs/down_treatandind.csv")

up
down











#-----------number of cells------------------------------------------------
dd2

cvt

dd3 <- left_join(cvt, dd2)
dd3

dd3 <- dd3 %>% dplyr::filter(dd3$MCls %in% res$MCls, !dd3$treat %in% "etOH")

dd4 <- dd3 %>% group_by(treat, MCls) %>% summarise(ncell_count = sum(ncell))

ncell_table <- data.frame(pivot_wider(dd4, names_from=treat, values_from=ncell_count))

rownames(ncell_table) <- ncell_table$MCls
ncell_table$MCls <- NULL

ncell_table

write.csv(ncell_table, paste0("./1_DiffPeak.outs/ncell_treat$ind.csv"), quote = FALSE, row.names = TRUE)

ncell <- read.csv("./1_DiffPeak.outs/ncell_treat$ind.csv")

ncell


#########################################################
# matrix_plots
#########################################################
up
class(up)
dim(up)

#up.df <- as.data.frame(up)
#library(plot.matrix)
## y<-matrix(up)
## y
## class(y)
## par(mar=c(5.1, 4.1, 4.1, 4.1)) # adapt margins
## plot(x)

library(reshape2)
#up
rownames(up)<- up$X
up$X <- NULL
#up
#colnames(up)
#rownames(up)
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
rownames(comb)<- comb$X
comb$X <- NULL
comb.2 <- subset(comb, select = -c(etOH))
comb.2 <- comb.2 %>% dplyr::filter(rownames(comb.2) %in% rownames(up))
#comb.2 <- comb[-c("Platelet", "Treg")]
#comb.2
comb.2.matrix <- as.matrix(comb.2)
comb.2.melt <- melt(comb.2.matrix)


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

figure <- ggarrange(comb.2.p, ncell.p, up.p, down.p + font("x.text", size = 12),
                    ncol = 2, nrow = 2)


png("./1_DiffPeak.outs/plot_all_treatandind.png", width=2000, height=2000, pointsize=14, res=175)
print(figure)
dev.off()







###
### 1. MA plots
figfn <- "./1_DiffPeak.outs/Figure1.1_MA_clusters_treatandind.png"
png(figfn, width=3500, height=4500, pointsize=16, res=225)
par(mar=c(4,4,2,2),mgp=c(2,1,0))
x <- matrix(1:63, 9, 7, byrow=T)
layout(x)
fn <- "./1_DiffPeak.outs/2.0_DESeq.results_clusters_treatandind.rds"
res <- read_rds(fn)%>%
       mutate(color=ifelse((p.adjusted<0.1)&(!is.na(p.adjusted)), T, F))
#MCls <- c("Bcell",  "Monocyte", "NK", "Tcell")
MCls  <- c("CD4Naive", "TCM", "NKcell", "Bcell",
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
figfn <- "./1_DiffPeak.outs/Figure1.2_qq_clusters_treatandind.png"
png(figfn, width=3500, height=4500, pointsize=16, res=225)
par(mar=c(4,4,2,2),mgp=c(2,1,0))
x <- matrix(1:63, 9, 7, byrow=T)
layout(x)
fn <- "./1_DiffPeak.outs/2.0_DESeq.results_clusters_treatandind.rds"
res <- read_rds(fn)%>%drop_na(p.value)
#MCls <- c("Bcell", "Monocyte", "NK", "Tcell")
MCls  <- c("CD4Naive", "TCM", "NKcell", "Bcell",
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





#------------------plot_loop--------------------------------
#colors()

nrow(res) #if qq plot has been called before this step, this already is with dropped NA values

res.drop.na <- res %>% drop_na(p.value)

nrow(res.drop.na)

max.y <- floor(max(-log10(res.drop.na$p.value)))/2

max.y

png("./1_DiffPeak.outs/Figure1.2_qq_clusters_together_treatandind.png", width=2000, height=5000, pointsize=16, res=200)
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
        legend("topleft", #inset=c(-0.2,0),
                          cex = 1,
                          pch = 19,
                          c("caff","nic", "zinc", "vitA", "vitD", "vitE", "water"),
                          fill=c("red","tan", "maroon3", "tan4", "seagreen4", "salmon3", "orange"))
    }                          
dev.off()








#----------plot function-individual plots----------------------------------------
res

MCls <- c("CD4Naive", "TCM", "NKcell", "Bcell",
                    "CD8Naive", "TEM", "Monocyte", "CTL", "dnT")
#MCls <- c("Bcell")

for (oneMCl in MCls){
        png(paste("./1_DiffPeak.outs/Figure1.2_qq_clusters_", oneMCl, "_together_y30.png", sep=""), pointsize=12, width=960, height=960, res=200)
            ##1
            res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="caffeine")%>% drop_na(p.value)
            plot(-log10(1:length(res2$p.value)/length(res2$p.value)),
                          -log10(sort(res2$p.value)),
                          main=oneMCl,
                          pch = 19,
                          ylim = c(0, 30),
                          ylab="observed -log10(p)",
                          xlab="expected -log10(p)",
                          col="red")
            ##2
            res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="nicotine")%>% drop_na(p.value)
            points(-log10(1:length(res2$p.value)/length(res2$p.value)),-log10(sort(res2$p.value)),  col="grey", pch = 19, ylim = c(0, 30))
            ##3
            res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="zinc")%>% drop_na(p.value)
            points(-log10(1:length(res2$p.value)/length(res2$p.value)),-log10(sort(res2$p.value)),  col="green", pch = 19, ylim = c(0, 30))
            ##4
            res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="vitA")%>% drop_na(p.value)
            points(-log10(1:length(res2$p.value)/length(res2$p.value)),-log10(sort(res2$p.value)),  col="blue", pch = 19, ylim = c(0, 30))
            ##5
            res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="vitD")%>% drop_na(p.value)
            points(-log10(1:length(res2$p.value)/length(res2$p.value)),-log10(sort(res2$p.value)),  col="pink", pch = 19, ylim = c(0, 30))
            ##6
            res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="vitE")%>% drop_na(p.value)
            points(-log10(1:length(res2$p.value)/length(res2$p.value)),-log10(sort(res2$p.value)),  col="yellow", pch = 19, ylim = c(0, 30))
            ##7
            res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="water")%>% drop_na(p.value)
            points(-log10(1:length(res2$p.value)/length(res2$p.value)),-log10(sort(res2$p.value)),  col="orange", pch = 19, ylim = c(0, 30))
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





###################################
### opening seurat object ###
###################################
combined <- read_rds("../1_processing/3_Clustering.outs/1_seurat.cluster.combined.rds")

combined[["ATAC"]]
combined[["RNA"]]

## #-------------------------------
## x <- combined@assays$ATAC@ranges
## x
## xrange <- ranges(x)
## xrange
## count <- combined@assays$ATAC@counts
## anno <- data.frame(rn=rownames(count), rnz=rowSums(count),
##                       chr=as.character(seqnames(x)), start=start(xrange), end=end(xrange))
## head(anno)
## autosome <- as.character(1:22)
## autosome
## annoSel <- anno%>%dplyr::filter(rnz>0, chr%in%autosome)


#---------------------------------------
# get gene annotations for hg38
#---------------------------------------
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

## chrs = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86, standard.chromosomes = TRUE)

head(annotation)

seqlevelsStyle(annotation) <- "UCSC"

head(annotation)

## autosome <- as.character(1:22)

## head(annotation)
## unique(annotation$seqnames)
## colnames(annotation)
## class(annotation)
## #annotation.2 <- as.data.frame(annotation)
## annotation.2 <- annotation%>%dplyr::filter(seqnames%in%autosome)
## Annotation(combined[["ATAC"]]) <- annotation







## ###################################
## ### Table of Differential peaks ###
## ###################################
## fn <- "./1_DiffPeak.outs/2.0_DESeq.results.rds"

## res.dropna.pvalue <- res %>%drop_na(p.value)

## res.dropna.pvalue %>%filter(p.adjusted<0.1,abs(estimate)>0.25)%>%
##    group_by(MCls, contrast)%>%
##    summarise(ngene=n(),.groups="drop")

## x <- res%>%filter(p.adjusted<0.1)%>%group_by(MCls)%>%nest()%>%
##     mutate(ngene=map(data, ~(.x)$gene))

## x





###########################
### if enriched in DEGs ###
###########################
#atac <- read_rds("../1_processing/4.2_Integrate.outs/3_scATAC.annot.rds")

combined <- read_rds("../1_processing/3_Clustering.outs/1_seurat.cluster.combined.rds")

combined

combined[["ATAC"]]

DefaultAssay(combined) <- "ATAC"

granges(combined[["ATAC"]])

#Annotation(combined) <- annotation

gene.ranges = genes(EnsDb.Hsapiens.v86)
#head(gene.ranges)

seqlevelsStyle(gene.ranges) <- 'UCSC'
#head(gene.ranges)

peakAnno <- ClosestFeature(combined,
                           regions=granges(combined),
                           annotation = gene.ranges,
                           sep = c(':', '-')
                           )

head(peakAnno)
colnames(peakAnno)

peakAnno <- peakAnno%>%dplyr::select(gene_name,query_region)%>%mutate(query_region=gsub(":", "-", query_region))

head(peakAnno)








# res object
#---------------------------------------------
### differential results
#res <- read_rds("./1_DiffPeak.outs/2.0_DESeq.results.rds")%>%as.data.frame()

res <- readRDS("./1_DiffPeak.outs/2.0_DESeq.results_clusters_treatandind.rds")

res

nrow(res)

#table(res$MCls)
#res <- res%>%mutate(MCls=gsub("NK", "NKcell", MCls))

res <- res%>%dplyr::rename("peak_region"="gene")%>%drop_na(p.value)

res

res <- res%>%
    left_join(peakAnno,by=c("peak_region"="query_region"))%>%
    mutate(comb=paste(MCls, contrast, sep="_"))

res

head(res$gene_name)


# DEGs
#---------------------------------------------
### previous identified DEGs
#fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew_output/Filter2/2_meta.rds"

fn <- "./1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind.rds"
resDE <- read_rds(fn)
#resDE
#resDE <- resDE %>% filter(qval<0.1,abs(beta)>0.5)

resDE <- resDE %>% dplyr::filter(p.adjusted<0.1,abs(estimate)>0.25)

#unique(resDE$gene)



# ENSEMBL to symbol
#-------------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)

x <- bitr(unique(resDE$gene),fromType="ENSEMBL",toType=c("SYMBOL"),OrgDb=org.Hs.eg.db)

#x <- bitr(unique(resDE$gene),fromType="GENENAME",toType=c("SYMBOL"),OrgDb=org.Hs.eg.db)

#x <- bitr(unique(resDE$gene),fromType="GENENAME",toType=c("SYMBOL"),EnsDb=EnsDb.Hsapiens.v86)

head(x)

#showMethods(keys)
#help("SYMBOL")

#resDE <- resDE%>%
#   inner_join(x,by=c("gene"="ENSEMBL"))%>%
#   mutate(comb=paste(MCls, contrast, sep="_"))


resDE <- resDE%>%mutate(comb=paste(MCls, contrast, sep="_"))
#head(resDE)


# trial
#---------------------------------------------
table(res$gene_name)

resDE.genes <- resDE$gene

head(resDE.genes)

new_df <- res %>% dplyr::filter(gene_name %in% resDE.genes)

new_df


# function
###-------------------------------------------
comb <- sort(unique(resDE$comb))
#comb

res2 <- res%>%dplyr::filter(comb=="TEM_nicotine")
res2

DEG <- resDE%>%dplyr::filter(comb=="TEM_nicotine")%>%dplyr::pull(gene)
DEG

res2 <- res2%>%mutate(is_DEG=ifelse(gene_name%in%DEG,1,0))
table(res2$is_DEG)

di <- res2%>%dplyr::filter(is_DEG==1)%>%arrange(p.value)

head(di$is_DEG)

ngene <- nrow(di)

ngene

di <- di%>%mutate(observed=-log10(p.value), expected=-log10(ppoints(ngene)))

di


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

dfNew <- map_dfr(comb, function(ii){
  res2 <- res%>%dplyr::filter(comb==ii)
  DEG <- resDE%>%dplyr::filter(comb==ii)%>%dplyr::pull(gene)
  res2 <- res2%>%mutate(is_DEG=ifelse(gene_name%in%DEG,1,0))
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

table(dfNew$contrast)



### plot
###-----------------------------------------------
#lab1 <- c(unique(dfNew$contrast))
#labl
lab1 <- c("caffeine"="caffeine", "nicotine"="nicotine", "vitA"="vitA", "vitD"="vitD", "vitE"="vitE", "water"="water", "zinc"="zinc")
p1 <- ggplot(dfNew, aes(x=expected,y=observed, colour=factor(is_DEG)))+
   ggrastr::rasterise(geom_point(size=0.3),dpi=300)+
   geom_abline(colour="red")+
   scale_colour_manual(values=c("0"="grey40","1"="green"),
                       labels=c("0"="Not DEG", "1"="DEG"),
                       guide=guide_legend(override.aes=list(size=1)))+
   facet_grid(MCls~contrast, scales="free", labeller=labeller(contrast=lab1))+
   xlab(bquote("Expected"~ -log[10]~"("~italic(plain(P))~")"))+
   ylab(bquote("Observed"~-log[10]~"("~italic(plain(P))~")"))+
   theme_bw()+
   theme(legend.title=element_blank(),strip.text=element_text(size=12))

figfn <- "./1_DiffPeak.outs/Figure1.3_qq.DEG.png"
png(figfn, width=750, height=750, res=120)
print(p1)
dev.off()
