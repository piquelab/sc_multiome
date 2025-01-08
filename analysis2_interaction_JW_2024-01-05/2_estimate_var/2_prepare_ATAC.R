###
library(Matrix)
library(tidyverse)
library(limma)
library(edgeR)
library(lme4)
library(DESeq2)
library(variancePartition)


## library(qqman)
## library(qvalue)
## library(biobroom)
## library(ashr)
## library(GenomicRanges)
## library(Seurat)
## library(SeuratDisk)
## library(SeuratData)
## library(Signac)
## library(clusterProfiler)
## library(org.Hs.eg.db)
## library(annotables)

library(ComplexHeatmap)
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(ggExtra)
library(RColorBrewer)
library(ggsci)
library(circlize)
library(viridis)
library(openxlsx)

rm(list=ls())

outdir <- "./2_vars_ATAC.outs/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)




#############################
### Normalize data 
#############################

YtX <- read_rds("../sc_multiome_data_MHB/2_Differential/1.2_DiffPeak.outs/1_YtX.sel_0.1_cn.rds")


###
### meta data
bti2 <- colnames(YtX)
#length(bti2) 
cvt0 <- str_split(bti2, "_", simplify=T)
cvt <- data.frame(bti=bti2, MCls=paste(cvt0[,1], cvt0[,2], sep="_"), treat=cvt0[,3], sampleID=cvt0[,4])

### filtering cell type_ind with 8 treatment
x <- cvt%>%group_by(MCls, sampleID)%>%mutate(ntreat=n())%>%ungroup()
btiSel <- x%>%dplyr::filter(ntreat==8)%>%pull(bti)

cvt2 <- cvt%>%dplyr::filter(bti%in%btiSel)

cvt2 <- cvt2%>%mutate(treat2=ifelse(treat%in%c("etOH", "water"), "control", treat))                   
YtX2 <- YtX[,cvt2$bti]



####
data <- YtX2
nsample <- ncol(data)

### create DGEList object, filter and calculate norm.factors
dge <- DGEList(data)
cpm <- cpm(dge) ## using original library size to calculate cpm
keep.exprs <- rowSums(cpm>0.1)>=(0.2*nsample)
dge <- dge[keep.exprs,]

dge <- calcNormFactors(dge) ## effect lib.size=lib.size*norm.factors 
## cpm <- cpm(dge) ### (counts/effect lib.size)*1e+6



### Normalized equal to cpm(log=TRUE)
v <- voom(dge)
geneExpr <- v$E

identical(cvt2$bti, colnames(geneExpr))

###
### output 
opfn <- paste(outdir, "1_cv.rds", sep="")
write_rds(cvt2, file=opfn)
###
opfn2 <- paste(outdir, "1_geneExpr.rds", sep="")
write_rds(geneExpr, file=opfn2)

###
### peaks for final analysis 
peaks <- rownames(geneExpr)
opfn <- paste(outdir, "peakAll.txt", sep="")
write.table(peaks, file=opfn, quote=F, row.names=F, col.names=F)
          
## form <- ~(1|MCls)+(1|treat)+(1|sampleID)
## cv <- cvt2[, -5]%>%column_to_rownames(var="bti")
## res <- fitExtractVarPartModel(geneExpr, form, cv)



######################################################################################################
### Estimate variance explained by cell types, treatments, individual and residual
######################################################################################################

###
## Example genes

fn <-"./2_vars_ATAC.outs/1_cv.rds"
cvt2 <- read_rds(fn)
###
fn2 <- "./2_vars_ATAC.outs/1_geneExpr.rds"
geneExpr <- read_rds(fn2)

## print(identical(cvt2$bti, colnames(geneExpr)))
      

## ### cell types
## m1 <- model.matrix(~0+MCls, data=cvt2)
## mm1 <- crossprod(m1)
## s1 <- mean(diag(mm1))
## ### treatment
## m2 <- model.matrix(~0+treat2, data=cvt2)
## mm2 <- crossprod(m2)
## s2 <- mean(diag(mm2))
## ### individual
## m3 <- model.matrix(~0+sampleID, data=cvt2)
## mm3 <- crossprod(m3)
## s3 <- mean(diag(mm3))

## ss <- c(s1, s2, s3, 1)
## names(ss) <- c("MCls", "treat2", "sampleID", "Residual")



## genes <- rownames(geneExpr)

## ### results from lmer 
## ii <- genes[1]
## DF <- data.frame(y=as.vector(geneExpr[ii,]), MCls=cvt2$MCls, treat2=cvt2$treat2, sampleID=cvt2$sampleID)
## lm0 <- lmer(y~(1|MCls)+(1|treat2)+(1|sampleID), data=DF)
## vc_df <- as.data.frame(VarCorr(lm0))%>%
##     dplyr::select(grp, vcov)%>%mutate(ss=ss[grp], vcov2=vcov*ss, prop=vcov2/sum(vcov2))
## df0 <- data.frame(t(vc_df$prop))
## df0$gene <- ii



fn <- "./2_vars_ATAC.outs/all_peaks.vars.txt"
res <- read.table(fn)
names(res) <- c("sampleID", "MCls", "treat2", "Residual", "peaks")
##
opfn <- paste(outdir, "all_peaks.vars.txt", sep="")
write.table(res, file=opfn, quote=FALSE, sep="\t", row.names=FALSE)
