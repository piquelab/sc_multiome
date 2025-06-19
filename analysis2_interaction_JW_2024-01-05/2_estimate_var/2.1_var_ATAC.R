###
library(Matrix)
library(tidyverse)
library(limma)
library(edgeR)
library(DESeq2)
library(variancePartition)
library(lme4)

rm(list=ls())

args <- commandArgs(trailingOnly=TRUE)
if ( length(args)>0){
    geneFile <- args[1]
}else{
    geneFile <- "splitPeak000"
}    

outdir <- "./2_vars_ATAC.outs/results/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


###########################################
### prepare input data
###########################################

fn <-"./2_vars_ATAC.outs/1_cv.rds"
cvt2 <- read_rds(fn)
###
fn2 <- "./2_vars_ATAC.outs/1_geneExpr.rds"
geneExpr <- read_rds(fn2)

print(identical(cvt2$bti, colnames(geneExpr)))
      
## genes <- rownames(geneExpr)
## write.table(genes, file="./1_vars_gene.outs/geneAll.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)



#####################################
### design matrix
#####################################

### cell types
m1 <- model.matrix(~0+MCls, data=cvt2)
mm1 <- tcrossprod(m1)
s1 <- mean(diag(mm1))
### treatment
m2 <- model.matrix(~0+treat2, data=cvt2)
mm2 <- tcrossprod(m2)
s2 <- mean(diag(mm2))
### individual
m3 <- model.matrix(~0+sampleID, data=cvt2)
mm3 <- tcrossprod(m3)
s3 <- mean(diag(mm3))

ss <- c(s1, s2, s3, 1)
names(ss) <- c("MCls", "treat2", "sampleID", "Residual")



#######################################################################################
#### Estimate variance explained by cell type, treatments and sampleID
#######################################################################################


fn <- paste0("./2_vars_ATAC.outs/", geneFile)
genes <- read.table(fn, header=F)$V1
ngene <- length(genes)

### results from lmer
res_lmer <- map_dfr(genes, function(ii){
###
    cat(ii, "\n")    
    DF <- data.frame(y=as.vector(geneExpr[ii,]), MCls=cvt2$MCls, treat2=cvt2$treat2, sampleID=cvt2$sampleID)
    lm0 <- lmer(y~(1|MCls)+(1|treat2)+(1|sampleID), data=DF)
    vc_df <- as.data.frame(VarCorr(lm0))%>%
       dplyr::select(grp, vcov)%>%mutate(ss=ss[grp], vcov2=vcov*ss, prop=vcov2/sum(vcov2))
###    
    res0 <- data.frame(t(vc_df$prop))
    res0$gene <- ii
    res0    
})

opfn <- paste(outdir, geneFile, "_results.txt", sep="")
write.table(res_lmer, file=opfn, row.names=FALSE, col.names=FALSE, quote=FALSE)

###
### END

