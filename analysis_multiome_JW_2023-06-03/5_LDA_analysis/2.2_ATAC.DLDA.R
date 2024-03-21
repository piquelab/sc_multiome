##
library(Matrix)
library(tidyverse)
###
library(data.table)
library(irlba)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(Signac)
library(SeuratWrappers)
library(SeuratObject)

rm(list=ls())

outdir <- "./2_DLDA.outs/"
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)


#### DARs

fn <- "../sc_multiome_data/2_Differential/1.2_DiffPeak.outs/2.0_DESeq.results_0.1_cn_sampleID+new_treat.rds"
res <- read_rds(fn)%>%as.data.frame()%>%mutate(comb=paste(MCls, contrast, sep="_"))

geneSel <- res%>%filter(p.adjusted<0.1, abs(estimate)>0.5)%>%pull(gene)%>%unique()

res2 <- res%>%filter(gene%in%geneSel)%>%
    mutate(is_sig=ifelse(abs(estimate)>0.5&p.adjusted<0.1, 1, 0))



contrast.list <- list("caffeine"=c("caffeine", "control"),
                      "nicotine"=c("nicotine", "control"),
                      "vitA"=c("vitA", "control"),
                      "vitD"=c("vitD", "control"),
                      "vitE"=c("vitE", "control"),
                      "zinc"=c("zinc", "control"))




##################################
### calculate DLDA
##################################
 

MCls <- sort(unique(res$MCls))
for (oneMCl in MCls){
    ##
    fn <- paste("./1_input_data/ATAC_outs/", oneMCl, ".atac.rds", sep="")
    sc <- read_rds(fn)
    DefaultAssay(sc) <- "ATAC"
    sc2 <- subset(sc, features=geneSel)
    ##
    
    data2 <- sc2[["ATAC"]]$data
    data2 <- as.matrix(data2)
    meta2 <- sc2@meta.data%>%dplyr::select(barcode=NEW_BARCODE, MCls, treat2)
    ###
    cat(oneMCl, identical(meta2$barcode, colnames(data2)), "\n")

    ####
    #### normalize data by each gene 
    sd0 <- apply(data2, 1, sd)
    iSel <- !is.na(sd0)&sd0>0
 
    data2 <- data2[iSel,]
    mu <- apply(data2, 1, mean)
    sd0 <- apply(data2, 1, sd)

    
    Xc <- sweep(data2, 1, mu, "-")
    Xc <- sweep(data2, 1, sd0, "/")
    gene <- rownames(Xc)

    ###
    ### for each treatments
    metaNew <- map_dfr(names(contrast.list), function(ii){
       ##       
       oneX <- contrast.list[[ii]] 
        
       ##
       res0 <- res2%>%filter(MCls==oneMCl, contrast==ii)
       s0 <- res0$is_sig
       names(s0) <- res0$gene

       s0 <- s0[gene]
       s0[is.na(s0)] <- 0
        
        
       ## treat2
       cellSel <- meta2%>%filter(treat2==oneX[1])%>%pull(barcode)
       x1 <- data2[, cellSel]
       mu1 <- apply(x1, 1, mean)
        
       ## contrast
       cellSel <- meta2%>%filter(treat2==oneX[2])%>%pull(barcode)
       x2 <- data2[, cellSel]
       mu2 <- apply(x2, 1, mean)

       diff <- as.matrix((mu1-mu2)*(1/sd0)*s0)

       zscore <- as.vector(crossprod(Xc, diff))
       meta2$zscore <- zscore
       meta2$LDA <- paste(oneMCl, ii, sep="_")
       rownames(meta2) <- NULL 
       meta2
    })

    metaNew <- as.data.frame(metaNew)

    opfn <- paste(outdir, "2_ATAC.", oneMCl, ".DLDA.rds", sep="")
    write_rds(metaNew, file=opfn)
}
    
        

        
