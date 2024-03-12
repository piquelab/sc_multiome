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


#### DEGs
fn <- "../sc_multiome_data/2_Differential/1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn_sampleID+new_treat.rds"
resDG <- read_rds(fn)%>%as.data.frame()%>%mutate(comb=paste(MCls, contrast, sep="_"))
DEG <- resDG%>%dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)%>%pull(gene)%>%unique()

resDG2 <- resDG%>%filter(gene%in%DEG)%>%
    mutate(is_sig=ifelse(abs(estimate)>0.5&p.adjusted<0.1, 1, 0))

###resDG <- resDG%>%filter(gene%in%DEG)


contrast.list <- list("caffeine"=c("caffeine", "control"),
                      "nicotine"=c("nicotine", "control"),
                      "vitA"=c("vitA", "control"),
                      "vitD"=c("vitD", "control"),
                      "vitE"=c("vitE", "control"),
                      "zinc"=c("zinc", "control"))




##################################
### calculate DLDA
##################################


MCls <- sort(unique(resDG$MCls))
for (oneMCl in MCls){
    ##
    fn <- paste("./1_input_data/", oneMCl, ".multiome.rds", sep="")
    sc2 <- read_rds(fn)
    ##
    data2 <- sc2[["SCT"]]$data
    meta2 <- sc2@meta.data%>%dplyr::select(barcode, MCls, treat2)

    identical(meta2$barcode, colnames(data2))
    cat(oneMCl, "\n")

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
       res0 <- resDG2%>%filter(MCls==oneMCl, contrast==ii)
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
       meta2
    })

    metaNew <- as.data.frame(metaNew)

    opfn <- paste(outdir, "1_RNA.", oneMCl, ".DLDA.rds", sep="")
    write_rds(metaNew, file=opfn)
}
    
        

        
