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
library(reshape2)
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
YtX_sel <- read_rds("./1_DiffRNA_2.outs/1_YtX.sel_clusters_0.1_cn.rds")
#YtX_sel <- read_rds("./1.2_DiffPeak.outs/1_YtX.sel_0.1_cn.rds")


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
#head(cvt)
cvt <- cvt%>%mutate(treat_IDs=paste(cvt$treat, ".", cvt$sampleID, sep=""))
#head(cvt)
cvt <- cvt%>%mutate(treat_MCls=paste(cvt$treat, ".", cvt$MCls, sep=""))
head(cvt)

MCls.IDs <- unique(cvt$MCls_IDs)
MCls.IDs

treat.IDs <- unique(cvt$treat_IDs)
treat.IDs

treat.MCls <- unique(cvt$treat_MCls)
treat.MCls








## # #################################################################
## # ### 2.1 call DESeq ###
## # #################################################################
## res.DE <- readRDS("./1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn.rds")

## res.DE.contrastetOH <- readRDS("./1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn_sampleID+treat_allvsetOH.rds")

## head(res.DE)
## table(res.DE$MCls)
## # #nrow(res.DE)
## # #nrow(res.DE %>% drop_na(p.adjusted))
## # #nrow(res.DE %>% drop_na(p.value))


## # #################################################################
## # ### #----comparisons all celltypes using resDE---- ##1st way
## # #################################################################
## head(res.DE)
## table(res.DE$contrast)
## table(res.DE$MCls)


## #-----------------MCls----------------------
## MCls <- unique(res.DE$MCls)
## MCls

## scorr_MCls <- c()
## i_MCls <- c()
## j_MCls <- c()


## for(i in MCls){
##     for(j in MCls){
##         cvt0 <- res.DE%>%dplyr::filter(MCls==i)
##         cvt1 <- res.DE%>%dplyr::filter(MCls==j)
##         ##
##         scorr <- as.numeric(cor.test(cvt0$estimate, cvt1$estimate, method = "spearman")$estimate)
##         print(scorr)
##         scorr_MCls <- c(scorr_MCls, scorr)
##         i_MCls <- c(i_MCls, i)
##         j_MCls <- c(j_MCls, j)
##         ##data1 <- plotPCA(vd1, intgroup=c("new_treat"), returnData=TRUE)
##         ##print(head(data1))
##         ##data1 <- data1 %>% mutate(sampleID=str_split(data$name, "_", simplify=T)[,4])
##         ##head(data1)
##     }
## }

## head(scorr_MCls)
## head(i_MCls)
## head(j_MCls)





## #-----------------treat----------------------
## treat <- unique(res.DE$contrast)
## treat

## scorr_treat <- c()
## i_treat <- c()
## j_treat <- c()

## for(i in treat){
##     for(j in treat){
##         cvt0 <- res.DE%>%dplyr::filter(contrast==i)
##         cvt1 <- res.DE%>%dplyr::filter(contrast==j)
##         ##
##         scorr <- as.numeric(cor.test(cvt0$estimate, cvt1$estimate, method = "spearman")$estimate)
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







## ## #----------------sampleIDS-----------------------
## ## sample_ID <- unique(res.DE$sample_ID)
## ## sample_ID

## ## scorr_sample_ID <- c()
## ## i_sample_ID <- c()
## ## j_sample_ID <- c()

## ## for(i in sample_ID){
## ##       for(j in sample_ID){
## ##               cvt0 <- res.DE%>%dplyr::filter(sample_ID==i)
## ##                   cvt1 <- res.DE%>%dplyr::filter(sample_ID==j)
## ##                   ###
## ##                   scorr <- as.numeric(cor.test(cvt0$estimate, cvt1$estimate, method = "spearman")$estimate)
## ##                   print(scorr)
## ##                   scorr_sample_ID <- c(scorr_sample_ID, scorr)
## ##                   i_sample_ID <- c(i_sample_ID, i)
## ##                   j_sample_ID <- c(j_sample_ID, j)
## ##                   #data1 <- plotPCA(vd1, intgroup=c("new_treat"), returnData=TRUE)
## ##                   #print(head(data1))
## ##                   #data1 <- data1 %>% mutate(sampleID=str_split(data$name, "_", simplify=T)[,4])
## ##                   #head(data1)
## ##                 }
## ##       }

## ## head(scorr_sample_ID)
## ## head(i_sample_ID)
## ## head(j_sample_ID)



## #----------------plot-----------------------
## #df <- data.frame(scorr_MCls, scorr_treat)

## #p <- boxplot(scorr_MCls, scorr_treat)
## png(paste0("./1_DiffRNA_2_3_corr.outs/corr_resDE.png"), width=500, height=500, pointsize=16, res=125)
## boxplot(scorr_MCls, scorr_treat)
## dev.off()







## # #########################################################################
## # ### #----comparisons all celltypes using resDE as in paper---- ##2nd way
## # #########################################################################
## #-----------------MCls----------------------
## treat <- unique(res.DE$contrast) 
## MCls <- unique(res.DE$MCls)
## MCls

## scorr_MCls <- c()
## pval_MCls <- c()
## k_MCls <- c()
## i_MCls <- c()
## j_MCls <- c()

## for(k in treat){
##     res.DE2 <- res.DE%>%dplyr::filter(contrast==k)
##     for(i in MCls){
##         res.DE3 <- res.DE2%>% dplyr::filter(MCls==i)
##         for(j in MCls){
##             res.DE4 <- res.DE2%>%dplyr::filter(MCls==j)
##             ##
##             scorr <- as.numeric(cor.test(res.DE3$estimate, res.DE4$estimate, method = "spearman")$estimate)
##             pval <- as.numeric(cor.test(res.DE3$estimate, res.DE4$estimate, method = "spearman")$p.value)
##             print(scorr)
##             scorr_MCls <- c(scorr_MCls, scorr)
##             pval_MCls <- c(pval_MCls, pval)
##             k_MCls <- c(k_MCls, k)
##             i_MCls <- c(i_MCls, i)
##             j_MCls <- c(j_MCls, j)
##         }
##     }
## }

## head(scorr_MCls)
## head(pval_MCls)
## head(k_MCls)
## head(i_MCls)
## head(j_MCls)

## df.MCls <- data.frame(treat=k_MCls, celltype1=i_MCls, celltype2=j_MCls, corr=scorr_MCls, pval=pval_MCls)
## head(df.MCls)
## nrow(df.MCls %>% dplyr::filter(pval>0.05))


## #-----------------treat----------------------
## #treat <- unique(res.DE$contrast)
## #treat

## scorr_treat <- c()
## pval_treat <- c()
## k_treat <- c()
## i_treat <- c()
## j_treat <- c()

## for(k in MCls){
##     res.DE2 <- res.DE%>%dplyr::filter(MCls==k)
##     for(i in treat){
##         res.DE3 <- res.DE2%>% dplyr::filter(contrast==i)
##         for(j in treat){
##             res.DE4 <- res.DE2%>%dplyr::filter(contrast==j)
##             ##
##             scorr <- as.numeric(cor.test(res.DE3$estimate, res.DE4$estimate, method = "spearman")$estimate)
##             pval <- as.numeric(cor.test(res.DE3$estimate, res.DE4$estimate, method = "spearman")$p.value)
##             print(scorr)
##             scorr_treat <- c(scorr_treat, scorr)
##             pval_treat <- c(pval_treat, pval)
##             k_treat <- c(k_treat, k)
##             i_treat <- c(i_treat, i)
##             j_treat <- c(j_treat, j)
##         }
##     }
## }

## head(scorr_treat)
## head(pval_treat)
## head(k_treat)
## head(i_treat)
## head(j_treat)


## df.treat <- data.frame(celltype=k_treat, treat1=i_treat, treat2=j_treat, corr=scorr_treat, pval=pval_treat)
## head(df.treat)
## nrow(df.treat %>% dplyr::filter(pval>0.05))


## #----------------plot-----------------------
## #p <- boxplot(scorr_MCls, scorr_treat)
## png(paste0("./1_DiffRNA_2_3_corr.outs/corr_resDE_2ndway_eachtreat.png"), width=500, height=500, pointsize=16, res=125)
## boxplot(scorr_MCls, scorr_treat)
## dev.off()

## head(res.DE)





## ######################################################################
## ## resDE - only treatments
## ######################################################################
## ##--treat - caff, nic, zinc vs water--##
## head(res.DE)

## MCls <- unique(res.DE$MCls)
## #treat <- unique(res2$treat)

## treat_resDE_vswater <- c("caffeine", "nicotine", "zinc")

## scorr_treat_resDE_vswater <- c()
## pval_treat_resDE_vswater <- c()
## k_treat_resDE_vswater <- c()
## i_treat_resDE_vswater <- c()
## #j_treat_resDE_vswater <- c()

## for(k in MCls){
##     res.DE2 <- res.DE%>%dplyr::filter(MCls==k)
##     for(i in treat_resDE_vswater){
##         res.DE3 <- res.DE2%>% dplyr::filter(contrast==i)
##         #for(j in treat_resDE_vswater){
##             res.DE4 <- res.DE2%>%dplyr::filter(contrast=="water")
##             ##
##             scorr <- as.numeric(cor.test(res.DE3$estimate, res.DE4$estimate, method = "spearman")$estimate)
##             pval <- as.numeric(cor.test(res.DE3$estimate, res.DE4$estimate, method = "spearman")$p.value)
##             print(scorr)
##             scorr_treat_resDE_vswater <- c(scorr_treat_resDE_vswater, scorr)
##             pval_treat_resDE_vswater <- c(pval_treat_resDE_vswater, pval)
##             k_treat_resDE_vswater <- c(k_treat_resDE_vswater, k)
##             i_treat_resDE_vswater <- c(i_treat_resDE_vswater, i)
## #            j_treat_resDE_vswater <- c(j_treat_resDE_vswater, j)
##     }
## }
## #}

## head(scorr_treat_resDE_vswater)
## head(pval_treat_resDE_vswater)
## head(k_treat_resDE_vswater)
## head(i_treat_resDE_vswater)
## #head(j_treat_resDE_vswater)


## df.treat_resDE_vswater <- data.frame(celltype=k_treat_resDE_vswater, treat1=i_treat_resDE_vswater, treat2="water", corr=scorr_treat_resDE_vswater, pval=pval_treat_resDE_vswater)
## head(df.treat_resDE_vswater)
## nrow(df.treat_resDE_vswater %>% dplyr::filter(pval>0.05))







## ##--treat - vits vs etOH--##
## table(res.DE$MCls)
## table(res.DE$contrast)

## MCls <- unique(res.DE$MCls)

## #treat <- unique(res.DE$contrast)
## #treat

## treat_resDE_vsetOH <- c("vitA", "vitD", "vitE")

## treat_resDE_vsetOH <- c("caffeine", "nicotine", "zinc")

## scorr_treat_resDE_vsetOH <- c()
## pval_treat_resDE_vsetOH <- c()
## k_treat_resDE_vsetOH <- c()
## i_treat_resDE_vsetOH <- c()
## #j_treat_resDE_vsetOH <- c()

## for(k in MCls){
##     res.DE2 <- res.DE.contrastetOH%>%dplyr::filter(MCls==k)
##     for(i in treat_resDE_vsetOH){
##         res.DE3 <- res.DE2%>% dplyr::filter(contrast==i)
##         #for(j in treat_resDE_vsetOH){
##             res.DE4 <- res.DE2%>%dplyr::filter(contrast=="water")
##             ##
##             scorr <- as.numeric(cor.test(res.DE3$estimate, res.DE4$estimate, method = "spearman")$estimate)
##             pval <- as.numeric(cor.test(res.DE3$estimate, res.DE4$estimate, method = "spearman")$p.value)
##             print(scorr)
##             scorr_treat_resDE_vsetOH <- c(scorr_treat_resDE_vsetOH, scorr)
##             pval_treat_resDE_vsetOH <- c(pval_treat_resDE_vsetOH, pval)
##             k_treat_resDE_vsetOH <- c(k_treat_resDE_vsetOH, k)
##             i_treat_resDE_vsetOH <- c(i_treat_resDE_vsetOH, i)
## #            j_treat_resDE_vsetOH <- c(j_treat_resDE_vsetOH, j)
##     }
## }
## #}

## head(scorr_treat_resDE_vsetOH)
## head(pval_treat_resDE_vsetOH)
## head(k_treat_resDE_vsetOH)
## head(i_treat_resDE_vsetOH)
## #head(j_treat_resDE_vsetOH)


## df.treat_resDE_vsetOH <- data.frame(celltype=k_treat_resDE_vsetOH, treat1=i_treat_resDE_vsetOH, treat2="etOH", corr=scorr_treat_resDE_vsetOH, pval=pval_treat_resDE_vsetOH)
## head(df.treat_resDE_vsetOH)
## nrow(df.treat_resDE_vsetOH %>% dplyr::filter(pval>0.05))







## ##--treat - water vs etOH--##
## MCls <- unique(res.DE$MCls)
## #treat <- unique(res.DE$treat)

## treat_resDE_cont <- c("caffeine", "nicotine", "zinc")

## scorr_treat <- c()
## pval_treat_resDE_cont <- c()
## k_treat_resDE_cont <- c()
## i_treat_resDE_cont <- c()
## #j_treat_resDE_cont <- c()

## for(k in MCls){
##     res.DE2 <- res.DE%>%dplyr::filter(MCls==k)
##     for(i in treat_resDE_cont){
##         res.DE3 <- res.DE2%>% dplyr::filter(contrast==i)
##         #for(j in treat_resDE_cont){
##             res.DE4 <- res.DE2%>%dplyr::filter(contrast=="water")
##             ##
##             scorr <- as.numeric(cor.test(res.DE3$estimate, res.DE4$estimate, method = "spearman")$estimate)
##             pval <- as.numeric(cor.test(res.DE3$estimate, res.DE4$estimate, method = "spearman")$p.value)
##             print(scorr)
##             scorr_treat_resDE_cont <- c(scorr_treat_resDE_cont, scorr)
##             pval_treat_resDE_cont <- c(pval_treat_resDE_cont, pval)
##             k_treat_resDE_cont <- c(k_treat_resDE_cont, k)
##             i_treat_resDE_cont <- c(i_treat_resDE_cont, i)
## #            j_treat_resDE_cont <- c(j_treat_resDE_cont, j)
##     }
## }
## #}

## head(scorr_treat_resDE_cont)
## head(pval_treat_resDE_cont)
## head(k_treat_resDE_cont)
## head(i_treat_resDE_cont)
## #head(j_treat_resDE_cont)


## df.treat_resDE_cont <- data.frame(celltype=k_treat_resDE_cont, treat1=i_treat_resDE_cont, treat2="water", corr=scorr_treat_resDE_cont, pval=pval_treat_resDE_cont)
## head(df.treat_resDE_cont)
## nrow(df.treat_resDE_cont %>% dplyr::filter(pval>0.05))






## #----------------plot-----------------------
## #df <- data.frame(scorr_MCls, scorr_treat)

## #p <- boxplot(scorr_MCls, scorr_treat)
## png(paste0("./1_DiffRNA_2_3_corr.outs/corr_resDE_treatments_allcontrastvsetOH_corrvswater.png"), width=500, height=500, pointsize=16, res=125)
## #boxplot(scorr_treat_resDE_vswater, scorr_treat_resDE_vsetOH, scorr_treat_resDE_cont)
## boxplot(scorr_treat_resDE_vswater,  scorr_treat_resDE_vsetOH)
## dev.off()
































# #################################################################
# ### #----comparisons all celltypes as defined in paper----
# #################################################################
head(cvt)

####
cvt.filtered <- map_dfr(MCls.IDs, function(oneX){
    time0 <- Sys.time()
    ##
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
YtX3 <- YtX[,bti3]
##ncol(YtX)


MCls <- unique(cvt.3$MCls)
MCls

#MCls <- "0_CD4Naive"
#MCls2 <- "2_NKcell"

#scorr_MCls <- c()
#i_MCls <- c()
#j_MCls <- c()


## res <- map_dfr(MCls, function(oneX){
## ##for(i in MCls){
##     ##for(j in MCls){
##         cvt0 <- cvt.3%>%dplyr::filter(MCls==oneX)
##         ##cvt1 <- cvt.3%>%dplyr::filter(MCls==j)
##         ##
##         YtX0 <- YtX.MCls[,cvt0$bti]
##         ##YtX1 <- YtX.MCls[,cvt1$bti]
##         ##
##         dds0 <- DESeqDataSetFromMatrix(YtX0, cvt0, ~ sampleID + vehicle + new_treat)
##         ##dds1 <- DESeqDataSetFromMatrix(YtX1, cvt1, ~ sampleID + vehicle + new_treat)
##         ##
##         vd0 <- vst(dds0)
##         ##vd1 <- vst(dds1)
##         ##
##         data0 <- melt(assay(vd0))
##         ##data1 <- melt(assay(vd1))
##         ##print(head(data1))
##         ##print(nrow(data1))
##         ##
##         ##scorr <- as.numeric(cor.test(data0$value, data1$value, method = "spearman")$estimate)
##         ##print(scorr)
##         ##scorr_MCls <- c(scorr_MCls, scorr)
##         ##i_MCls <- c(i_MCls, i)
##         ##j_MCls <- c(j_MCls, j)
##         ##
##         ##data1 <- plotPCA(vd1, intgroup=c("new_treat"), returnData=TRUE)
##         ##print(head(data1))
##         ##data1 <- data1 %>% mutate(sampleID=str_split(data$name, "_", simplify=T)[,4])
##         ##head(data1)
## ##    }
## ##}
## })



res <- map_dfr(MCls, function(oneX){
    cvt0 <- cvt.3%>%dplyr::filter(MCls==oneX)
    YtX0 <- YtX3[,cvt0$bti]
    dds0 <- DESeqDataSetFromMatrix(YtX0, cvt0, ~ sampleID + vehicle + new_treat)
    vd0 <- vst(dds0)
    data0 <- melt(assay(vd0))
})

head(res)
nrow(res)

res2 <- res %>% mutate(MCls=paste0(str_split(res$Var2, "_", simplify=T)[,1], "_", str_split(res$Var2, "_", simplify=T)[,2]), treat=str_split(res$Var2, "_", simplify=T)[,3], sampleID=str_split(res$Var2, "_", simplify=T)[,4]) 
head(res2)
nrow(res2)

table(res2$MCls)
table(res2$treat)
table(res2$sampleID)



##--MCls--##
#MCls <- unique(res2$MCls)
#treat <- unique(res2$treat)
sampleID <- unique(res2$sampleID)
sampleID

#sampleID[1]
#sampleID <- "AL-018"

celltypes <- unique(res2$MCls)
celltypes
celltypes[7]



scorr_MCls <- c()
pval_MCls <- c()
i_MCls <- c()
j_MCls <- c()
k_MCls <- c()
l_MCls <- c()


#for(i in 1:6){
#    print(i)
#    }

for(i in sampleID){
    res3 <- res2 %>% dplyr::filter(sampleID==i)
    for (j in unique(res3$treat)){
        res4 <- res3 %>% dplyr::filter(treat==j)
        celltypes <- unique(res4$MCls)
        #print(celltypes)
        #print(length(celltypes))
        #for (k in unique(res4$MCls)){
        for (k in 1:length(celltypes)){
            print(k)
            res5 <- res4 %>% dplyr::filter(MCls==celltypes[k])
            #print(nrow(res5))
            if(k != length(celltypes)){
                for (l in (k+1):length(celltypes)){
                    #print(length(celltypes))
                    print(l)
                    res6 <- res4 %>% dplyr::filter(MCls==celltypes[l])
                    scorr <- as.numeric(cor.test(res5$value, res6$value, method = "spearman")$estimate)
                    pval <- as.numeric(cor.test(res5$value, res6$value, method = "spearman")$p.value)
                    print(scorr)
                    scorr_MCls <- c(scorr_MCls, scorr)
                    pval_MCls <- c(pval_MCls, pval)
                    i_MCls <- c(i_MCls, i)
                    j_MCls <- c(j_MCls, j)
                    k_MCls <- c(k_MCls, celltypes[k])
                    l_MCls <- c(l_MCls, celltypes[l])
                }
            }
        }
    }
}

head(scorr_MCls)
head(pval_MCls)
head(i_MCls)
head(j_MCls)
head(k_MCls)
head(l_MCls)
                
df.MCls <- data.frame(sampleID=i_MCls, treat=j_MCls, celltype1=k_MCls, celltype2=l_MCls, corr=scorr_MCls, pval=pval_MCls)
head(df.MCls)
nrow(df.MCls %>% dplyr::filter(pval>0.05))




##--treat--##
#MCls <- unique(res2$MCls)
#treat <- unique(res2$treat)
sampleID <- unique(res2$sampleID)
sampleID


scorr_treat <- c()
pval_treat <- c()
i_treat <- c()
j_treat <- c()
k_treat <- c()
l_treat <- c()


for(i in sampleID){
    res3 <- res2 %>% dplyr::filter(sampleID==i)
    for (j in unique(res3$MCls)){
        res4 <- res3 %>% dplyr::filter(MCls==j)
        treats <- unique(res4$treat)
        #for (k in unique(res4$treat)){
        for (k in 1:length(treats)){
            #print(k)
            res5 <- res4 %>% dplyr::filter(treat==treats[k])
            #print(nrow(res5))
            if(k != length(treats)){
                for (l in (k+1):length(treats)){
                    #print(l)
                    res6 <- res4 %>% dplyr::filter(treat==treats[l])
                    ##res5 <- res4 %>% dplyr::filter(treat==k)
                    ##for (l in unique(res4$treat)){
                    ##    res6 <- res4 %>% dplyr::filter(treat==l)
                    scorr <- as.numeric(cor.test(res5$value, res6$value, method = "spearman")$estimate)
                    pval <- as.numeric(cor.test(res5$value, res6$value, method = "spearman")$p.value)
                    print(scorr)
                    scorr_treat <- c(scorr_treat, scorr)
                    pval_treat <- c(pval_treat, pval)
                    i_treat <- c(i_treat, i)
                    j_treat <- c(j_treat, j)
                    k_treat <- c(k_treat, treats[k])
                    l_treat <- c(l_treat, treats[l])
                }
            }
        }
    }
}

head(scorr_treat)
head(pval_treat)
head(i_treat)
head(j_treat)
head(k_treat)
head(l_treat)


df.treat <- data.frame(sampleID=i_treat, celltype=j_treat, treat1=k_treat, treat2=l_treat, corr=scorr_treat, pval=pval_treat)
head(df.treat)

nrow(df.treat %>% dplyr::filter(pval>0.05))





##--sampleID--##
MCls <- unique(res2$MCls)
MCls
#treat <- unique(res2$treat)
#sampleID <- unique(res2$sampleID)


scorr_sampleID <- c()
pval_sampleID <- c()
i_sampleID <- c()
j_sampleID <- c()
k_sampleID <- c()
l_sampleID <- c()

#head(res2)

for(i in MCls){
    res3 <- res2 %>% dplyr::filter(MCls==i)
    for (j in unique(res3$treat)){
        res4 <- res3 %>% dplyr::filter(treat==j)
        sampleIDs <- unique(res4$sampleID)
        ##for (k in unique(res4$sampleID)){
        for (k in 1:length(sampleIDs)){
            #print(k)
            res5 <- res4 %>% dplyr::filter(sampleID==sampleIDs[k])
            #print(nrow(res5))
            if(k != length(sampleIDs)){
                for (l in (k+1):length(sampleIDs)){
                    #print(l)
                    res6 <- res4 %>% dplyr::filter(sampleID==sampleIDs[l])
                    ##res5 <- res4 %>% dplyr::filter(sampleID==k)
                    ##for (l in unique(res4$sampleID)){
                    ##res6 <- res4 %>% dplyr::filter(sampleID==l)
                    scorr <- as.numeric(cor.test(res5$value, res6$value, method = "spearman")$estimate)
                    pval <- as.numeric(cor.test(res5$value, res6$value, method = "spearman")$p.value)
                    print(scorr)
                    scorr_sampleID <- c(scorr_sampleID, scorr)
                    pval_sampleID <- c(pval_sampleID, pval)
                    i_sampleID <- c(i_sampleID, i)
                    j_sampleID <- c(j_sampleID, j)
                    k_sampleID <- c(k_sampleID, sampleIDs[k])
                    l_sampleID <- c(l_sampleID, sampleIDs[l])
                }
            }
        }
    }
}

head(scorr_sampleID)
head(pval_sampleID)
head(i_sampleID)
head(j_sampleID)
head(k_sampleID)
head(l_sampleID)


df.sampleID <- data.frame(sampleID=i_sampleID, celltype=j_sampleID, sampleID1=k_sampleID, sampleID2=l_sampleID, corr=scorr_sampleID, pval=pval_sampleID)
head(df.sampleID)
nrow(df.sampleID %>% dplyr::filter(pval>0.05))



#----------------plot-----------------------
#df <- data.frame(scorr_MCls, scorr_treat)

#p <- boxplot(scorr_MCls, scorr_treat)
#png(paste0("./1_DiffRNA_2_3_corr.outs/20230117_corr_paper_no_repetation.png"), width=500, height=500, pointsize=16, res=125)
##boxplot(scorr_MCls, scorr_treat, scorr_sampleID)
png(paste0("./1_DiffRNA_2_3_corr.outs/20230117_corr_paper_no_repetation_ATAC.png"), width=600, height=600, res=125)
boxplot(scorr_sampleID, scorr_treat,
        main = "DEseq_ATAC_output_correlation",
        names=c("individuals", "treatments"))
dev.off()

t.test(scorr_sampleID, scorr_treat)

end






# ###############################################################################
# ### #----comparisons all celltypes as defined in paper - only treatments----
# ###############################################################################
##--treat - caff, nic, zinc vs water--##
#MCls <- unique(res2$MCls)
#treat <- unique(res2$treat)
sampleID <- unique(res2$sampleID)
sampleID

#unique(res2$treat)

treat_vswater <- c("caffeine", "nicotine", "zinc", "water")


scorr_treat_vswater <- c()
pval_treat_vswater <- c()
i_treat_vswater <- c()
j_treat_vswater <- c()
k_treat_vswater <- c()
l_treat_vswater <- c()


for(i in sampleID){
    res3 <- res2 %>% dplyr::filter(sampleID==i)
    for (j in unique(res3$MCls)){
        res4 <- res3 %>% dplyr::filter(MCls==j)
        treats <- treat_vswater
        for (k in 1:length(treats)){
            res5 <- res4 %>% dplyr::filter(treat==treats[k])
            if(k != length(treats)){
                for (l in (k+1):length(treats)){
                    res6 <- res4 %>% dplyr::filter(treat==treats[l])
                    scorr <- as.numeric(cor.test(res5$value, res6$value, method = "spearman")$estimate)
                pval <- as.numeric(cor.test(res5$value, res6$value, method = "spearman")$p.value)
                print(scorr)
                scorr_treat_vswater <- c(scorr_treat_vswater, scorr)
                pval_treat_vswater <- c(pval_treat_vswater, pval)
                i_treat_vswater <- c(i_treat_vswater, i)
                j_treat_vswater <- c(j_treat_vswater, j)
                k_treat_vswater <- c(k_treat_vswater, treats[k])
                l_treat_vswater <- c(l_treat_vswater, treats[l])
                }
            }
        }
    }
}

head(scorr_treat_vswater)
head(pval_treat_vswater)
head(i_treat_vswater)
head(j_treat_vswater)
head(k_treat_vswater)
head(l_treat_vswater)


df.treat_vswater <- data.frame(sampleID=i_treat_vswater, celltype=j_treat_vswater, treat_vswater1=k_treat_vswater, treat_vswater2=l_treat_vswater, corr=scorr_treat_vswater, pval=pval_treat_vswater)
head(df.treat_vswater)
nrow(df.treat_vswater %>% dplyr::filter(pval>0.05))







##--treat - vits vs etOH--##
#MCls <- unique(res2$MCls)
#treat <- unique(res2$treat)
sampleID <- unique(res2$sampleID)
sampleID

#unique(res2$treat)

treat_vsetOH <- c("vitA", "vitD", "vitE", "etOH")

scorr_treat_vsetOH <- c()
pval_treat_vsetOH <- c()
i_treat_vsetOH <- c()
j_treat_vsetOH <- c()
k_treat_vsetOH <- c()
l_treat_vsetOH <- c()


for(i in sampleID){
    res3 <- res2 %>% dplyr::filter(sampleID==i)
    for (j in unique(res3$MCls)){
        res4 <- res3 %>% dplyr::filter(MCls==j)
        treats <- treat_vsetOH
        for (k in 1:length(treats)){
            res5 <- res4 %>% dplyr::filter(treat==treats[k])
            if(k != length(treats)){
                for (l in (k+1):length(treats)){
                    res6 <- res4 %>% dplyr::filter(treat==treats[l])  
                    scorr <- as.numeric(cor.test(res5$value, res6$value, method = "spearman")$estimate)
                    pval <- as.numeric(cor.test(res5$value, res6$value, method = "spearman")$p.value)
                    print(scorr)
                    scorr_treat_vsetOH <- c(scorr_treat_vsetOH, scorr)
                    pval_treat_vsetOH <- c(pval_treat_vsetOH, pval)
                    i_treat_vsetOH <- c(i_treat_vsetOH, i)
                    j_treat_vsetOH <- c(j_treat_vsetOH, j)
                    k_treat_vsetOH <- c(k_treat_vsetOH, treats[k])
                    l_treat_vsetOH <- c(l_treat_vsetOH, treats[l])
                }
            }
        }
    }
}

head(scorr_treat_vsetOH)
head(pval_treat_vsetOH)
head(i_treat_vsetOH)
head(j_treat_vsetOH)
head(k_treat_vsetOH)
head(l_treat_vsetOH)


df.treat_vsetOH <- data.frame(sampleID=i_treat_vsetOH, celltype=j_treat_vsetOH, treat_vsetOH1=k_treat_vsetOH, treat_vsetOH2=l_treat_vsetOH, corr=scorr_treat_vsetOH, pval=pval_treat_vsetOH)
head(df.treat_vsetOH)
nrow(df.treat_vsetOH %>% dplyr::filter(pval>0.05))







##--treat - caff,zinc, nic vs vits--##
#MCls <- unique(res2$MCls)
#treat <- unique(res2$treat)
sampleID <- unique(res2$sampleID)
sampleID

#unique(res2$treat)

treat_1 <- c("caffeine", "zinc", "nicotine")
treat_2 <- c("vitA", "vitD", "vitE")


scorr_treat_wsvsvits <- c()
pval_treat_wsvsvits <- c()
i_treat_wsvsvits <- c()
j_treat_wsvsvits <- c()
k_treat_wsvsvits <- c()
l_treat_wsvsvits <- c()


for(i in sampleID){
    res3 <- res2 %>% dplyr::filter(sampleID==i)
    for (j in unique(res3$MCls)){
        res4 <- res3 %>% dplyr::filter(MCls==j)
        for (k in treat_1){
            res5 <- res4 %>% dplyr::filter(treat==k)
            for (l in treat_2){
                res6 <- res4 %>% dplyr::filter(treat==l)  
                scorr <- as.numeric(cor.test(res5$value, res6$value, method = "spearman")$estimate)
                pval <- as.numeric(cor.test(res5$value, res6$value, method = "spearman")$p.value)
                print(scorr)
                scorr_treat_wsvsvits <- c(scorr_treat_wsvsvits, scorr)
                pval_treat_wsvsvits <- c(pval_treat_wsvsvits, pval)
                i_treat_wsvsvits <- c(i_treat_wsvsvits, i)
                j_treat_wsvsvits <- c(j_treat_wsvsvits, j)
                k_treat_wsvsvits <- c(k_treat_wsvsvits, k)
                l_treat_wsvsvits <- c(l_treat_wsvsvits, l)
            }
        }
    }
}


head(scorr_treat_wsvsvits)
head(pval_treat_wsvsvits)
head(i_treat_wsvsvits)
head(j_treat_wsvsvits)
head(k_treat_wsvsvits)
head(l_treat_wsvsvits)


df.treat_wsvsvits <- data.frame(sampleID=i_treat_wsvsvits, celltype=j_treat_wsvsvits, treat_wsvsvits1=k_treat_wsvsvits, treat_wsvsvits2=l_treat_wsvsvits, corr=scorr_treat_wsvsvits, pval=pval_treat_wsvsvits)
head(df.treat_wsvsvits)
nrow(df.treat_wsvsvits %>% dplyr::filter(pval>0.05))








##--treat - water vs etOH--##
#MCls <- unique(res2$MCls)
#treat_cont <- unique(res2$treat)
sampleID <- unique(res2$sampleID)
sampleID

#unique(res2$treat)

treat_cont <- c("water")

scorr_treat_cont <- c()
pval_treat_cont <- c()
i_treat_cont <- c()
j_treat_cont <- c()
k_treat_cont <- c()
l_treat_cont <- c()


for(i in sampleID){
    res3 <- res2 %>% dplyr::filter(sampleID==i)
    for (j in unique(res3$MCls)){
        res4 <- res3 %>% dplyr::filter(MCls==j)
        for (k in treat_cont){
            res5 <- res4 %>% dplyr::filter(treat==k)
            #for (l in unique(res4$treat_cont)){
            res6 <- res4 %>% dplyr::filter(treat=="etOH")
            scorr <- as.numeric(cor.test(res5$value, res6$value, method = "spearman")$estimate)
            pval <- as.numeric(cor.test(res5$value, res6$value, method = "spearman")$p.value)
            print(scorr)
            scorr_treat_cont <- c(scorr_treat_cont, scorr)
            pval_treat_cont <- c(pval_treat_cont, pval)
            i_treat_cont <- c(i_treat_cont, i)
            j_treat_cont <- c(j_treat_cont, j)
            k_treat_cont <- c(k_treat_cont, k)
            l_treat_cont <- c(l_treat_cont, "etOH")
        }
    }
}
#}


head(scorr_treat_cont)
head(pval_treat_cont)
head(i_treat_cont)
head(j_treat_cont)
head(k_treat_cont)
head(l_treat_cont)


df.treat_cont <- data.frame(sampleID=i_treat_cont, celltype=j_treat_cont, treat_cont1=k_treat_cont, treat_cont2=l_treat_cont, corr=scorr_treat_cont, pval=pval_treat_cont)
head(df.treat_cont)
nrow(df.treat_cont %>% dplyr::filter(pval>0.05))






#----------------plot-----------------------
library("tidyverse")

#df <- data.frame(scorr_MCls, scorr_treat)

#p <- boxplot(scorr_MCls, scorr_treat)
png(paste0("./1_DiffRNA_2_3_corr.outs/corr_paper_all.png"), width=1500, height=1500, pointsize=16, res=225)
label=c("celltypes","individuals","treatments", "water_soluble_treatments", "etOH_soluble_treatments", "ws_vs_etOHs ", "water")
par(cex.axis=0.5) # is for x-axis
boxplot(scorr_MCls, scorr_sampleID, scorr_treat, scorr_treat_vswater, scorr_treat_vsetOH, scorr_treat_wsvsvits, scorr_treat_cont, names=label)
dev.off()



end

## # #################################################################
## # ### #----comparisons all celltypes----
## # #################################################################
## head(cvt)

## ##
## MCls <- unique(cvt$MCls)
## MCls

## limit <- length(MCls)
## limit

## cvt.filtered <- map_dfr(treat.IDs, function(oneX){
##     time0 <- Sys.time()
##     ##
##     cvt.2 <- cvt%>%dplyr::filter(treat==str_split(oneX, "\\.", simplify=T)[1], sampleID==str_split(oneX, "\\.", simplify=T)[2])
##     ##
##     cvt.2 <- cvt.2%>%mutate(count=ifelse(nrow(cvt.2)>=limit,1,0))
##     ##
##     time1 <- Sys.time()
##     elapsed <- difftime(time1, time0, units="mins")
##     cat(oneX, elapsed, "Done\n")
##     cvt.2
## })
## #cvt.filtered

## cvt.2 <- cvt.filtered %>% dplyr::filter(count==1)
## head(cvt.2)

## cvt.3 <- cvt.2 %>% mutate(vehicle=treat)
## #cvt.3 <- cvt.3 %>% mutate(vehicle=gsub("etOH", "etOH", vehicle))
## #cvt.3 <- cvt.3 %>% mutate(vehicle=gsub("water", "etOH", vehicle))
## cvt.3 <- cvt.3 %>% mutate(vehicle=gsub("caffeine", "water", vehicle))
## cvt.3 <- cvt.3 %>% mutate(vehicle=gsub("nicotine", "water", vehicle))
## cvt.3 <- cvt.3 %>% mutate(vehicle=gsub("zinc", "water", vehicle))
## cvt.3 <- cvt.3 %>% mutate(vehicle=gsub("vitA", "etOH", vehicle))
## cvt.3 <- cvt.3 %>% mutate(vehicle=gsub("vitD", "etOH", vehicle))
## cvt.3 <- cvt.3 %>% mutate(vehicle=gsub("vitE", "etOH", vehicle))

## cvt.3 <- cvt.3 %>% mutate(new_treat=treat)
## cvt.3 <- cvt.3 %>% mutate(new_treat=gsub("etOH", "control", new_treat))
## cvt.3 <- cvt.3 %>% mutate(new_treat=gsub("water", "control", new_treat))

## head(cvt.3)

## comb.2 <- table(cvt.2$MCls, cvt.2$treat)
## comb.2

## ##write.csv(comb.2, paste0("./1_DiffRNA_2.outs/combinations_2_", reso, "_cn.csv"), quote = FALSE, row.names = TRUE)
## ##comb.2 <- read.csv("./1_DiffRNA_2.outs/combinations_2.csv")
## ##comb.2

## bti3 <- cvt.2$bti
## bti3

## ##ncol(YtX)
## YtX.MCls <- YtX[,bti3]
## ##ncol(YtX)


## ##
## MCls <- unique(cvt.3$MCls)
## MCls

## scorr_MCls <- c()
## i_MCls <- c()
## j_MCls <- c()

## for(i in MCls){
##         for(j in MCls){
##                     cvt0 <- cvt.3%>%dplyr::filter(MCls==i)
##                             ##print(head(cvt0))
##                             ##print("###############new line###############")
##                             ##print(table(cvt0$vehicle, cvt0$treat))
##                             cvt1 <- cvt.3%>%dplyr::filter(MCls==j)
##                             #print(head(cvt1))
##                             YtX0 <- YtX.MCls[,cvt0$bti]
##                             YtX1 <- YtX.MCls[,cvt1$bti]
##                             dds0 <- DESeqDataSetFromMatrix(YtX0, cvt0, ~ treat)
##                             dds1 <- DESeqDataSetFromMatrix(YtX1, cvt1, ~ treat)
##                             ##dds1 <- DESeqDataSetFromMatrix(YtX1, cvt1, ~ sampleID + vehicle + treat)
##                             vd0 <- vst(dds0)
##                             #print(head(assay(vd0)))
##                             data0 <- melt(assay(vd0))
##                             ##print(head(data0))
##                             ##print(nrow(data0))
##                             #print(table(data0$Var2))
##                             #print(colnames(vd0))
##                             #print(head(vd0))
##                             #data0 <- plotPCA(vd0, intgroup=c("new_treat"), returnData=TRUE)
##                             #print(head(data0))
##                             #data0 <- data0 %>% mutate(sampleID=str_split(data$name, "_", simplify=T)[,4])
##                             #head(data0)
##                             vd1 <- vst(dds1)
##                             #print(head(vd1))
##                             data1 <- melt(assay(vd1))
##                             ##print(head(data1))
##                             ##print(nrow(data1))
##                             ###
##                             ###
##                             ###
##                             scorr <- as.numeric(cor.test(data0$value, data1$value, method = "spearman")$estimate)
##                             print(scorr)
##                             scorr_MCls <- c(scorr_MCls, scorr)
##                             i_MCls <- c(i_MCls, i)
##                             j_MCls <- c(j_MCls, j)
##                             ###
##                             #data1 <- plotPCA(vd1, intgroup=c("new_treat"), returnData=TRUE)
##                             #print(head(data1))
##                             #data1 <- data1 %>% mutate(sampleID=str_split(data$name, "_", simplify=T)[,4])
##                             #head(data1)
##                             }
##             }

## ##

## head(scorr_MCls)
## head(i_MCls)
## head(j_MCls)

## df <- data.frame(i_MCls, j_MCls, scorr_MCls)

## head(df)









## # #################################################################
## # ### #----comparisons across all treats----
## # #################################################################
## head(cvt)

## ##
## treat <- unique(cvt$treat)
## treat

## limit <- length(treat)
## limit

## cvt.filtered <- map_dfr(MCls.IDs, function(oneX){
##           time0 <- Sys.time()
##                   ###
##                   cvt.2 <- cvt%>%dplyr::filter(MCls==str_split(oneX, "\\.", simplify=T)[1], sampleID==str_split(oneX, "\\.", simplify=T)[2])
##                   ##
##                   cvt.2 <- cvt.2%>%mutate(count=ifelse(nrow(cvt.2)>=limit,1,0))
##                   ##
##                   time1 <- Sys.time()
##                   elapsed <- difftime(time1, time0, units="mins")
##                   cat(oneX, elapsed, "Done\n")
##                   cvt.2
##                 })
## #cvt.filtered

## cvt.2 <- cvt.filtered %>% dplyr::filter(count==1)
## head(cvt.2)

## cvt.3 <- cvt.2 %>% mutate(vehicle=treat)
## #cvt.3 <- cvt.3 %>% mutate(vehicle=gsub("etOH", "etOH", vehicle))
## #cvt.3 <- cvt.3 %>% mutate(vehicle=gsub("water", "etOH", vehicle))
## cvt.3 <- cvt.3 %>% mutate(vehicle=gsub("caffeine", "water", vehicle))
## cvt.3 <- cvt.3 %>% mutate(vehicle=gsub("nicotine", "water", vehicle))
## cvt.3 <- cvt.3 %>% mutate(vehicle=gsub("zinc", "water", vehicle))
## cvt.3 <- cvt.3 %>% mutate(vehicle=gsub("vitA", "etOH", vehicle))
## cvt.3 <- cvt.3 %>% mutate(vehicle=gsub("vitD", "etOH", vehicle))
## cvt.3 <- cvt.3 %>% mutate(vehicle=gsub("vitE", "etOH", vehicle))

## cvt.3 <- cvt.3 %>% mutate(new_treat=treat)
## cvt.3 <- cvt.3 %>% mutate(new_treat=gsub("etOH", "control", new_treat))
## cvt.3 <- cvt.3 %>% mutate(new_treat=gsub("water", "control", new_treat))

## head(cvt.3)

## comb.2 <- table(cvt.2$MCls, cvt.2$treat)
## comb.2

## ##write.csv(comb.2, paste0("./1_DiffRNA_2.outs/combinations_2_", reso, "_cn.csv"), quote = FALSE, row.names = TRUE)
## ##comb.2 <- read.csv("./1_DiffRNA_2.outs/combinations_2.csv")
## ##comb.2

## bti3 <- cvt.2$bti
## bti3

## ##ncol(YtX)
## YtX.treat <- YtX[,bti3]
## ##ncol(YtX)



## ##
## treats <- unique(cvt.3$treat)
## treats

## scorr_treats <- c()
## i_treats <- c()
## j_treats <- c()

## for(i in treats){
##         for(j in treats){
##                     cvt0 <- cvt.3%>%dplyr::filter(treat==i)
##                             cvt1 <- cvt.3%>%dplyr::filter(treat==j)
##                             YtX0 <- YtX[,cvt0$bti]
##                             YtX1 <- YtX[,cvt1$bti]
##                             dds0 <- DESeqDataSetFromMatrix(YtX0, cvt0, ~ sampleID + vehicle + new_treat)
##                             dds1 <- DESeqDataSetFromMatrix(YtX1, cvt1, ~ sampleID + vehicle + new_treat)
##                             vd0 <- vst(dds0)
##                             #print(head(assay(vd0)))
##                             data0 <- melt(assay(vd0))
##                             #print(head(data0))
##                             #print(nrow(data0))
##                             #print(table(data0$Var2))
##                             #print(colnames(vd0))
##                             #print(head(vd0))
##                             #data0 <- plotPCA(vd0, intgroup=c("new_treat"), returnData=TRUE)
##                             #print(head(data0))
##                             #data0 <- data0 %>% mutate(sampleID=str_split(data$name, "_", simplify=T)[,4])
##                             #head(data0)
##                             vd1 <- vst(dds1)
##                             #print(head(vd1))
##                             data1 <- melt(assay(vd1))
##                             #print(head(data1))
##                             #print(nrow(data1))
##                             ###
##                             ###
##                             ###
##                             scorr <- as.numeric(cor.test(data0$value, data1$value, method = "spearman")$estimate)
##                             print(scorr)
##                             scorr_treats <- c(scorr_treats, scorr)
##                             i_treats <- c(i_treats, i)
##                             j_treats <- c(j_treats, j)
##                             ###
##                             #data1 <- plotPCA(vd1, intgroup=c("new_treat"), returnData=TRUE)
##                             #print(head(data1))
##                             #data1 <- data1 %>% mutate(sampleID=str_split(data$name, "_", simplify=T)[,4])
##                             #head(data1)
##                             }
##             }

## ##

## head(scorr_treats)
## head(i_treats)
## head(j_treats)









## # #################################################################
## # ### #----comparisons across all individuals----
## # #################################################################
## head(cvt)

## ##
## MCls <- unique(cvt$MCls)
## MCls
## limit <- length(MCls)
## limit

## cvt.filtered <- map_dfr(treat.IDs, function(oneX){
##           time0 <- Sys.time()
##                   ###
##                   cvt.2 <- cvt%>%dplyr::filter(treat==str_split(oneX, "\\.", simplify=T)[1], sampleID==str_split(oneX, "\\.", simplify=T)[2])
##                   ##
##                   cvt.2 <- cvt.2%>%mutate(count=ifelse(nrow(cvt.2)>=limit,1,0))
##                   ##
##                   time1 <- Sys.time()
##                   elapsed <- difftime(time1, time0, units="mins")
##                   cat(oneX, elapsed, "Done\n")
##                   cvt.2
##                 })
## #cvt.filtered

## cvt.2 <- cvt.filtered %>% dplyr::filter(count==1)
## head(cvt.2)

## cvt.3 <- cvt.2 %>% mutate(vehicle=treat)
## #cvt.3 <- cvt.3 %>% mutate(vehicle=gsub("etOH", "etOH", vehicle))
## #cvt.3 <- cvt.3 %>% mutate(vehicle=gsub("water", "etOH", vehicle))
## cvt.3 <- cvt.3 %>% mutate(vehicle=gsub("caffeine", "water", vehicle))
## cvt.3 <- cvt.3 %>% mutate(vehicle=gsub("nicotine", "water", vehicle))
## cvt.3 <- cvt.3 %>% mutate(vehicle=gsub("zinc", "water", vehicle))
## cvt.3 <- cvt.3 %>% mutate(vehicle=gsub("vitA", "etOH", vehicle))
## cvt.3 <- cvt.3 %>% mutate(vehicle=gsub("vitD", "etOH", vehicle))
## cvt.3 <- cvt.3 %>% mutate(vehicle=gsub("vitE", "etOH", vehicle))

## cvt.3 <- cvt.3 %>% mutate(new_treat=treat)
## cvt.3 <- cvt.3 %>% mutate(new_treat=gsub("etOH", "control", new_treat))
## cvt.3 <- cvt.3 %>% mutate(new_treat=gsub("water", "control", new_treat))

## head(cvt.3)

## comb.2 <- table(cvt.2$MCls, cvt.2$treat)
## comb.2

## ##write.csv(comb.2, paste0("./1_DiffRNA_2.outs/combinations_2_", reso, "_cn.csv"), quote = FALSE, row.names = TRUE)
## ##comb.2 <- read.csv("./1_DiffRNA_2.outs/combinations_2.csv")
## ##comb.2

## bti3 <- cvt.2$bti
## bti3

## ##ncol(YtX)
## YtX <- YtX[,bti3]
## ##ncol(YtX)



## ##
## treats <- unique(cvt.3$treat)
## treats

## scorr_treats <- c()
## i_treats <- c()
## j_treats <- c()

## for(i in treats){
##         for(j in treats){
##                     cvt0 <- cvt.3%>%dplyr::filter(treat==i)
##                             cvt1 <- cvt.3%>%dplyr::filter(treat==j)
##                             YtX0 <- YtX[,cvt0$bti]
##                             YtX1 <- YtX[,cvt1$bti]
##                             dds0 <- DESeqDataSetFromMatrix(YtX0, cvt0, ~ sampleID + vehicle + new_treat)
##                             dds1 <- DESeqDataSetFromMatrix(YtX1, cvt1, ~ sampleID + vehicle + new_treat)
##                             vd0 <- vst(dds0)
##                             #print(head(assay(vd0)))
##                             data0 <- melt(assay(vd0))
##                             print(head(data0))
##                             print(nrow(data0))
##                             #print(table(data0$Var2))
##                             #print(colnames(vd0))
##                             #print(head(vd0))
##                             #data0 <- plotPCA(vd0, intgroup=c("new_treat"), returnData=TRUE)
##                             #print(head(data0))
##                             #data0 <- data0 %>% mutate(sampleID=str_split(data$name, "_", simplify=T)[,4])
##                             #head(data0)
##                             vd1 <- vst(dds1)
##                             #print(head(vd1))
##                             data1 <- melt(assay(vd1))
##                             print(head(data1))
##                             print(nrow(data1))
##                             ###
##                             ###
##                             #scorr <- as.numeric(cor.test(data0$value, data1$value, method = "spearman")$estimate)
##                             #print(scorr)
##                             #scorr_treats <- c(scorr_treats, scorr)
##                             #i_treats <- c(i_treats, i)
##                             #j_treats <- c(j_treats, j)
##                             ###
##                             #data1 <- plotPCA(vd1, intgroup=c("new_treat"), returnData=TRUE)
##                             #print(head(data1))
##                             #data1 <- data1 %>% mutate(sampleID=str_split(data$name, "_", simplify=T)[,4])
##                             #head(data1)
##                             }
##             }

## ##

## head(scorr_treats)
## head(i_treats)
## head(j_treats)











