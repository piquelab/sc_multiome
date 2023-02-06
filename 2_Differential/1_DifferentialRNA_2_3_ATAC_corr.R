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



#################################################################################
#################################################################################
##    ATAC    ##
#################################################################################
#################################################################################
####################################
### sys_args ###
####################################
reso = 0.1

####################################
### 1. Generate pseudo-bulk data ###
####################################
YtX_sel <- read_rds("./1.2_DiffPeak.outs/1_YtX.sel_0.1_cn.rds")

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




# #################################################################
# ### 2.1 call DESeq ###
# #################################################################
res.DE.ATAC <- readRDS("./1_DiffPeak.outs/2.0_DESeq.results_clusters_treatandind_filtered_0.1_cn.rds")
head(res.DE.ATAC)
table(res.DE.ATAC$MCls)


######################################################################
##  resDEATAC - only treatments
######################################################################
##--treat_ - caff, nic, zinc vs water--##
head(res.DE.ATAC)

MCls <- unique(res.DE.ATAC$MCls)

treat_vswater <- c("caffeine", "nicotine", "zinc")

scorr_treat_vswater <- c()
pval_treat_vswater <- c()
k_treat_vswater <- c()
i_treat_vswater <- c()
#j_treat_vswater <- c()

for(k in MCls){
    res.DE.ATAC2 <- res.DE.ATAC%>%dplyr::filter(MCls==k)
    for(i in treat_vswater){
        res.DE.ATAC3 <- res.DE.ATAC2%>% dplyr::filter(contrast==i)
        ##for(j in treat_vswater){
        res.DE.ATAC4 <- res.DE.ATAC2%>%dplyr::filter(contrast=="water")
        ##
        scorr <- as.numeric(cor.test(res.DE.ATAC3$estimate, res.DE.ATAC4$estimate, method = "spearman")$estimate)
        pval <- as.numeric(cor.test(res.DE.ATAC3$estimate, res.DE.ATAC4$estimate, method = "spearman")$p.value)
        print(scorr)
        scorr_treat_vswater <- c(scorr_treat_vswater, scorr)
        pval_treat_vswater <- c(pval_treat_vswater, pval)
        k_treat_vswater <- c(k_treat_vswater, k)
        i_treat_vswater <- c(i_treat_vswater, i)
        ##            j_treat_vswater <- c(j_treat_vswater, j)
    }
}
#}

head(scorr_treat_vswater)
head(pval_treat_vswater)
head(k_treat_vswater)
head(i_treat_vswater)
#head(j_treat_vswater)


df.treat_vswater <- data.frame(celltype=k_treat_vswater, treat1=i_treat_vswater, treat2="water", corr=scorr_treat_vswater, pval=pval_treat_vswater)
head(df.treat_vswater)
nrow(df.treat_vswater %>% dplyr::filter(pval>0.05))







##--treat_ - vits vs etOH--##
table(res.DE.ATAC$MCls)
table(res.DE.ATAC$contrast)

MCls <- unique(res.DE.ATAC$MCls)

#treat_ <- unique(res.DE.ATAC$contrast)
#treat_

treat_vsetOH <- c("vitA", "vitD", "vitE")

scorr_treat_vsetOH <- c()
pval_treat_vsetOH <- c()
k_treat_vsetOH <- c()
i_treat_vsetOH <- c()
#j_treat_vsetOH <- c()

for(k in MCls){
        res.DE.ATAC2 <- res.DE.ATAC%>%dplyr::filter(MCls==k)
        for(i in treat_vsetOH){
            res.DE.ATAC3 <- res.DE.ATAC2%>% dplyr::filter(contrast==i)
            ##for(j in treat_vsetOH){
            res.DE.ATAC4 <- res.DE.ATAC2%>%dplyr::filter(contrast=="etOH")
            ##
            scorr <- as.numeric(cor.test(res.DE.ATAC3$estimate, res.DE.ATAC4$estimate, method = "spearman")$estimate)
            pval <- as.numeric(cor.test(res.DE.ATAC3$estimate, res.DE.ATAC4$estimate, method = "spearman")$p.value)
            print(scorr)
            scorr_treat_vsetOH <- c(scorr_treat_vsetOH, scorr)
            pval_treat_vsetOH <- c(pval_treat_vsetOH, pval)
            k_treat_vsetOH <- c(k_treat_vsetOH, k)
            i_treat_vsetOH <- c(i_treat_vsetOH, i)
            ##            j_treat_vsetOH <- c(j_treat_vsetOH, j)
        }
}
#}

head(scorr_treat_vsetOH)
head(pval_treat_vsetOH)
head(k_treat_vsetOH)
head(i_treat_vsetOH)
#head(j_treat_vsetOH)


df.treat_vsetOH <- data.frame(celltype=k_treat_vsetOH, treat1=i_treat_vsetOH, treat2="water", corr=scorr_treat_vsetOH, pval=pval_treat_vsetOH)
head(df.treat_vsetOH)
nrow(df.treat_vsetOH %>% dplyr::filter(pval>0.05))







##--treat_ - water vs etOH--##
MCls <- unique(res.DE.ATAC$MCls)
#treat_ <- unique(res.DE.ATAC$treat_)

treat_cont <- c("caffeine", "nicotine", "zinc")

scorr_treat_ <- c()
pval_treat_cont <- c()
k_treat_cont <- c()
i_treat_cont <- c()
#j_treat_cont <- c()

for(k in MCls){
        res.DE.ATAC2 <- res.DE.ATAC%>%dplyr::filter(MCls==k)
            for(i in treat_cont){
                res.DE.ATAC3 <- res.DE.ATAC2%>% dplyr::filter(contrast==i)
                ##for(j in treat_cont){
                res.DE.ATAC4 <- res.DE.ATAC2%>%dplyr::filter(contrast=="water")
                ##
                scorr <- as.numeric(cor.test(res.DE.ATAC3$estimate, res.DE.ATAC4$estimate, method = "spearman")$estimate)
                pval <- as.numeric(cor.test(res.DE.ATAC3$estimate, res.DE.ATAC4$estimate, method = "spearman")$p.value)
                print(scorr)
                scorr_treat_cont <- c(scorr_treat_cont, scorr)
                pval_treat_cont <- c(pval_treat_cont, pval)
                k_treat_cont <- c(k_treat_cont, k)
                i_treat_cont <- c(i_treat_cont, i)
                ##            j_treat_cont <- c(j_treat_cont, j)
            }
}
#}

head(scorr_treat_cont)
head(pval_treat_cont)
head(k_treat_cont)
head(i_treat_cont)
#head(j_treat_cont)


df.treat_cont <- data.frame(celltype=k_treat_cont, treat1=i_treat_cont, treat2="water", corr=scorr_treat_cont, pval=pval_treat_cont)
head(df.treat_cont)
nrow(df.treat_cont %>% dplyr::filter(pval>0.05))






#----------------plot-----------------------
#df <- data.frame(scorr_MCls, scorr_treat_)

#p <- boxplot(scorr_MCls, scorr_treat_)
png(paste0("./1_DiffRNA_2_3_corr.outs/corr_resDE_treat_ments_ATAC.png"), width=500, height=500, pointsize=16, res=125)
#boxplot(scorr_treat_vswater, scorr_treat_vsetOH, scorr_treat_cont)
boxplot(scorr_treat_vswater)
dev.off()



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


res <- map_dfr(MCls, function(oneX){
      cvt0 <- cvt.3%>%dplyr::filter(MCls==oneX)
        YtX0 <- YtX3[,cvt0$bti]
        dds0 <- DESeqDataSetFromMatrix(YtX0, cvt0, ~ sampleID + vehicle + new_treat)
        vd0 <- vst(dds0)
        data0 <- melt(assay(vd0))
      })
head(res)


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
png(paste0("./1_DiffRNA_2_3_corr.outs/corr_paper_no_repetation.png"), width=500, height=500, pointsize=16, res=125)
boxplot(scorr_MCls, scorr_treat, scorr_sampleID)
dev.off()











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
#df <- data.frame(scorr_MCls, scorr_treat)

#p <- boxplot(scorr_MCls, scorr_treat)
png(paste0("./1_DiffRNA_2_3_corr.outs/corr_paper_ATAC_all.png"), width=500, height=500, pointsize=16, res=125)
label=c("celltypes","individuals","treatments", "water_soluble_treatments", "etOH_soluble_treatments", "ws_vs_etOHs ", "water")
par(cex.axis=0.5) # is for x-axis
boxplot(scorr_MCls, scorr_sampleID, scorr_treat, scorr_treat_vswater, scorr_treat_vsetOH, scorr_treat_wsvsvits, scorr_treat_cont, names=label)
#boxplot(scorr_MCls, scorr_treat, scorr_sampleID, scorr_treat_vswater, scorr_treat_vsetOH, scorr_treat_cont)
dev.off()
