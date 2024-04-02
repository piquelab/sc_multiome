###
library(Matrix)
library(tidyverse)
library(qqman)
library(qvalue)
##
library(DESeq2)
library(biobroom)
library(ashr)
library(GenomicRanges)

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

outdir <- "./2_inter_ATAC.outs/"
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=T, recursive=T)



library("BiocParallel")
register(MulticoreParam(5)) 
 

##########################
# input data 
#########################

###
### read data
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


##################################################################################
### run DESeq2 with full model agains reduce model w/o interaction term
##################################################################################

dds <- DESeqDataSetFromMatrix(
    countData=YtX2,
    colData=cvt2,
    design=~MCls+treat2+MCls:treat2+sampleID)
 
system.time(dds <- DESeq(dds, test="LRT", reduced=~MCls+treat2+sampleID, parallel=T))
###
opfn <- paste(outdir, "1_reduce_ATAC.dds.rds", sep="")
write_rds(dds, file=opfn)


###
### results 
res <- results(dds)%>%tidy()
opfn <- paste(outdir, "1.2_lrt_results.rds", sep="")
write_rds(res, file=opfn)





###############################################
### run DESeq with full model
###############################################

cvt2 <- cvt2%>%mutate(new_gr=paste(MCls, treat2, sep="_"))
dds <- DESeqDataSetFromMatrix(
    countData=YtX2,
    colData=cvt2,
    design=~new_gr+sampleID)

system.time(dds_full <- DESeq(dds, parallel=T))

opfn2 <- paste(outdir, "2_full_ATAC.dds.rds", sep="")
write_rds(dds_full, file=opfn2)



#######################################################
### extract results for defined contrasts 
#######################################################


### bash 2.1_submit.sh 


###
### combine results

files <- list.files("./2_inter_ATAC.outs/dds_results/", ".*.rds")
###
res_comb <- map_dfr(files, function(ii){
    ##
    cat(ii, "\n")
    fn <- paste("./2_inter_ATAC.outs/dds_results/", ii, sep="")
    res0 <- read_rds(file=fn)
    res0
})

###
### output 
opfn <- paste(outdir, "2.2_each_pair.results.rds", sep="")
write_rds(res_comb, file=opfn)
    








