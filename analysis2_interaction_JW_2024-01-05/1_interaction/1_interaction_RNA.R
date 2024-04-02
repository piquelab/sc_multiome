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

###
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

outdir <- "./1_inter_RNA.outs/"
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=T, recursive=T)


library("BiocParallel")
register(MulticoreParam(4))


###################################################
#### Prepare input data for DESeq2 
###################################################


###
### read data
YtX <- read_rds("../sc_multiome_data_MHB/2_Differential/1_DiffRNA_2.outs/1_YtX.sel_clusters_0.1_cn.rds")


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


#####################################################################################
### run DESeq2 full model against reduced model (w/0 interaction terms)
#####################################################################################

dds <- DESeqDataSetFromMatrix(
    countData=YtX2,
    colData=cvt2,
    design=~MCls+treat2+MCls:treat2+sampleID)
system.time(dds <- DESeq(dds, test="LRT", reduced=~MCls+treat2+sampleID, parallel=T))

### output 
opfn <- paste(outdir, "1_reduce_RNA.dds.rds", sep="")
write_rds(dds, file=opfn)

###
### results-lrt
res <- results(dds)%>%tidy() ## as.data.frame()%>%rownames_to_column(var="gene")
opfn <- paste(outdir, "1.2_lrt_results.rds", sep="")
write_rds(res, file=opfn)



#######################################################################
### run DESeq with ~new_gr+sampleID
#######################################################################

cvt2 <- cvt2%>%mutate(new_gr=paste(MCls, treat2, sep="_"))
dds <- DESeqDataSetFromMatrix(
    countData=YtX2,
    colData=cvt2,
    design=~new_gr+sampleID)
system.time(dds_full <- DESeq(dds, parallel=T))
## output
opfn2 <- paste(outdir, "2_full_RNA.dds.rds", sep="")
write_rds(dds_full, file=opfn2)


#########################################################
### Create  contrast list and extract results 
#########################################################

###
### contrast-1, extract results from dds

MCls <- sort(unique(cvt2$MCls))
treats <- sort(unique(cvt2$treat2))[-2]
##
n_MCls <- length(MCls)
n_trs <- length(treats)

contrast_df <- data.frame(contrast_nn=rep("new_gr", each=n_trs*n_MCls),
    con1=paste(rep(MCls, each=n_trs), rep(treats, times=n_MCls), sep="_"),
    con0 =paste(rep(MCls, each=n_trs), rep("control", each=n_trs*n_MCls), sep="_"))
###
opfn <- paste(outdir, "contrast_list.txt", sep="")
write.table(contrast_df, file=opfn, row.names=F, quote=F)


###
### contrast-2, compare the marginal effects across cell-type
MCls <- sort(unique(cvt2$MCls))[-1]
treats <- sort(unique(cvt2$treat2))[-2]
n_MCls <- length(MCls)
n_trs <- length(treats)

contrast_df2 <- data.frame(con1=paste(rep(MCls, times=n_trs), rep(treats, each=n_MCls), sep="_"),
   con0=paste(rep("0_CD4Naive", n_trs*n_MCls), rep(treats, each=n_MCls), sep="_"))
###
opfn <- paste(outdir, "contrast_list2.txt", sep="")
write.table(contrast_df2, file=opfn, row.names=F, quote=F)




###
### extract results using 1.1_submit.sh


###
### combine results

fns <- list.files("./1_inter_RNA.outs/dds_results/", ".*.rds")

res_comb <- map_dfr(fns, function(ii){
    ###
    cat(ii, "\n")
    fn <- paste("./1_inter_RNA.outs/dds_results/", ii, sep="")
    res0 <- read_rds(fn)
    res0
})    

opfn <- paste(outdir, "2.2_each_pair.results.rds", sep="")
write_rds(res_comb, file=opfn)


## ###
## system.time(res2 <- contrast.list%>%map(~results(dds_full, contrast=.x)))

## ## parallel
## library(parallel)
## num_mc <- detectCores()
 
## system.time(res2_ls <- mclapply(contrast.list, function(x){
##     ##
##     res0 <- results(dds_full, contrast=x)
##     res0
## }, mc.cores=10))    

## res2_comb <- tibble(treats=names(contrast.list), data=map(res2_ls, tidy))%>%unnest(data)
## opfn3 <- paste(outdir, "2.2_each_pair.results.rds", sep="")
## write_rds(res2_comb, file=opfn3)
 











