###
library(tidyverse)
##
library(cowplot) ##lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(grid)
library(gridExtra)
library(ggExtra)
library(RColorBrewer)
library(ggsci)
library(viridis)
library(ComplexHeatmap)
library(circlize)
library(openxlsx)


###
###

##########################
### By JW, Apr-05,2023 ###
##########################


#######################
### cell type peaks ###
#######################

outdir2 <- "./Peak_bed/"
if (!file.exists(outdir2)) dir.create(outdir2, showWarnings=F, recursive=T)

## combined <- read_rds("../sc_multiome_data/1_processing/3_Clustering.outs/1_seurat.cluster.combined.mitofilt.rds")

fn <- "../sc_multiome_data/3_motif/1_motif.outs/1_scATAC.motif.rds"
atac <- read_rds(fn)
meta <- atac@meta.data
##
meta2 <- meta%>%mutate(MCl2=paste(wsnn_res.0.1, MCls, sep="_"))
atac <- AddMetaData(atac, meta2)

MCls <- c("0_CD4Naive", "1_TCM", "2_NKcell", "3_TEM", "4_Bcell", "5_CD8Naive", "6_Monocyte", "7_dnT")

### threshold 0.01
for (oneMCl in MCls){
   ##
   atac2 <- subset(atac, MCl2==oneMCl)
   x <- atac2@assays$ATAC@counts
   rms <- rowMeans(x>0) 
   sub0 <- rms[rms>0.01]
   cat(oneMCl, length(sub0), "\n") 
   ##
   cvt <- str_split(names(sub0), "-", simplify=T)
   opfn <- paste(outdir2, oneMCl, "_peak.bed", sep="")
   write.table(cvt, file=opfn, quote=F, row.names=F, col.names=F, sep="\t")
}    
##    
    

###
### summary

library(openxlsx)
### threshold 0.01
summ <- NULL
###
for (oneMCl in MCls){
   ##
   atac2 <- subset(atac, MCl2==oneMCl)
   x <- atac2@assays$ATAC@counts
   summ2 <- map_dfr(c(0.005, 0.01, 0.02, 0.05, 0.08, 0.1), function(ii){
   ###    
      rms <- rowMeans(x>0) 
      cat(oneMCl, sum(rms>0), "\n")
      ###    
      tmp <- data.frame(MCls=oneMCl, th=ii, ntotal=sum(rms>0), npeak=sum(rms>ii))
      tmp 
   })                  
   summ <- rbind(summ, summ2)
}    
##    

summ <- summ%>%mutate(prop=npeak/ntotal)

tmp <- summ%>%
    pivot_wider(id_cols=c(MCls, ntotal),
                names_from=th, names_prefix="th_", values_from=npeak)
opfn <- "./Peak_bed/0_summary_peaks.xlsx"
write.xlsx(tmp, file=opfn, overwrite=T)




#######################
#### response motif ###
#######################


outdir2 <- "./Response_motif/"
if (!file.exists(outdir2)) dir.create(outdir2, showWarnings=F, recursive=T)

fn <- "../sc_multiome_data/3_motif/2_motif.activities.outs/3_motif.diff.results_macs2_0.1_cn.rds"
res <- read_rds(fn)
res <- res%>%mutate(comb=paste(MCls, contrast, sep="_"))
###
opfn1 <- paste(outdir2, "1_diff_motif.results.rds", sep="")
write_rds(res, file=opfn1)





###
### threshold  across conditions
th0 <- quantile(abs(res$beta), probs=0.9)
res2 <- res%>%dplyr::filter(qval<0.1, abs(beta)>th0)%>%dplyr::select(gene, motif, comb)
opfn2 <- paste(outdir2, "2_response_motif.txt", sep="")
write.table(res2, file=opfn2, quote=F, row.names=F, col.names=T, sep="\t")



###
### default each condition

### 10%
fn <- "./Response_motif/1_diff_motif.results.rds"
res <- read_rds(fn)
##
comb2 <- sort(unique(res$comb))
res2 <- map_dfr(comb2, function(ii){
   ##
   tmp <- res%>%dplyr::filter(comb==ii)
   th0 <- quantile(abs(tmp$beta), probs=0.9)
   tmp2 <- tmp%>%dplyr::filter(qval<0.1, abs(beta)>th0)%>%
       dplyr::select(motif_ID=gene, motif_name=motif, comb)
   cat(ii, nrow(tmp2), "\n") 
   tmp2
})    

opfn2 <- paste(outdir2, "2.2_response_motif_th0.1.txt", sep="")
write.table(res2, file=opfn2, quote=F, row.names=F, col.names=T, sep="\t")


### 20%
##
comb2 <- sort(unique(res$comb))
res2 <- map_dfr(comb2, function(ii){
   ##
   tmp <- res%>%dplyr::filter(comb==ii)
   th0 <- quantile(abs(tmp$beta), probs=0.8)
   tmp2 <- tmp%>%dplyr::filter(qval<0.1, abs(beta)>th0)%>%
       dplyr::select(motif_ID=gene, motif_name=motif, comb)
   cat(ii, nrow(tmp2), "\n") 
   tmp2
})    
 
opfn2 <- paste(outdir2, "2.3_response_motif_th0.2.txt", sep="")
write.table(res2, file=opfn2, quote=F, row.names=F, col.names=T, sep="\t")



###
### summary #response motif

summ <- map_dfr(comb2, function(ii){
   ##
   tmp <- res%>%dplyr::filter(comb==ii)
   th0 <- quantile(abs(tmp$beta), probs=0.9)
   tmp2 <- tmp%>%dplyr::filter(qval<0.1, abs(beta)>th0)%>%
       dplyr::select(motif_ID=gene, motif_name=motif, comb)
   df <- data.frame(comb=ii, th=th0, nmotif=nrow(tmp2))
   df 
})
names(summ) <- c("comb", "th_0.1", "nmotif_0.1")

###
### th 0.2
summ2 <- map_dfr(comb2, function(ii){
   ##
   tmp <- res%>%dplyr::filter(comb==ii)
   th0 <- quantile(abs(tmp$beta), probs=0.8)
   tmp2 <- tmp%>%dplyr::filter(qval<0.1, abs(beta)>th0)%>%
       dplyr::select(motif_ID=gene, motif_name=motif, comb)
   df <- data.frame(comb=ii, th=th0, nmotif=nrow(tmp2))
   df 
})
names(summ2) <- c("comb", "th_0.2", "nmotif_0.2")
summ <- summ%>%inner_join(summ2, by="comb")

###
cvt <- str_split(summ$comb, "_", simplify=T)
summ <- summ%>%
   mutate(MCls=paste(cvt[,1], cvt[,2], sep="_"), treats=cvt[,3])%>%
   arrange(treats)%>%filter(!grepl("8_MAIT|water", comb))

    

opfn <- paste(outdir2, "2_summary.xlsx", sep="")
write.xlsx(summ, file=opfn, overwrite=T)
