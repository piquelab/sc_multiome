##
library(tidyverse)
library(data.table)
library(qvalue)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(RColorBrewer)
library(viridis)
library(openxlsx)

rm(list=ls())
 
outdir <- "./3.0_compare.outs/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


####
fn <- "../../torus_input/zzz_torus.annot.gz"
anno <- fread(fn, header=T, data.table=F)

x <- gsub("_d$", "", colnames(anno)[-c(1:2)])
x2 <- read.table("annot_ls.txt")$V1
comb <- unique(c(x,x2))


cvt <- str_split(comb, "_", simplify=T)
cvt2 <-  data.frame("comb"=comb, MCls=paste(cvt[,1], cvt[,2], sep="_"), treats=cvt[,3])


### current analysis
est2 <- map_dfr(comb, function(ii){
   ##
   fn <- paste("./torus_output/", ii, "_WBL.est", sep="")
   if ( file.exists(fn)){
      ### 
      x <- read.table(fn)
      ##kk <- ifelse(grepl("peaking", ii), 2, 3)
      kk <- 3
      x2 <- x[kk,]
      x2$comb <- ii
   }else{
      x2 <- NULL
   }   
   x2
})
 
est2 <- est2%>%dplyr::rename("odds"="V2", "CI_lower"="V3", "CI_upper"="V4")
plotDF <- cvt2%>%full_join(est2, by="comb")%>%arrange(MCls)


###
### color cell type
col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
  "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
   "6_Monocyte"="#984ea3", "7_dnT"="black") ##, "88_Peak"="#7570b3", "88_Response"="#e7298a")

## col_lab <- c("0_CD4Naive"="0_CD4Naive", "1_TCM"="1_TCM", "2_NKcell"="2_NKcell",
##   "3_TEM"="3_TEM", "4_Bcell"="4_Bcell", "5_CD8Naive"="5_CD8Naive",
##    "6_Monocyte"="6_Monocyte", "7_dnT"="7_dnT", "88_Peak"="Peak", "88_Response"="Response")

###
### color treats
## col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4", "vitD"="seagreen4",
##           "vitE"="salmon3", "zinc"="maroon3", "peaking"="#7570b3")
 
p1_new <- ggplot(plotDF, aes(x=odds, y=comb, color=factor(MCls)))+
    geom_errorbarh(aes(xmax=CI_upper, xmin=CI_lower), size=0.5, height=0.2)+
    geom_point(shape=19, size=0.5)+
    geom_vline(aes(xintercept=0), size=0.25, linetype="dashed")+
    ##scale_y_discrete(labels=ylab)+
    scale_x_continuous("log odds ratio", breaks=seq(-2, 3, by=1), limits=c(-2, 3))+
    scale_colour_manual(values=col1)+
    ggtitle("GTEx WBL")+
    theme_bw()+
    theme(legend.position="none",
          ##legend.key=element_blank(),
          legend.key.size=grid::unit(0.8, "lines"),
          legend.text=element_text(size=10),
          plot.title=element_text(hjust=0.5, size=12),
          axis.title.x=element_text(size=12),
          axis.title.y=element_blank(),
          axis.text.y=element_text(size=7))
          ##axis.title.x=element_text(size=10),
          ###strip.text=element_text(size=12),
          ##plot.title=element_text(hjust=0.5, size=12),
          ## axis.text.x=element_text(size=9, angle=-90, hjust=0, vjust=0.5),
          ###axis.text=element_text(size=9))


####
#### old analysis
est2 <- map_dfr(comb, function(ii){
   ##
   fn <- paste("../../analysis3_each/torus_output/", ii, "_WBL.est", sep="")
   if ( file.exists(fn)){
      ### 
      x <- read.table(fn)
      ##kk <- ifelse(grepl("peaking", ii), 2, 3)
      kk <- 3
      x2 <- x[kk,]
      x2$comb <- ii
   }else{
      x2 <- NULL
   }   
   x2
})
 
est2 <- est2%>%dplyr::rename("odds"="V2", "CI_lower"="V3", "CI_upper"="V4")
plotDF2 <- cvt2%>%full_join(est2, by="comb")%>%arrange(MCls)
 

p2_old <- ggplot(plotDF2, aes(x=odds, y=comb, color=factor(MCls)))+
    geom_errorbarh(aes(xmax=CI_upper, xmin=CI_lower), size=0.5, height=0.2)+
    geom_point(shape=19, size=0.5)+
    geom_vline(aes(xintercept=0), size=0.25, linetype="dashed")+
    ##scale_y_discrete(labels=ylab)+
    scale_x_continuous("log odds ratio", breaks=seq(-2, 3, by=1), limits=c(-2, 3))+
    scale_colour_manual(values=col1)+
    ggtitle("GTEx WBL(Old)")+
    theme_bw()+
    theme(legend.title=element_blank(),
          ##legend.key=element_blank(),
          legend.key.size=grid::unit(0.8, "lines"),
          legend.text=element_text(size=9),
          plot.title=element_text(hjust=0.5, size=12),
          axis.title.x=element_text(size=12),
          axis.title.y=element_blank(),
          axis.text.y=element_text(size=7))


###
###
pcomb <- plot_grid(p1_new, p2_old, rel_widths=c(1, 1.4))
figfn <- paste(outdir, "Figure1_new_old_forest.png", sep="")
ggsave(figfn, pcomb, device="png", width=1000, height=480, units="px", dpi=120)



################################################################################
#### compare number Response motif and number of genetic variants annotated ####
################################################################################

###
### response motif
fn <- "/nfs/rprdata/julong/sc_multiome/genetic_variant_annot/4_SNPAnnot_correct/2.2_response_motif_th0.1.txt"
res <- read.table(fn, header=T)

fn <- "/nfs/rprdata/julong/sc_multiome/genetic_variant_annot/Response_motif/2.2_response_motif_th0.1.txt"
res_old <- read.table(fn, header=T)%>%
    filter(!grepl("8_MAIT|water", comb))

comb <- sort(unique(c(res$comb, res_old$comb)))
x <- str_split(comb, "_", simplify=T)

cvt <- data.frame("comb"=comb, MCls=paste(x[,1], x[,2], sep="_"), treats=x[,3])


###
### annotation
fn  <- "../../eQTL_results/snpList.txt"
snpSel <- read.table(fn)$V1

fn <- "/nfs/rprdata/julong/sc_multiome/genetic_variant_annot/4_SNPAnnot_correct/zzz_th0.1_multiConditions_torus.annot.gz"
anno <- fread(fn, header=T, data.table=F)
anno <- anno%>%filter(SNP%in%snpSel)
nsnp <- colSums(anno[,-c(1,2)])
names(nsnp) <- gsub("_d$", "", names(nsnp))


fn <- "/nfs/rprdata/julong/sc_multiome/genetic_variant_annot/4_SNPAnnot.outs/zzz_th0.1_multiConditions_torus.annot.gz"
anno_old <- fread(fn, header=T, data.table=F)
anno_old <- anno_old%>%filter(SNP%in%snpSel)
nsnp_old <- colSums(anno_old[,-c(1,2)])
names(nsnp_old) <- gsub("_d$", "", names(nsnp_old))



###
### compare current and old
summ <- map_dfr(comb, function(ii){
   ###
   print(ii) 
   motifs <- res%>%filter(comb==ii)%>%pull(motif_ID)%>%unique()
   motifs_old <- res_old%>%filter(comb==ii)%>%pull(motif_ID)%>%unique()
   motifs_shared <- intersect(motifs, motifs_old)
   ###
   n_snps <- as.numeric(ifelse(is.na(nsnp[ii]), 0, nsnp[ii]))
   n_snps_old <- as.numeric(ifelse(is.na(nsnp_old[ii]), 0, nsnp_old[ii]))

    
   ### SNPs
   ii2 <- paste(ii, "_d", sep="") 
   snps <- NULL 
   if ( n_snps!=0){
       aa <- anno[,c("SNP", ii2)]
       snps <- aa[aa[,2]==1,]%>%pull(SNP)%>%unique()
   }
   ##
   snps_old <- NULL
   if ( n_snps_old!=0){
       aa <- anno_old[,c("SNP", ii2)]
       snps_old <- aa[aa[,2]==1,]%>%pull(SNP)%>%unique()
   }    
   snps_shared <- intersect(snps, snps_old)
    
   df2 <- data.frame(comb=ii, n_motifs=length(motifs), "n_snps"=n_snps,
      n_motifs_old=length(motifs_old), "n_snps_old"=n_snps_old,
      olap_motifs=length(motifs_shared), olap_snps=length(snps_shared))
   df2
})
    
###
opfn <- paste(outdir, "1_compare_SNPs_old.xlsx", sep="")
write.xlsx(summ, file=opfn, overwrite=T)


