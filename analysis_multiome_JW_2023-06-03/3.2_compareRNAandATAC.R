###
library(Matrix)
library(tidyverse)
library(data.table)
library(irlba)
##
library(cowplot)
library(RColorBrewer)
library(ComplexHeatmap)
library(viridis)
library(circlize)
library(ggrepel)
library(ggrastr)
library(openxlsx)

library(MASS)
library("rainbow", lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages")
library("fds", lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages")
library(fda, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages")
library(CCA, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages")

library("GGally")

rm(list=ls())



outdir <- "./3_compare_plots.outs/"
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)

###
### clean script, in the final we used union of DARs to calculate LFC of gene chromatin accessibility 


#########################################################################
### Box plots of correlation between LFC on gene expression and peaks ###
#########################################################################


### 1.2_*, median union DARs ## used for final analysis


##########################################
### correlation of LFC on RNA and ATAC ###
##########################################

#### DEGs
fn <- "./sc_multiome_data/2_Differential/1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn_sampleID+new_treat.rds"
resDG <- read_rds(fn)%>%as.data.frame()%>%mutate(comb=paste(MCls, contrast, sep="_"))
DEG <- resDG%>%dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)%>%pull(gene)%>%unique()
###resDG <- resDG%>%filter(gene%in%DEG)


### DARs
fn <- "./sc_multiome_data/2_Differential/1.2_DiffPeak.outs/2.0_DESeq.results_0.1_cn_sampleID+new_treat.rds"
resDP <- read_rds(fn)%>%as.data.frame()%>%mutate(comb=paste(MCls, contrast, sep="_"))
resDP <- resDP%>%dplyr::rename("peak"="gene")
DP <- resDP%>%dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)%>%pull(peak)%>%unique()

###
fn <- "./sc_multiome_data/2_Differential/2.2_compareRNAandATAC.outs/1_annot.signac_macs2_0.1_cn.rds"
peakAnno <- read_rds(fn)
peakAnno <- peakAnno%>%dplyr::select(peak=query_region, gene_id=gene_name, dtss=distance)%>%
    mutate(peak=gsub(":", "-", peak), dtss=abs(dtss))

 

####################################################################
### (2) another option, only calcualte median value of LFC of DP ###
####################################################################
comb <- sort(unique(resDP$comb))
dfcomb <- map_dfr(comb, function(ii){
   ###
   cat(ii, "\n") 
   x <- resDP%>%filter(comb==ii)%>%dplyr::select(estimate, peak)
   x2 <- x%>%left_join(peakAnno, by="peak")
   x2 <- x2 %>%filter(dtss<1e+05, peak%in%DP) ## 100 kb
   ###
   x3 <- x2%>%group_by(gene_id)%>%summarize(LFC_ATAC=median(estimate, na.rm=T),.groups="drop")
   x3 <- x3%>%filter(gene_id%in%DEG) 
   ###
    
   df2 <- resDG%>%dplyr::filter(comb==ii, gene%in%DEG)%>%
       dplyr::select(comb, MCls, contrast, gene, LFC_gene=estimate)%>%
       inner_join(x3, by=c("gene"="gene_id"))
   ### 
   df2
})


opfn <- paste(outdir, "1.2_LFC_gene_ATAC_comb.rds", sep="")
write_rds(dfcomb, file=opfn)



####
#### box plot show correlation

fn <- paste(outdir, "1.2_LFC_gene_ATAC_comb.rds", sep="")
df <- read_rds(fn)%>%drop_na(LFC_gene, LFC_ATAC) 

dfcorr <- df%>%group_by(MCls, contrast)%>%
    summarize(rr=cor(LFC_gene, LFC_ATAC, method="spearman"), .groups="drop")%>%as.data.frame()

rmat <- dfcorr%>%pivot_wider(id_cols=MCls, names_from=contrast, values_from=rr)%>%as.data.frame()


###
### setting colors
col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
  "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
   "6_Monocyte"="#984ea3", "7_dnT"="black")
col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
       "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3")

### group by cell types
p1 <- ggplot(dfcorr, aes(x=MCls, y=rr, fill=MCls))+
   ##geom_violin()+ 
   geom_boxplot(outlier.shape=NA, lwd=0.25)+
   stat_summary(fun=median, color="grey", geom="point", shape=23, size=2, stroke=0.9)+
   geom_jitter(width=0.2, size=0.5)+     
   scale_fill_manual(values=col1)+
   ylab("Spearman's correlation")+ylim(-0.1, 0.5)+
   theme_bw()+
   theme(legend.position="none",
         axis.title.x=element_blank(),
         axis.title.y=element_text(size=9),
         axis.text.x=element_text(angle=45, hjust=1, size=9),
         axis.text.y=element_text(size=9))

###
### group by treat
p2 <- ggplot(dfcorr, aes(x=contrast, y=rr, fill=contrast))+
   ##geom_violin()+ 
   geom_boxplot(outlier.shape=NA, lwd=0.25)+
   stat_summary(fun=median, color="grey", geom="point", shape=23, size=2, stroke=0.9)+
   geom_jitter(width=0.2, size=0.5)+     
   scale_fill_manual(values=col2)+
   ylab("Spearman's correlation")+ylim(-0.1, 0.5)+
   theme_bw()+
   theme(legend.position="none",
         axis.title.x=element_blank(),
         axis.title.y=element_text(size=9),
         axis.text.x=element_text(angle=45, hjust=1, size=9),
         axis.text.y=element_text(size=9))

comb_plots <- plot_grid(p1, p2, nrow=2, align="v", rel_heights=c(1, 0.9))
figfn <- paste(outdir, "Figure1.2_corr_gene_ATAC_box.png", sep="")
ggsave(figfn, comb_plots, device="png", width=320, height=560, units="px", dpi=120)


#####################################################################
### Scatter plots of SCC against #DEGs to compare which is better ###
#####################################################################

### setting colors
col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
  "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
   "6_Monocyte"="#984ea3", "7_dnT"="black")
col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
       "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3")


###
### spearman's correlation, option-2 
fn <- paste(outdir, "1.2_LFC_gene_ATAC_comb.rds", sep="")
df <- read_rds(fn)%>%drop_na(LFC_gene, LFC_ATAC) 
dfcorr <- df%>%group_by(MCls, contrast)%>%
    summarize(rr=cor(LFC_gene, LFC_ATAC, method="spearman"), .groups="drop")%>%as.data.frame()
dfcorr <- dfcorr%>%mutate(comb=paste(MCls, contrast, sep="_"))



###
### Number of DEGs
fn <- "./sc_multiome_data/2_Differential/1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn_sampleID+new_treat.rds"
resDEG <- read_rds(fn)%>%as.data.frame()%>%mutate(comb=paste(MCls, contrast, sep="_"))
resDEG2 <- resDEG%>%dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)

summ <-  resDEG2%>%group_by(MCls, contrast)%>%
    summarise(ngene=n(), .groups="drop")%>%as.data.frame()%>%mutate(comb=paste(MCls, contrast, sep="_"))
summ2 <- summ%>%dplyr::select(comb, ngene)




###
### scatter of plots of SCC of LFC of gene and peaks against #DEGs 

plotDF <- dfcorr%>%inner_join(summ2, by="comb")

corr <- cor.test(plotDF$rr, plotDF$ngene, method="spearman")
rr <- round(as.numeric(corr$estimate), digits=3)
eq <- deparse(bquote(italic(rho)==.(rr)~",***"))


p1 <- ggplot(plotDF, aes(x=ngene, y=rr, color=MCls))+
   geom_point(size=1.5)+
   annotate("text", label=eq, x=700, y=0.45, size=3, parse=T)+ 
   scale_color_manual(values=col1)+
   xlab("Number of DEGs")+xlim(0, 2000)+
   ylab("SCC of LFC on gene expression and peaks")+ylim(-0.1, 0.5)+ 
   theme_bw()+
   theme(legend.position="none",
         legend.title=element_blank(),
         legend.text=element_text(size=8),
         legend.key.size=unit(0.4, "cm"),
         axis.title=element_text(size=9),
         axis.text=element_text(size=8),
         plot.title=element_text(hjust=0.5, size=10))


###
### scatter plots against #DARs
 
### DARs
fn <- "./sc_multiome_data/2_Differential/1.2_DiffPeak.outs/2.0_DESeq.results_0.1_cn_sampleID+new_treat.rds"
resDP <- read_rds(fn)%>%as.data.frame()%>%mutate(comb=paste(MCls, contrast, sep="_"))
resDP2 <- resDP%>%dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)
summ2 <- resDP2%>%group_by(comb)%>%summarise(ngene=n(),.groups="drop")



plotDF2 <- dfcorr%>%inner_join(summ2, by="comb")
 
corr <- cor.test(plotDF2$rr, plotDF2$ngene, method="spearman")
rr <- round(as.numeric(corr$estimate), digits=3)
eq <- deparse(bquote(italic(rho)==.(rr)~",***"))

 
p2 <- ggplot(plotDF2, aes(x=ngene, y=rr, color=MCls))+
   geom_point(size=1.5)+
   annotate("text", label=eq, x=1500, y=0.45, size=3, parse=T)+ 
   scale_color_manual(values=col1)+
   xlab("Number of DARs")+xlim(0, 5000)+
   ylab("SCC of LFC on gene expression and peaks")+ylim(-0.1, 0.5)+ 
   theme_bw()+
   theme(##legend.position="none",
         legend.title=element_blank(),
         legend.text=element_text(size=8),
         legend.key.size=unit(0.4, "cm"),
         axis.title=element_text(size=9),
         axis.text=element_text(size=8),
         plot.title=element_text(hjust=0.5, size=10))

###
###
comb_plots <- plot_grid(p1, p2, ncol=2, align="h", rel_widths=c(1, 1.5))
figfn <- paste(outdir, "Figure1.3_scatter_SCCvsDEGs_DARs.png", sep="")
ggsave(figfn, comb_plots, device="png", width=800, height=380, units="px", dpi=120)


###




####################################################################
#### Enrichment analysis to examine if DEGs are enriched in DARs ###
####################################################################

rm(list=ls())

outdir <- "./3_compare_plots.outs/"


### DEG results
fn <- "./sc_multiome_data/2_Differential/1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn_sampleID+new_treat.rds"
resDG <- read_rds(fn)%>%as.data.frame()%>%mutate(comb=paste(MCls, contrast, sep="_"))


#### DP results
fn <- "./sc_multiome_data/2_Differential/1.2_DiffPeak.outs/2.0_DESeq.results_0.1_cn_sampleID+new_treat.rds"
resDP <- read_rds(fn)%>%as.data.frame()%>%mutate(comb=paste(MCls, contrast, sep="_"))
resDP <- resDP%>%dplyr::rename("peak"="gene")

###
fn <- "./sc_multiome_data/2_Differential/2.2_compareRNAandATAC.outs/1_annot.signac_macs2_0.1_cn.rds"
peakAnno <- read_rds(fn)
peakAnno <- peakAnno%>%dplyr::select(peak=query_region, gene_id=gene_name, dtss=distance)%>%
    mutate(peak=gsub(":", "-", peak), dtss=abs(dtss))



comb <- sort(unique(resDG$comb))
cvt <- str_split(comb, "_", simplify=T)
cvt <- data.frame(comb=comb, MCls=paste(cvt[,1], cvt[,2], sep="_"), contrast=cvt[,3])
DFcomb <- NULL
for (i in 1:nrow(cvt)){
   ### 
   ii <- cvt$comb[i]
   oneMCl <- cvt$MCls[i]
   contrast <- cvt$contrast[i]
    
   resDP2 <- resDP%>%
       dplyr::filter(comb==ii, p.adjusted<0.1, abs(estimate)>0.5)%>%
       left_join(peakAnno, by="peak")%>%dplyr::filter(dtss<1e+5)
   geneDP <- unique(resDP2$gene_id)

   cat(ii, "\n")
    
   if ( nrow(resDP2)>50){
       ###
       BG <- resDG%>%dplyr::filter(MCls==oneMCl)%>%pull(gene)%>%unique()
       DEG <- resDG%>%dplyr::filter(comb==ii, p.adjusted<0.1, abs(estimate)>0.5)%>%pull(gene)%>%unique()

       ## contingency table
       interest.in.DARs <- length(intersect(DEG, geneDP))
       interest.not.DARs <- length(DEG)-interest.in.DARs
       ###
       notDEG <- setdiff(BG, DEG)
       not.interest.in.DARs <- length(intersect(notDEG, geneDP))
       not.interest.not.DARs <- length(notDEG)-not.interest.in.DARs
       ###
       df2 <- data.frame("interest.in.DARs"=interest.in.DARs, "interest.not.DARs"=interest.not.DARs,
           "not.interest.in.DARs"=not.interest.in.DARs, "not.interest.not.DARs"=not.interest.not.DARs,
           comb=ii, MCls=oneMCl, "contrast"=contrast)
       dmat <- matrix(unlist(df2[1,1:4]), 2, 2)

       ### fisher exact test
       enrich <- fisher.test(dmat)
       ##enrich2 <- fisher.test(dmat, alternative="greater")       
       df2 <- df2%>%mutate(odds=enrich$estimate, pval=enrich$p.value,
           lower=enrich$conf.int[1], upper=enrich$conf.int[2])
       
       ###
       DFcomb <- rbind(DFcomb, df2)
    }### END IF       
}

DFcomb <- DFcomb%>%mutate(log_odds=log(odds), log_lower=log(lower), log_upper=log(upper))
DFcomb <- DFcomb%>%mutate(se=abs(log_odds-log_lower)/1.96, FDR=p.adjust(pval, method="BH"))

opfn <- paste(outdir, "2_enrich_Fisher.txt", sep="")
write.table(DFcomb, file=opfn, quote=F, row.names=F)

##
x2 <- DFcomb%>%filter(odds>1, FDR<0.05)%>%dplyr::select(comb, MCls, contrast, FDR, log_odds, se)


### next time skip the above enrich procedure and directly read data. 
fn <- "./3_compare_plots.outs/2_enrich_Fisher.txt"
DFcomb <- read.table(fn, header=T)
plotDF <- DFcomb
## plotDF <- cvt%>%dplyr::select(comb)%>%left_join(DFcomb, by="comb")


### setting color
col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
  "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
   "6_Monocyte"="#984ea3", "7_dnT"="black")
col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
       "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3")


### Fprestplots
plotDF <- plotDF%>%
    mutate(log_odds=ifelse(is.infinite(log_odds)|is.infinite(log_upper)|is.infinite(log_lower), NA, log_odds),
           log_lower=ifelse(is.na(log_odds), NA, log_lower),
           log_upper=ifelse(is.na(log_odds), NA, log_upper))

plotDF <- plotDF%>%drop_na(log_odds, log_lower, log_upper)

p <- ggplot(plotDF, aes(x=log_odds, y=comb))+
   geom_errorbarh(aes(xmax=log_upper, xmin=log_lower, colour=MCls),
       size=0.5, height=0.2)+
   geom_point(aes(colour=MCls), shape=19, size=1.5)+
   scale_colour_manual(values=col1)+
   geom_vline(aes(xintercept=0), size=0.25, linetype="dashed")+ 
   xlab("log odds ratio")+xlim(-3, 6)+
   theme_bw()+
   theme(##plot.title=element_text(hjust=0.5, size=14),
         axis.title.y=element_blank(),
         axis.title.x=element_text(size=10),
         axis.text.x=element_text(size=10),
         axis.text.y=element_text(size=8),
         legend.title=element_blank(),
         legend.text=element_text(size=8),
         legend.key.size=unit(0.4, "cm"))
         ##legend.position="none")
 
figfn <- paste(outdir, "Figure2_enrich_forest.png", sep="")
ggsave(figfn, p, device="png", width=620, height=550, units="px", dpi=120)



#######################################################################################
### A pair of heatmap shows the LFC on gene expression or peaks of each top 10 DEGs ###
#######################################################################################

### we don't use the script to generate any plots for paper. But the script is still useful for the similar plots 

## rm(list=ls())

## outdir <- "./3_compare_plots.outs/2_comb/"
## if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)


## ### DEG results
## fn <- "./sc_multiome_data/2_Differential/1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn_sampleID+new_treat.rds"
## resDG <- read_rds(fn)%>%as.data.frame()%>%
##     mutate(comb=paste(MCls, contrast, sep="_"),
##            is_sig=ifelse(p.adjusted<0.1&abs(estimate)>0.5, 1, 0))

## res2 <- resDG%>%dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5) 
## DEG <- unique(res2$gene)
## top10 <- res2%>%group_by(comb)%>%slice_max(order_by=abs(statistic), n=10)%>%ungroup()%>%pull(gene)%>%unique()
## ## length(intersect(top10, DEG))



## ### gene annotation of peaks 
## fn <- "./sc_multiome_data/2_Differential/2.2_compareRNAandATAC.outs/1_annot.signac_macs2_0.1_cn.rds"
## peakAnno <- read_rds(fn)
## peakAnno <- peakAnno%>%dplyr::select(peak=query_region, gene=gene_name, dtss=distance)%>%
##     mutate(peak=gsub(":", "-", peak), dtss=abs(dtss))

## peakAnno2 <- peakAnno%>%dplyr::filter(gene%in%top10, dtss<1e+5) 
## geneSel <- intersect(peakAnno2$gene, top10)

## ## x <- peakAnno2%>%group_by(gene)%>%summarize(ny=n(),.groups="drop")


## #### DP results
## fn <- "./sc_multiome_data/2_Differential/1.2_DiffPeak.outs/2.0_DESeq.results_0.1_cn_sampleID+new_treat.rds"
## resDP <- read_rds(fn)%>%as.data.frame()%>%
##     mutate(comb=paste(MCls, contrast, sep="_"), is_sig=ifelse(p.adjusted<0.1&abs(estimate)>0.5, 1, 0))
## resDP <- resDP%>%dplyr::rename("peak"="gene")
## peakSel <- resDP%>%filter(is_sig==1)%>%pull(peak)%>%unique()


######################
### heatmap of DEG ###
######################
 
## df2 <- resDG%>%dplyr::filter(gene%in%geneSel)
## dfmat <- df2%>%pivot_wider(id_cols=gene, names_from=comb, values_from=statistic)
## mat <- dfmat%>%column_to_rownames(var="gene")%>%as.matrix()

## ## mat2 <- mat

## rnz <- rowSums(is.na(mat))
## mat2 <- mat[rnz==0,]


## ###
## ### get colnames and re-order by treats
## rn <- colnames(mat2)
## x <- str_split(rn, "_", simplify=T)
## cvt <- data.frame(comb=rn, MCls=paste(x[,1], x[,2], sep="_"), contrast=x[,3])
## cvt <- cvt%>%arrange(contrast)


## ################################
## ### select the top 100 genes ###
## ################################

## x <- mat2
## vx <- diag(var(t(x)))
## vx <- sort(vx, T)

## geneSel2 <- names(vx[1:80])
 
## ## xx <- x%*%t(x)
## ## xx.eigen <- eigen(xx)
## ## dd <- xx.eigen$values
## ## U <- xx.eigen$vectors

## ### reorder column and select some variable gene
## mat2 <- mat2[geneSel2, cvt$comb]



## ### color for heatmap value
## y <- as.vector(mat2)
## ##y0 <- y[abs(y)<2] ### quantile(abs(y), probs=0.99)
## mybreak <- c(min(y,na.rm=T), seq(-10, 10, length.out=98), max(y,na.rm=T))

## quantile(abs(y), probs=c(0.9,0.95,0.99))
 
## mycol <- colorRamp2(mybreak, colorRampPalette(rev(brewer.pal(n=7,name="RdBu")))(100))

 
## ###
## ### annotation columns
## col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
##   "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
##    "6_Monocyte"="#984ea3", "7_dnT"="black")
## col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
##        "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3")


## x <- str_split(colnames(mat2), "_", simplify=T)
## df_col <- data.frame(celltype=paste(x[,1], x[,2], sep="_"), contrast=x[,3])
## col_ha <- HeatmapAnnotation(df=df_col, col=list(celltype=col1, contrast=col2),
##     annotation_name_gp=gpar(fontsize=10),
##     annotation_legend_param=list(
##   celltype=list(title_gp=gpar(fontsize=9), labels_gp=gpar(fontsize=9), title="celltype",
##                 grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm")),
##   contrast=list(title_gp=gpar(fontsize=9), labels_gp=gpar(fontsize=9), title="contrast",
##                 grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm"))))
##   ##show_legend=c(F,F),
##   ##simple_anno_size=unit(0.3, "cm"))


## ###
## ### main heatmap

## fn <- "./3_compare_plots.outs/3.0_gene_reorder.xlsx"
## cl_sort <- read.xlsx(fn)
## rn <- cl_sort$gene

## mat3 <- mat2[rn,]
## p1 <- Heatmap(mat3, col=mycol,
##    cluster_rows=F, cluster_columns=F,
##    row_split=c(rep("1", 6), rep("2", 35), rep("3", 21), rep("4", 18)), 
##    show_row_names=T, row_names_gp=gpar(fontsize=6),
##    show_column_names=T, column_names_gp=gpar(fontsize=6),
##    show_row_dend=F, show_column_dend=F,
##    ##column_names_rot=-45,   
##    top_annotation=col_ha,
##    heatmap_legend_param=list(title=bquote(~italic(Z)~"-score"),
##       title_gp=gpar(fontsize=9),
##       at=seq(-12, 12, by=4), 
##       labels_gp=gpar(fontsize=9),
##       grid_width=grid::unit(0.45, "cm"),
##       legend_height=grid::unit(7.8, "cm")))
 
## ###
## figfn <- paste(outdir, "Figure3.1_zscore_gene.heatmap.png", sep="")
## ## ggsave(figfn, plot=p2, device="png", width=680, height=700, units="px", dpi=300) ## ggsave 
## png(figfn, width=850, height=950,res=120)
## set.seed(0)
## p1 <- draw(p1)
## dev.off()


## ##############
## ### is_sig ###
## ##############
 
## rn <- rownames(mat3)
## comb <- colnames(mat3)
## df2 <- resDG%>%dplyr::filter(gene%in%rn)
## is_sig <- df2%>%pivot_wider(id_cols=gene, names_from=comb, values_from=is_sig)%>%
##     column_to_rownames(var="gene")%>%as.matrix()
## ###
## is_sig <- is_sig[rn,comb]
## is_sig[is.na(is_sig)] <- 0
## mat3_sig <- mat3*is_sig

## p1_b <- Heatmap(mat3_sig, col=mycol,
##    cluster_rows=F, cluster_columns=F,
##    row_split=c(rep("1", 6), rep("2", 35), rep("3", 21), rep("4", 18)), 
##    show_row_names=T, row_names_gp=gpar(fontsize=6),
##    show_column_names=T, column_names_gp=gpar(fontsize=6),
##    show_row_dend=F, show_column_dend=F,
##    ##column_names_rot=-45,   
##    top_annotation=col_ha,
##    heatmap_legend_param=list(title=bquote(~italic(Z)~"-score"),
##       title_gp=gpar(fontsize=9),
##       at=seq(-12, 12, by=4), 
##       labels_gp=gpar(fontsize=9),
##       grid_width=grid::unit(0.45, "cm"),
##       legend_height=grid::unit(7.8, "cm")))
 
## ###
## figfn <- paste(outdir, "Figure3.2_zscore_gene_isDEG.heatmap.png", sep="")
## ## ggsave(figfn, plot=p2, device="png", width=680, height=700, units="px", dpi=300) ## ggsave 
## png(figfn, width=850, height=950,res=120)
## set.seed(0)
## p1_b <- draw(p1_b)
## dev.off()




## ### cluster gene and re-order gene
## hmap <- Heatmap(mat2, cluster_rows=T, row_km=4, cluster_columns=F)
## set.seed(0)
## hmap <- draw(hmap)
## cl <- row_order(hmap)


## geneSel <- rownames(mat2)
## DF_cl <- NULL
## for (i in names(cl)){
##    cl_tmp <- data.frame(cluster=i, gene=geneSel[cl[[i]]])
##    DF_cl <- rbind(DF_cl, cl_tmp)
## }
## ##
## newCL <- c("4"="1", "3"="2", "2"="3", "1"="4")
## DF_cl <- DF_cl%>%mutate(cluster2=newCL[cluster])%>%arrange(cluster2)

## opfn <- paste(outdir, "3.0_gene_reorder.xlsx", sep="")
## write.xlsx(DF_cl, file=opfn, overwrite=T)
 


## ########################
## #### Heatmap for DPs ###
## ########################

## fn <- "3_compare_plots.outs/3.0_gene_reorder.xlsx"
## DF_cl <- read.xlsx(fn)
## geneSel <- DF_cl$gene

## resDP2 <- resDP%>%inner_join(peakAnno2, by="peak")
## resDP2 <- resDP2%>%group_by(comb, gene)%>%
##     slice_max(order_by=abs(statistic), n=1)%>%ungroup()%>%as.data.frame()
 
## ##df2 <- resDF%>%dplyr::filter(gene%in%top10)
## mat <- resDP2%>%
##     pivot_wider(id_cols=gene, names_from=comb, values_from=statistic)%>%
##     column_to_rownames(var="gene")%>%as.matrix()
 
## mat2 <- mat[geneSel,]

## ## rnz <- rowSums(is.na(mat))
## ## mat2 <- as.matrix(mat[rnz==0,])
## ## rownames(mat2) <- dfmat$gene[rnz==0]


## ###
## ### re-order by treats
## rn <- colnames(mat2)
## x <- str_split(rn, "_", simplify=T)
## cvt <- data.frame(comb=rn, MCls=paste(x[,1], x[,2], sep="_"), contrast=x[,3])
## cvt <- cvt%>%arrange(contrast)

## mat2 <- mat2[, cvt$comb]

## rnz <- rowSums(is.na(mat2))
## mat2 <- mat2[rnz==0,]


## ### color for heatmap value
## y <- as.vector(mat2)
## ##y0 <- y[abs(y)<2] ### quantile(abs(y), probs=0.99)
## mybreak <- c(min(y), seq(-10, 10, length.out=98), max(y)) 
## mycol <- colorRamp2(mybreak, colorRampPalette(rev(brewer.pal(n=7,name="RdBu")))(100))

## quantile(abs(y), probs=c(0.9, 0.95, 0.99))
 
## ###
## ### annotation columns
## col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
##   "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
##    "6_Monocyte"="#984ea3", "7_dnT"="black")
## col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
##        "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3")


## x <- str_split(colnames(mat2), "_", simplify=T)
## df_col <- data.frame(celltype=paste(x[,1], x[,2], sep="_"), contrast=x[,3])
## col_ha <- HeatmapAnnotation(df=df_col, col=list(celltype=col1, contrast=col2),
##     annotation_name_gp=gpar(fontsize=10),
##     annotation_legend_param=list(
##   celltype=list(title_gp=gpar(fontsize=9), labels_gp=gpar(fontsize=9), title="celltype",
##                 grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm")),
##   contrast=list(title_gp=gpar(fontsize=9), labels_gp=gpar(fontsize=9), title="contrast",
##                 grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm"))),
##   show_legend=c(F,F))

##   ## simple_anno_size=unit(0.3, "cm"))


## ###
## ### main heatmap  
## p2 <- Heatmap(mat2, col=mycol,
##    cluster_rows=F, cluster_columns=F,
##    row_split=c(rep("1", 6), rep("2", 35), rep("3", 21), rep("4", 18)), 
##    show_row_names=T, row_names_gp=gpar(fontsize=6),
##    show_column_names=T, column_names_gp=gpar(fontsize=6),
##    show_row_dend=F,
##    ##column_names_rot=-45,   
##    top_annotation=col_ha,
##    heatmap_legend_param=list(title=bquote(~italic(Z)~"-score"),
##       title_gp=gpar(fontsize=9),
##       at=seq(-12, 12, by=4), 
##       labels_gp=gpar(fontsize=9),
##       grid_width=grid::unit(0.45, "cm"),
##       legend_height=grid::unit(7.8, "cm")))


   
## ###
## figfn <- paste(outdir, "Figure4.1_zscore_peak.heatmap.png", sep="")
## ## ggsave(figfn, plot=p2, device="png", width=680, height=700, units="px", dpi=300) ## ggsave 
## png(figfn, width=700, height=950,res=120)
## p2 <- draw(p2)
## dev.off()


## #######################
## ### if significance ###
## #######################
## resDP2 <- resDP2%>%mutate(is_sig=ifelse(p.adjusted<0.1&abs(estimate)>0.5, 1, 0))

## rn <- rownames(mat2)
## comb <- colnames(mat2)
## df2 <- resDP2%>%dplyr::filter(gene%in%rn)
## is_sig <- df2%>%pivot_wider(id_cols=gene, names_from=comb, values_from=is_sig)%>%
##     column_to_rownames(var="gene")%>%as.matrix()
## ###
## is_sig <- is_sig[rn,comb]
## is_sig[is.na(is_sig)] <- 0
## mat2_sig <- mat2*is_sig


## p2_b <- Heatmap(mat2_sig, col=mycol,
##    cluster_rows=F, cluster_columns=F,
##    row_split=c(rep("1", 6), rep("2", 35), rep("3", 21), rep("4", 18)), 
##    show_row_names=T, row_names_gp=gpar(fontsize=6),
##    show_column_names=T, column_names_gp=gpar(fontsize=6),
##    show_row_dend=F,
##    ##column_names_rot=-45,   
##    top_annotation=col_ha,
##    heatmap_legend_param=list(title=bquote(~italic(Z)~"-score"),
##       title_gp=gpar(fontsize=9),
##       at=seq(-12, 12, by=4), 
##       labels_gp=gpar(fontsize=9),
##       grid_width=grid::unit(0.45, "cm"),
##       legend_height=grid::unit(7.8, "cm")))


## ###
## figfn <- paste(outdir, "Figure4.2_zscore_peak_isDAR.heatmap.png", sep="")
## ## ggsave(figfn, plot=p2, device="png", width=680, height=700, units="px", dpi=300) ## ggsave 
## png(figfn, width=700, height=950,res=120)
## p2_b <- draw(p2_b)
## dev.off()





###############################################
### combine LFC on gene expression and ATAC ### 
###############################################

### Data from the above part

rm(list=ls())

outdir2 <- "./3_compare_plots.outs/2_comb/"
if ( !file.exists(outdir2) ) dir.create(outdir2, showWarnings=F, recursive=T)


### DEG results
fn <- "./sc_multiome_data/2_Differential/1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn_sampleID+new_treat.rds"
resDG <- read_rds(fn)%>%as.data.frame()%>%
    mutate(comb=paste(MCls, contrast, sep="_"),
           is_sig=ifelse(p.adjusted<0.1&abs(estimate)>0.5, 1, 0))

DEG <- resDG%>%dplyr::filter(is_sig==1)%>%pull(gene)%>%unique() 
## geneSel <- DEG


top10 <- resDG%>%filter(is_sig==1)%>%
    group_by(comb)%>%slice_max(order_by=abs(statistic), n=10)%>%ungroup()%>%pull(gene)%>%unique()
geneSel <- intersect(top10, DEG)




### gene annotation of peaks 
fn <- "./sc_multiome_data/2_Differential/2.2_compareRNAandATAC.outs/1_annot.signac_macs2_0.1_cn.rds"
peakAnno <- read_rds(fn)
peakAnno <- peakAnno%>%dplyr::select(peak=query_region, gene=gene_name, dtss=distance)%>%
    mutate(peak=gsub(":", "-", peak), dtss=abs(dtss))

## peakAnno2 <- peakAnno%>%dplyr::filter(gene%in%DEG, dtss<1e+5) 
## geneSel <- intersect(peakAnno2$gene, top10)

## x <- peakAnno2%>%group_by(gene)%>%summarize(ny=n(),.groups="drop")


#### DP results
fn <- "./sc_multiome_data/2_Differential/1.2_DiffPeak.outs/2.0_DESeq.results_0.1_cn_sampleID+new_treat.rds"
resDP <- read_rds(fn)%>%as.data.frame()%>%
    mutate(comb=paste(MCls, contrast, sep="_"), is_sig=ifelse(p.adjusted<0.1&abs(estimate)>0.5, 1, 0))
resDP <- resDP%>%dplyr::rename("peak"="gene")
peakSel <- resDP%>%filter(is_sig==1)%>%pull(peak)%>%unique()


### LFC on gene 
df2 <- resDG%>%dplyr::filter(gene%in%geneSel)
dfmat <- df2%>%pivot_wider(id_cols=gene, names_from=comb, values_from=statistic)
mat <- dfmat%>%column_to_rownames(var="gene")%>%as.matrix()
comb <- paste0("1_", colnames(mat))
colnames(mat) <- comb


### LFC on gene chromatin 
resDP2 <- resDP%>%inner_join(peakAnno, by="peak")%>%
     filter(abs(dtss)<1e+05, gene%in%geneSel)%>%mutate(zscore=ifelse(peak%in%peakSel, statistic, NA))    

### median value
resDP3 <- resDP2%>%group_by(comb, gene)%>%
    summarise(zscore=median(zscore, na.rm=T),.groups="drop")%>%ungroup()%>%as.data.frame()

### maximum value
## resDP3 <- resDP2%>%group_by(comb, gene)%>%
##     slice_max(order_by=abs(statistic), n=1)%>%ungroup()%>%as.data.frame()%>%
##     dplyr::rename(zscore=statistic)
 

mat2 <- resDP3%>%
    pivot_wider(id_cols=gene, names_from=comb, values_from=zscore)%>%
    column_to_rownames(var="gene")%>%as.matrix()

rn <- intersect(rownames(mat), rownames(mat2))
mat2 <- mat2[rn,]
mat <- mat[rn,]

comb <- paste0("2_", colnames(mat2))
colnames(mat2) <- comb

identical(rownames(mat), rownames(mat2))


mat3 <- cbind(mat, mat2)


###
### data frame of colnames 
x <- str_split(colnames(mat3), "_", simplify=T)

x2 <- data.frame(comb=colnames(mat3), Feature=as.numeric(x[,1]),
   MCls=paste(x[,2], x[,3], sep="_"), contrast=x[,4])
x2 <- x2%>%arrange(Feature, contrast)

mat3 <- mat3[, x2$comb]

rnz <- rowSums(is.na(mat3))
mat3 <- mat3[rnz==0,]

opfn <- paste(outdir, "comb2_zscore_RNA_ATAC.rds", sep="")
write_rds(mat3, file=opfn)





#############
### plots ###
#############

rm(list=ls())

outdir2 <- "./3_compare_plots.outs/2_comb/"
if ( !file.exists(outdir2) ) dir.create(outdir2, showWarnings=F, recursive=T)

fn <- paste(outdir2, "comb2_zscore_RNA_ATAC_149genes.rds", sep="")
mat3 <- read_rds(fn)
 

### color for heatmap value
y <- as.vector(mat3)
##y0 <- y[abs(y)<2] ### quantile(abs(y), probs=0.99)
mybreak <- c(min(y), seq(-8, 8, length.out=98), max(y)) 
mycol <- colorRamp2(mybreak, colorRampPalette(rev(brewer.pal(n=7,name="RdBu")))(100))

quantile(abs(y), probs=c(0.9, 0.95, 0.99))


###
### annotation columns
col0 <- c("1"="#8c510a", "2"="#d8b365")
col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
  "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
   "6_Monocyte"="#984ea3", "7_dnT"="black")
col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
       "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3")


x <- str_split(colnames(mat3), "_", simplify=T)
x2 <- data.frame(comb=colnames(mat3), Feature=as.numeric(x[,1]),
   MCls=paste(x[,2], x[,3], sep="_"), contrast=x[,4])
x2 <- x2%>%arrange(Feature, contrast)


df_col <- x2%>%dplyr::select(MCls, contrast) ##%>%mutate(Feature=as.character(Feature))

col_ha <- HeatmapAnnotation(df=df_col, col=list(MCls=col1, contrast=col2),
    annotation_name_gp=gpar(fontsize=10),
    annotation_legend_param=list(       
  MCls=list(title_gp=gpar(fontsize=9), labels_gp=gpar(fontsize=9), title="celltype",
                grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm")),
  contrast=list(title_gp=gpar(fontsize=9), labels_gp=gpar(fontsize=9), title="contrast",
                grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm"))))

comb <- gsub("^[12]_", "", colnames(mat3))
colnames(mat3) <- comb


## show_legend=c(F,F))
### main plots
p2 <- Heatmap(mat3, col=mycol,
   cluster_rows=T, cluster_columns=F, row_km=7, column_split=rep(c("1 RNA", "2 ATAC"), each=48),
   show_row_names=T, row_names_gp=gpar(fontsize=4),
   show_column_names=T, column_names_gp=gpar(fontsize=6),
   show_row_dend=T, row_dend_width=grid::unit(2, "cm"),
   show_column_dend=F,
   column_names_rot=-45,   
   top_annotation=col_ha,
   heatmap_legend_param=list(title=bquote(~italic(Z)~"-score"),
      title_gp=gpar(fontsize=9),
      at=seq(-8, 8, by=4), 
      labels_gp=gpar(fontsize=9),
      grid_width=grid::unit(0.45, "cm"),
      legend_height=grid::unit(7.8, "cm")))
 
###
figfn <- paste(outdir2, "Figure1.2_top10_cl7_zscore_comb2.heatmap.png", sep="")
## ggsave(figfn, plot=p2, device="png", width=680, height=700, units="px", dpi=300) ## ggsave 
png(figfn, width=1600, height=1200,res=120)
set.seed(0)
p2 <- draw(p2)
dev.off()


###
### cluster id

hmap <- Heatmap(mat3, cluster_rows=T, row_km=7, cluster_columns=F)
set.seed(0)
hmap <- draw(hmap)
cl <- row_order(hmap)
lengths(cl)

geneSel <- rownames(mat3)
DF_cl <- NULL
for (i in names(cl)){
   cl_tmp <- data.frame(cluster=i, gene=geneSel[cl[[i]]])
   DF_cl <- rbind(DF_cl, cl_tmp)
}

opfn <- paste(outdir2, "1.2_cl7_gene_reorder.xlsx", sep="")
write.xlsx(DF_cl, file=opfn, overwrite=T)
 







###########################################################################
### Scatter plots of LFC on gene expression and chromatin accessibility ###
###########################################################################


###
### plot data

fn <- paste(outdir, "1.2_LFC_gene_ATAC_comb.rds", sep="")
df2 <- read_rds(fn)%>%drop_na(LFC_gene, LFC_ATAC) 


###
### First example
ii <- "6_Monocyte_vitD"
df3 <- df2%>%dplyr::filter(comb==ii)


###
### if it is DEG or not 
fn <- "./sc_multiome_data/2_Differential/1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn_sampleID+new_treat.rds"
resDG <- read_rds(fn)%>%as.data.frame()%>%
    mutate(comb=paste(MCls, contrast, sep="_"),
           is_sig=ifelse(p.adjusted<0.1&abs(estimate)>0.5, 1, 0),
           is_sig=ifelse(is.na(is_sig), 0, is_sig))

resDG2 <- resDG%>%filter(comb==ii)%>%dplyr::select(gene, is_sig)

### combine plot data with 
df3 <- df3%>%left_join(resDG2, by="gene")


###
### add annotation text
corr <- cor.test(df3$LFC_ATAC, df3$LFC_gene, method="spearman")
rr <- round(as.numeric(corr$estimate), digits=3)
pval <- corr$p.value

symb <- case_when(pval<=1e-03~"***",
          pval>1e-03&pval<=0.01~"**",
          pval>1e-01&pval<=0.05~"*",
          TRUE~"NS")
eq <- deparse(bquote(italic(rho)==.(rr)~.(symb)))


### modification 
df3 <- df3%>%mutate(LFC_ATAC=ifelse(LFC_ATAC<(-5), -5, LFC_ATAC),
                    LFC_ATAC=ifelse(LFC_ATAC>(5), 5, LFC_ATAC),
                    LFC_gene=ifelse(LFC_gene<(-5), 5, LFC_gene),
                    LFC_gene=ifelse(LFC_gene>5, 5, LFC_gene))
df3 <- as.data.frame(df3)

y0 <- 5
x0 <- -3.7
  
p0 <- ggplot(df3, aes(x=LFC_ATAC, y=LFC_gene))+
   geom_point(aes(color=factor(is_sig)), size=0.3)+
   scale_color_manual(values=c("grey", "red"), labels=c("Not DEG", "DEG"),
      guide=guide_legend(override.aes=list(size=1.5)))+                
   annotate("text", label=eq, x=x0, y=y0, size=3.5, parse=T)+ 
   scale_x_continuous("LFC on gene chromatin accessibility", limits=c(-6,6))+
   scale_y_continuous("LFC on gene expression", limits=c(-6,6))+
   geom_vline(xintercept=0, linetype="dashed", size=0.4)+
   geom_hline(yintercept=0, linetype="dashed", size=0.4)+ 
   ## geom_smooth(method="lm", formula=y~x, size=0.5, se=F)+ 
   ggtitle("vitamin D in Monocyte")+
   theme_bw()+ 
   theme(legend.position=c(0.2, 0.8),
         legend.title=element_blank(),
         legend.text=element_text(size=8),
         legend.key.size=unit(0.3, "cm"),
         legend.background=element_blank(),
         legend.key=element_blank(),
         legend.box.background=element_blank(),
         axis.title=element_text(size=10),
         axis.text=element_text(size=10),
         plot.title=element_text(hjust=0.5, size=12))        

figfn <- paste(outdir, "Figure1.4_", ii, "_scatter_LFC.png", sep="")
ggsave(figfn, p0, device="png", width=380, height=380, units="px", dpi=120)


###
### another example

###
### plot data, combine the information about if it is DEG or not

ii <- "6_Monocyte_zinc"
df3 <- df2%>%dplyr::filter(comb==ii)
resDG2 <- resDG%>%filter(comb==ii)%>%dplyr::select(gene, is_sig)

df3 <- df3%>%left_join(resDG2, by="gene")


## df3 <- df3%>%mutate(gr2=ifelse(grepl("^MT", gene), 1, 0))

###
### add annotation text
corr <- cor.test(df3$LFC_ATAC, df3$LFC_gene, method="spearman")
rr <- round(as.numeric(corr$estimate), digits=3)
pval <- corr$p.value

symb <- case_when(pval<=1e-03~"***",
          pval>1e-03&pval<=0.01~"**",
          pval>1e-01&pval<=0.05~"*",
          TRUE~"NS")
eq <- deparse(bquote(italic(rho)==.(rr)~.(symb)))

y0 <- 4.2
x0 <- -3.8

 
### modification plot data
df3 <- df3%>%mutate(LFC_ATAC=ifelse(LFC_ATAC<(-5), -5, LFC_ATAC),
                    LFC_ATAC=ifelse(LFC_ATAC>(5), 5, LFC_ATAC),
                    LFC_gene=ifelse(LFC_gene<(-5), 5, LFC_gene),
                    LFC_gene=ifelse(LFC_gene>5, 5, LFC_gene))
df3 <- as.data.frame(df3)


###
### annotate genes 
anno2 <- df3%>%filter(LFC_gene==5, is_sig>0, grepl("^MT", gene))

###
### main figures
p0 <- ggplot(df3, aes(x=LFC_ATAC, y=LFC_gene))+
   geom_point(aes(color=factor(is_sig)), size=0.3)+
   scale_color_manual(values=c("grey", "red"), labels=c("Not DEG", "DEG"),
      guide=guide_legend(override.aes=list(size=1.5)))+                
   annotate("text", label=eq, x=x0, y=y0, size=3.5, parse=T)+
   geom_text_repel(data=anno2,
       aes(x=LFC_ATAC, y=LFC_gene, label=gene),
       max.overlaps=30, box.padding=0.3, color="maroon3", fontface="italic", size=3)+ 
   scale_x_continuous("LFC on gene chromatin accessibility", limits=c(-6, 6))+
   scale_y_continuous("LFC on gene expression", limits=c(-6,6))+
   geom_vline(xintercept=0, linetype="dashed", size=0.4)+
   geom_hline(yintercept=0, linetype="dashed", size=0.4)+ 
   ## geom_smooth(method="lm", formula=y~x, size=0.5, se=F)+ 
   ggtitle("zinc in Monocyte")+
   theme_bw()+ 
   theme(legend.position="none",
         legend.title=element_blank(),
         legend.text=element_text(size=8),
         legend.key.size=unit(0.4, "cm"),
         legend.background=element_blank(),
         legend.key=element_blank(),
         legend.box.background=element_blank(),
         axis.title=element_text(size=10),
         axis.text=element_text(size=10),
         plot.title=element_text(hjust=0.5, size=12))        

figfn <- paste(outdir, "Figure1.4_", ii, "_scatter_LFC.png", sep="")
ggsave(figfn, p0, device="png", width=380, height=380, units="px", dpi=120)




################################
#### Heatmap of correlation ####
################################


rm(list=ls())

outdir <- "./3_compare_plots.outs/"
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)



fn <- paste(outdir, "1.2_LFC_gene_ATAC_comb.rds", sep="")
df <- read_rds(fn)%>%drop_na(LFC_gene, LFC_ATAC) 

dfcorr <- df%>%group_by(MCls, contrast)%>%
    summarize(rr=cor(LFC_gene, LFC_ATAC, method="spearman"), .groups="drop")%>%as.data.frame()

dfmat <- dfcorr%>%pivot_wider(id_cols=MCls, names_from=contrast, values_from=rr)%>%as.data.frame()

mat <- dfmat%>%column_to_rownames(var="MCls")
mat <- as.matrix(mat)


###
### significance
sig <- df%>%group_by(MCls, contrast)%>%
    summarize(pval=cor.test(LFC_gene, LFC_ATAC,  method="spearman")$p.value, .groups="drop")%>%as.data.frame()
sig <- sig%>%mutate(FDR=p.adjust(pval, method="BH"), is_sig=ifelse(FDR<0.05, 1, 0))

sig_mat <- sig%>%pivot_wider(id_cols=MCls, names_from=contrast, values_from=is_sig)%>%as.data.frame()
 
sig_mat <- sig_mat%>%column_to_rownames(var="MCls")
sig_mat <- as.matrix(sig_mat)

mat2 <- mat*sig_mat

mat2[mat2==0] <- NA


### setting color
## mybreak <- seq(0, 1, length.out=20)
## col0 <- brewer.pal(n=9,name="Reds")
## mycol <- colorRamp2(mybreak, colorRampPalette(col0)(20))

mycol <- colorRamp2(seq(-1, 1, length.out=50), colorRampPalette(rev(brewer.pal(n=7,name="RdBu")))(50))



## mycol <- colorRamp2(seq(0, 1, length.out=12), colorRampPalette(c("white", "red"))(12))

p <- Heatmap(mat2, name="SCC", na_col="grey90", 
    col=mycol, cluster_rows=F, cluster_columns=F,
    row_names_gp=gpar(fontsize=10), column_names_gp=gpar(fontsize=10), ##column_names_rot=-45,
    heatmap_legend_param=list(at=seq(-1, 1, by=0.5),
        grid_width=grid::unit(0.38, "cm"), legend_height=grid::unit(4, "cm"),
        title_gp=gpar(fontsize=8), labels_gp=gpar(fontsize=8)),
    cell_fun=function(j, i, x, y, width, height, fill){
       grid.text(round(mat[i,j],digits=3), x, y, gp=gpar(fontsize=10))
    })

figfn <- paste(outdir, "Figure1.5_corr_heatmap.png", sep="")
###ggsave(figfn, p, width=520, height=520, units="px",dpi=120)
png(figfn, height=520, width=520, res=120)
print(p)
dev.off()    
 


#############################
### correlation for genes ###
#############################

rm(list=ls())

outdir2 <- "./3_compare_plots.outs/2_comb/"
if ( !file.exists(outdir2) ) dir.create(outdir2, showWarnings=F, recursive=T)

### gene cluster
fn <- "./3_compare_plots.outs/2_comb/1.2_cl7_gene_reorder.xlsx"
df_cl <- read.xlsx(fn)
geneSel <- unique(df_cl$gene)

###
fn <- "./3_compare_plots.outs/1.2_LFC_gene_ATAC_comb.rds"
df2 <- read_rds(fn)

###
### LFC of gene
x1 <- df2%>%pivot_wider(id_cols=gene, names_from=comb, values_from=LFC_gene, values_fill=NA)%>%
    column_to_rownames(var="gene")%>%as.matrix()
x1 <- x1[geneSel,]


###
### LFC of gene chromatin
x2 <- df2%>%pivot_wider(id_cols=gene, names_from=comb, values_from=LFC_ATAC, values_fill=NA)%>%
    column_to_rownames(var="gene")%>%as.matrix()
x2 <- x2[geneSel,]

sum(rowSums(is.na(x2))>0)


identical(rownames(x1), rownames(x2))
identical(colnames(x1), colnames(x2))


###
### correlation, SCC

rr <- diag(cor(t(x1), t(x2), method="spearman"))
df_rr <- data.frame(gene=names(rr), rr=as.numeric(rr))
df_rr <- df_rr%>%inner_join(df_cl, by="gene")

opfn <- paste(outdir2, "1.2_corr.xlsx", sep="")
write.xlsx(df_rr, file=opfn)




###
### box plot show correlation

fn <- "3_compare_plots.outs/2_comb/1.2_corr.xlsx"
df_rr <- read.xlsx(fn)%>%mutate(cluster2=paste("CL", cluster, sep=""))


p <- ggplot(df_rr, aes(x=cluster2, y=rr, color=cluster2))+
   geom_boxplot(outlier.shape=NA, lwd=0.25)+
   geom_jitter(width=0.2, size=0.5)+ 
   stat_summary(fun=median, geom="point", shape=23, size=1.8, stroke=0.9)+ylim(-1,1)+
   ggtitle("SCC of LFC between RNA and ATAC")+  
   theme_bw()+
   theme(legend.position="none",
         axis.title=element_blank(),
         axis.text=element_text(size=10),
         plot.title=element_text(size=10,hjust=0.5))

figfn <- paste(outdir2, "Figure2.0_corr.box.png", sep="")
ggsave(figfn, p, width=450, height=350, units="px", dpi=120)




####################
### CCA analysis ###
####################


rm(list=ls())

outdir2 <- "./3_compare_plots.outs/2_comb/"
if ( !file.exists(outdir2) ) dir.create(outdir2, showWarnings=F, recursive=T)


### DEG results
fn <- "./sc_multiome_data/2_Differential/1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn_sampleID+new_treat.rds"
resDG <- read_rds(fn)%>%as.data.frame()%>%
    mutate(comb=paste(MCls, contrast, sep="_"),
           is_sig=ifelse(p.adjusted<0.1&abs(estimate)>0.5, 1, 0))


top10 <- resDG%>%filter(is_sig==1)%>%
    group_by(comb)%>%slice_max(order_by=abs(statistic), n=10)%>%ungroup()%>%pull(gene)%>%unique()


fn <- "./3_compare_plots.outs/2_comb/comb2_zscore_RNA_ATAC.rds"
mat <- read_rds(fn)
rnz <- rowSums(is.na(mat))
mat <- mat[rnz==0,]   ## 2975 genes without any missing values

geneSel <- intersect(rownames(mat), top10)

###
mat2 <- mat[geneSel,]
###
opfn <- paste(outdir2, "comb2_zscore_RNA_ATAC_149genes.rds", sep="")
write_rds(mat2, file=opfn)



####
### top 10 genes

###
### input data
fn <- paste(outdir2, "comb2_zscore_RNA_ATAC_149genes.rds", sep="")
x <- read_rds(fn)
comb2 <- gsub("^[12]_", "", colnames(x))
colnames(x) <- comb2


###
x1 <- x[,1:48]
x1 <- scale(x1)
###
x2 <- x[,49:96]
x2 <- scale(x2)


###
### cca analysis
cc0 <- cc(x1, x2)
cc2 <- comput(x1, x2, cc0)



####
####

###
rna <- as.data.frame(cc2$xscores)
names(rna) <- paste("s", 1:48, sep="")
###
peak <- as.data.frame(cc2$yscores)
names(peak) <- paste("s", 1:48, sep="")


identical(rownames(rna), rownames(peak))



###
### custom function 
my_box <- function(data, mapping, shape,...){
   ##
   p <- ggplot(data=data, mapping=mapping)+
       geom_boxplot(outlier.shape=NA, lwd=0.25)+coord_flip()
   ## geom_jitter(data=data, mapping=mapping, width=0.2, size=0.5, shape=23)+coord_flip()## + 
   ## stat_summary(fun=median, geom="point", shape=23, size=1.8, stroke=0.9) ###ylim(-1,1)+
   p
}    


fn <- "./3_compare_plots.outs/2_comb/1.2_corr.xlsx"
df_cl <- read.xlsx(fn)%>%mutate(cluster2=paste("CL", cluster, sep=""))


###
### The correlation of score from motif and TF-genes

df2 <- data.frame(x=1:48, rr=cc0$cor)
###
p0 <- ggplot(df2, aes(x=x, y=rr))+
   geom_point(size=1.5, shape=1, color="#dd1c77")+
   guides(color=guide_legend(override.aes=list(size=2)))+ 
   xlab("Canonical variable")+ylab("Correlation (RNA-ATAC)")+
   theme_bw()+
   theme(axis.title=element_text(size=9),
         axis.text=element_text(size=9))

figfn <- paste(outdir2, "Figure2.1_corr_score.scatter.png", sep="")
ggsave(figfn, p0, width=420, height=420, units="px", dpi=120)





###
### plot data for RNA

plotDF <- rna[,1:5]%>%rownames_to_column(var="gene")%>%
    inner_join(df_cl, by="gene")%>%
    mutate(cluster2=paste("CL", cluster, sep=""))


###
###
p0 <- ggpairs(plotDF, columns=2:6, aes(color=cluster2),
              upper=list(continuous=wrap("cor", size=2)),
              diag=list(continuous=my_box),
              lower=list(continuous=wrap("points", size=0.8)))+
    theme_bw()+
    theme(strip.text=element_text(size=12))
    

figfn <- paste(outdir2, "Figure2.2_RNA_score.pair.png", sep="")
ggsave(figfn, plot=p0, width=800, height=720, units="px", dpi=120)



###
### plot data for peak

plotDF2 <- peak[,1:5]%>%rownames_to_column(var="gene")%>%
    inner_join(df_cl, by="gene")%>%
    mutate(cluster2=paste("CL", cluster, sep=""))


###
###
p2 <- ggpairs(plotDF2, columns=2:6, aes(color=cluster2),
              upper=list(continuous=wrap("cor", size=2)),
              diag=list(continuous=my_box),
              lower=list(continuous=wrap("points", size=0.8)))+
    theme_bw()+
    theme(strip.text=element_text(size=12))
    
figfn <- paste(outdir2, "Figure2.2_peak_score.pair.png", sep="")
ggsave(figfn, plot=p2, width=780, height=720, units="px", dpi=120)




##############################
### Heatmap score from cca ###
##############################

m1 <- cc2$xscores
m2 <- cc2$yscores

colnames(m1) <- paste("S", 1:48, sep="")
colnames(m2) <- paste("S", 1:48, sep="")

mat <- cbind(m1, m2)



#### setting colors
y <- as.vector(mat)
##y0 <- y[abs(y)<2] ### quantile(abs(y), probs=0.99)
mybreak <- c(min(y), seq(-3, 3, length.out=98), max(y)) 
mycol <- colorRamp2(mybreak, colorRampPalette(rev(brewer.pal(n=7,name="RdBu")))(100))

quantile(abs(y), probs=c(0.9, 0.95, 0.99))

### main plots
p2 <- Heatmap(mat, col=mycol,
   cluster_rows=T, row_km=6, cluster_columns=F,
   column_split=rep(c("1 RNA", "2 ATAC"), each=48),
   show_row_dend=T, row_dend_width=grid::unit(2, "cm"),   
   show_row_names=T, row_names_gp=gpar(fontsize=6),
   show_column_names=T, column_names_gp=gpar(fontsize=6),
   ## show_column_names=F, column_names_gp=gpar(fontsize=6),
   ###top_annotation=col_ha,
   heatmap_legend_param=list(title="score",
      title_gp=gpar(fontsize=9),
      at=seq(-4.5, 4.5, by=1.5), 
      labels_gp=gpar(fontsize=9),
      grid_width=grid::unit(0.45, "cm"),
      legend_height=grid::unit(7.8, "cm")))
 
###
figfn <- paste(outdir2, "Figure2.3_score_comb2.heatmap.png", sep="")
## ggsave(figfn, plot=p2, device="png", width=680, height=700, units="px", dpi=300) ## ggsave 
png(figfn, width=1600, height=1200,res=120)
set.seed(0)
p2 <- draw(p2)
dev.off()



########################
### loading from cca ###
########################

m1 <- t(cc2$corr.X.xscores)
m2 <- t(cc2$corr.Y.yscores)

rownames(m1) <- paste("S", 1:48, sep="")
rownames(m2) <- paste("S", 1:48, sep="")


###
### loading of RNA

#### setting colors
y <- as.vector(m1)
##y0 <- y[abs(y)<2] ### quantile(abs(y), probs=0.99)
## mybreak <- c(min(y), seq(-1, 1, length.out=98), max(y))
mybreak <- seq(-1, 1, length.out=100)
mycol <- colorRamp2(mybreak, colorRampPalette(rev(brewer.pal(n=7,name="RdBu")))(100))

## quantile(abs(y), probs=c(0.9, 0.95, 0.99))


###
### annotation columns
col0 <- c("1"="#8c510a", "2"="#d8b365")
col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
  "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
   "6_Monocyte"="#984ea3", "7_dnT"="black")
col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
       "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3")


cvt <- str_split(colnames(m1), "_", simplify=T)
cvt2 <- data.frame(comb=colnames(m1),
   MCls=paste(cvt[,1], cvt[,2], sep="_"), contrast=cvt[,3])
cvt2 <- cvt2%>%arrange(contrast)


df_col <- cvt2%>%dplyr::select(MCls, contrast) ##%>%mutate(Feature=as.character(Feature))

col_ha <- HeatmapAnnotation(df=df_col, col=list(MCls=col1, contrast=col2),
    annotation_name_gp=gpar(fontsize=10),
    annotation_legend_param=list(       
  MCls=list(title_gp=gpar(fontsize=9), labels_gp=gpar(fontsize=9), title="celltype",
                grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm")),
  contrast=list(title_gp=gpar(fontsize=9), labels_gp=gpar(fontsize=9), title="contrast",
                grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm"))))

### plot data
m1a <- m1[,cvt2$comb]


### main plots
p1 <- Heatmap(m1a, col=mycol,
   cluster_rows=F, cluster_columns=F,
   show_column_dend=F, column_dend_height=grid::unit(2, "cm"),   
   show_row_names=T, row_names_gp=gpar(fontsize=7),
   show_column_names=T, column_names_gp=gpar(fontsize=6),
   top_annotation=col_ha,
   heatmap_legend_param=list(title="loading",
      title_gp=gpar(fontsize=9),
      at=seq(-1, 1, by=0.5), 
      labels_gp=gpar(fontsize=9),
      grid_width=grid::unit(0.45, "cm"),
      legend_height=grid::unit(7.8, "cm")))
 
###
figfn <- paste(outdir2, "Figure2.4_loading_RNA.heatmap.png", sep="")
png(figfn, width=720, height=720,res=120)
set.seed(0)
p1 <- draw(p1)
dev.off()



###
### loading of peak

#### setting colors
y <- as.vector(m2)
##y0 <- y[abs(y)<2] ### quantile(abs(y), probs=0.99)
## mybreak <- c(min(y), seq(-1, 1, length.out=98), max(y))
mybreak <- seq(-1, 1, length.out=100)
mycol <- colorRamp2(mybreak, colorRampPalette(rev(brewer.pal(n=7,name="RdBu")))(100))

## quantile(abs(y), probs=c(0.9, 0.95, 0.99))


###
### annotation columns
col0 <- c("1"="#8c510a", "2"="#d8b365")
col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
  "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
   "6_Monocyte"="#984ea3", "7_dnT"="black")
col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
       "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3")


cvt <- str_split(colnames(m2), "_", simplify=T)
cvt2 <- data.frame(comb=colnames(m2),
   MCls=paste(cvt[,1], cvt[,2], sep="_"), contrast=cvt[,3])
cvt2 <- cvt2%>%arrange(contrast)


df_col <- cvt2%>%dplyr::select(MCls, contrast) ##%>%mutate(Feature=as.character(Feature))

col_ha <- HeatmapAnnotation(df=df_col, col=list(MCls=col1, contrast=col2),
    annotation_name_gp=gpar(fontsize=10),
    annotation_legend_param=list(       
  MCls=list(title_gp=gpar(fontsize=9), labels_gp=gpar(fontsize=9), title="celltype",
                grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm")),
  contrast=list(title_gp=gpar(fontsize=9), labels_gp=gpar(fontsize=9), title="contrast",
                grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm"))))



### plot data
m2a <- m2[,cvt2$comb]

### main plots
p2 <- Heatmap(m2a, col=mycol,
   cluster_rows=F, cluster_columns=F,
   show_column_dend=F, column_dend_height=grid::unit(2, "cm"),   
   show_row_names=T, row_names_gp=gpar(fontsize=7),
   show_column_names=T, column_names_gp=gpar(fontsize=6),
   top_annotation=col_ha,
   heatmap_legend_param=list(title="loading",
      title_gp=gpar(fontsize=9),
      at=seq(-1, 1, by=0.5), 
      labels_gp=gpar(fontsize=9),
      grid_width=grid::unit(0.45, "cm"),
      legend_height=grid::unit(7.8, "cm")))
 
###
figfn <- paste(outdir2, "Figure2.4_loading_peak.heatmap.png", sep="")
## ggsave(figfn, plot=p2, device="png", width=680, height=700, units="px", dpi=300) ## ggsave 
png(figfn, width=720, height=720,res=120)
set.seed(0)
p2 <- draw(p2)
dev.off()





################################
### scatter plots of loading ###
################################

### custom function 
my_box <- function(data, mapping,...){
   ##
   p <- ggplot(data=data, mapping=mapping)+
       geom_boxplot(outlier.shape=NA, lwd=0.25)+coord_flip()
   ## geom_jitter(data=data, mapping=mapping, width=0.2, size=0.5, shape=23)+coord_flip()## + 
   ## stat_summary(fun=median, geom="point", shape=23, size=1.8, stroke=0.9) ###ylim(-1,1)+
   p
}


m1 <- as.data.frame(cc2$corr.X.xscores)
m2 <- as.data.frame(cc2$corr.Y.yscores)

colnames(m1) <- paste("S", 1:48, sep="")
colnames(m2) <- paste("S", 1:48, sep="")


col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
  "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
   "6_Monocyte"="#984ea3", "7_dnT"="black")
col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
       "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3")


###
### plot data for motif

plotDF <- m1[,1:5]%>%as.data.frame()%>%rownames_to_column(var="comb")
x <- str_split(plotDF$comb, "_", simplify=T)
plotDF <- plotDF%>%mutate(MCls=paste(x[,1], x[,2], sep="_"), contrast=x[,3])


###
###
p0 <- ggpairs(plotDF, columns=2:6, aes(color=contrast),
              upper=list(continuous=wrap("cor", size=2)),
              diag=list(continuous=my_box),
              lower=list(continuous=wrap("points", size=0.8)))+
    scale_color_manual(values=col2)+
    theme_bw()+
    theme(strip.text=element_text(size=12))
    

figfn <- paste(outdir2, "Figure2.5_RNA_loading.pair.png", sep="")
ggsave(figfn, plot=p0, width=800, height=720, units="px", dpi=120)



###
### plot data for TF-regulated genes

plotDF2 <- m2[,1:5]%>%as.data.frame()%>%rownames_to_column(var="comb")
x <- str_split(plotDF2$comb, "_", simplify=T)
plotDF2 <- plotDF2%>%mutate(MCls=paste(x[,1], x[,2], sep="_"), contrast=x[,3])


###
###
p2 <- ggpairs(plotDF2, columns=2:6, aes(color=contrast),
              upper=list(continuous=wrap("cor", size=2)),
              diag=list(continuous=my_box),
              lower=list(continuous=wrap("points", size=0.8)))+
    scale_color_manual(values=col2)+
    theme_bw()+
    theme(strip.text=element_text(size=12))
    

figfn <- paste(outdir2, "Figure2.5_peak_loading.pair.png", sep="")
ggsave(figfn, plot=p2, width=780, height=720, units="px", dpi=120)







    







## xs2 <- cc0$scores$xscores
## ys2 <- cc0$scores$yscores

###

