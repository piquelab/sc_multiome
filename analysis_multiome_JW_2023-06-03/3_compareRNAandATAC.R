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

rm(list=ls())

outdir <- "./3_compare_plots.outs/"
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)




#########################################################################
### Box plots of correlation between LFC on gene expression and peaks ###
#########################################################################

### 1.0_*, median all the peaks
### 1.2_*, median union DARs ## used for final analysis
### 1.3_*, median union DARs in the cell type


###############################################
### (1). correlation of LFC on RNA and ATAC ###
###############################################

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

 
comb <- sort(unique(resDP$comb))
dfcomb <- map_dfr(comb, function(ii){
   ###
   cat(ii, "\n") 
   x <- resDP%>%filter(comb==ii)%>%dplyr::select(estimate, peak)
   x2 <- x%>%left_join(peakAnno, by="peak")
   x2 <- x2 %>%filter(dtss<1e+05) ## 100 kb
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


## 1.0, LFC of gene chromatin accessbility calculated by median of all the peaks
opfn <- paste(outdir, "1.0_LFC_gene_ATAC_comb.rds", sep="")
write_rds(dfcomb, file=opfn)



####
#### box plot show correlation


fn <- paste(outdir, "1.0_LFC_gene_ATAC_comb.rds", sep="")
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


## dfcorr <- dfcorr%>%mutate(rr=round(rr,3))
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

## figfn <- paste(outdir, "Figure1.0_corr_gene_ATAC_box.png", sep="")
## ggsave(figfn, p, device="png", width=350, height=320, units="px", dpi=120)

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
figfn <- paste(outdir, "Figure1.0_corr_gene_ATAC_box.png", sep="")
ggsave(figfn, comb_plots, device="png", width=320, height=560, units="px", dpi=120)



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

## figfn <- paste(outdir, "Figure1.0_corr_gene_ATAC_box.png", sep="")
## ggsave(figfn, p, device="png", width=350, height=320, units="px", dpi=120)

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



##############################################
#### (3) cell type specific DEGs and DARs ####
##############################################

comb <- sort(unique(resDP$comb))
cvt <- str_split(comb, "_", simplify=T)
cvt <- data.frame(comb=comb, MCls=paste(cvt[,1], cvt[,2], sep="_"), contrast=cvt[,3])

###
dfcomb <- map_dfr(1:nrow(cvt), function(i){
   ###
   ii <- cvt$comb[i]
   oneMCl <- cvt$MCls[i]
   DEG <- resDG%>%dplyr::filter(MCls==oneMCl, p.adjusted<0.1, abs(estimate)>0.5)%>%pull(gene)%>%unique()
   DP <- resDP%>%dplyr::filter(MCls==oneMCl, p.adjusted<0.1, abs(estimate)>0.5)%>%pull(peak)%>%unique()
   cat(ii, length(DEG), length(DP), "\n")
    
   ### 
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

opfn <- paste(outdir, "1.3_LFC_gene_ATAC_comb.rds", sep="")
write_rds(dfcomb, file=opfn)



####
#### box plot show correlation

fn <- paste(outdir, "1.3_LFC_gene_ATAC_comb.rds", sep="")
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
   ylab("Spearman's correlation")+ylim(-0.25, 0.7)+
   theme_bw()+
   theme(legend.position="none",
         axis.title.x=element_blank(),
         axis.title.y=element_text(size=9),
         axis.text.x=element_text(angle=45, hjust=1, size=9),
         axis.text.y=element_text(size=9))

## figfn <- paste(outdir, "Figure1.0_corr_gene_ATAC_box.png", sep="")
## ggsave(figfn, p, device="png", width=350, height=320, units="px", dpi=120)

###
### group by treat
p2 <- ggplot(dfcorr, aes(x=contrast, y=rr, fill=contrast))+
   ##geom_violin()+ 
   geom_boxplot(outlier.shape=NA, lwd=0.25)+
   stat_summary(fun=median, color="grey", geom="point", shape=23, size=2, stroke=0.9)+
   geom_jitter(width=0.2, size=0.5)+     
   scale_fill_manual(values=col2)+
   ylab("Spearman's correlation")+ylim(-0.25, 0.7)+
   theme_bw()+
   theme(legend.position="none",
         axis.title.x=element_blank(),
         axis.title.y=element_text(size=9),
         axis.text.x=element_text(angle=45, hjust=1, size=9),
         axis.text.y=element_text(size=9))

comb_plots <- plot_grid(p1, p2, nrow=2, align="v", rel_heights=c(1, 0.9))
figfn <- paste(outdir, "Figure1.3_corr_gene_ATAC_box.png", sep="")
ggsave(figfn, comb_plots, device="png", width=320, height=560, units="px", dpi=120)



###
### Scatter plots of SCC against #DEGs to compare which is better

### setting colors
col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
  "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
   "6_Monocyte"="#984ea3", "7_dnT"="black")
col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
       "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3")


fn <- "./sc_multiome_data/2_Differential/1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn_sampleID+new_treat.rds"
res <- read_rds(fn)%>%as.data.frame()%>%mutate(comb=paste(MCls, contrast, sep="_"))
res2 <- res%>%dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)

summ <-  res2%>%group_by(MCls, contrast)%>%
    summarise(ngene=n(), .groups="drop")%>%as.data.frame()%>%mutate(comb=paste(MCls, contrast, sep="_"))
summ2 <- summ%>%dplyr::select(comb, ngene)

 
###
### spearman's correlation, option-2 
fn <- paste(outdir, "1.2_LFC_gene_ATAC_comb.rds", sep="")
df <- read_rds(fn)%>%drop_na(LFC_gene, LFC_ATAC) 
dfcorr <- df%>%group_by(MCls, contrast)%>%
    summarize(rr=cor(LFC_gene, LFC_ATAC, method="spearman"), .groups="drop")%>%as.data.frame()
dfcorr <- dfcorr%>%mutate(comb=paste(MCls, contrast, sep="_"))


plotDF <- dfcorr%>%inner_join(summ2, by="comb")

corr <- cor.test(plotDF$rr, plotDF$ngene, method="spearman")
rr <- round(as.numeric(corr$estimate), digits=3)
eq <- deparse(bquote(italic(rho)==.(rr)~",***"))


p1 <- ggplot(plotDF, aes(x=ngene, y=rr, color=MCls))+
   geom_point(size=1.5)+
   annotate("text", label=eq, x=350, y=0.45, parse=T)+ 
   scale_color_manual(values=col1)+
   ggtitle("Option 1")+ 
   xlab("Number of DEGs")+xlim(0, 2000)+
   ylab("SCC of LFC on gene expression and peaks")+ylim(-0.1, 0.5)+ 
   theme_bw()+
   theme(##legend.position="none",
         legend.title=element_blank(),
         legend.text=element_text(size=8),
         legend.key.size=unit(0.4, "cm"),
         axis.title=element_text(size=9),
         axis.text=element_text(size=9),
         plot.title=element_text(hjust=0.5, size=10))


## p1 <- p1+theme(plot.title=element_blank())
## figfn <- paste(outdir, "Figure1.4.1_compare_scatter.png", sep="")
## ggsave(figfn, p1, device="png", width=600, height=380, units="px", dpi=120)

 

###
### option-3
fn <- paste(outdir, "1.3_LFC_gene_ATAC_comb.rds", sep="")
df <- read_rds(fn)%>%drop_na(LFC_gene, LFC_ATAC) 
dfcorr <- df%>%group_by(MCls, contrast)%>%
    summarize(rr=cor(LFC_gene, LFC_ATAC, method="spearman"), .groups="drop")%>%as.data.frame()
dfcorr <- dfcorr%>%mutate(comb=paste(MCls, contrast, sep="_"))


plotDF2 <- dfcorr%>%inner_join(summ2, by="comb")

corr <- cor.test(plotDF2$rr, plotDF2$ngene, method="spearman")
rr <- round(as.numeric(corr$estimate), digits=3)
eq <- deparse(bquote(italic(rho)==.(rr)~",***"))

p2 <- ggplot(plotDF2, aes(x=ngene, y=rr, color=MCls))+
   geom_point(size=1.5)+
   annotate("text", label=eq, x=300, y=0.6, parse=T)+ 
   scale_color_manual(values=col1)+
   xlab("Number of DEGs")+xlim(0, 2000)+
   ylab("SCC of LFC on gene expression and ATAC")+ylim(-0.25, 0.7)+
   ggtitle("Option 2")+ 
   theme_bw()+
   theme(legend.title=element_blank(),
         legend.text=element_text(size=8),
         legend.key.size=unit(0.4, "cm"),
         axis.title=element_text(size=10),
         axis.text=element_text(size=9),
         plot.title=element_text(hjust=0.5, size=10))

###
###
comb_plots <- plot_grid(p1, p2, ncol=2, align="h", rel_widths=c(1, 1.5))
figfn <- paste(outdir, "Figure1.4_compare_scatter.png", sep="")
ggsave(figfn, comb_plots, device="png", width=720, height=380, units="px", dpi=120)



###
### scatter plots against #DARs

### DARs
fn <- "./sc_multiome_data/2_Differential/1.2_DiffPeak.outs/2.0_DESeq.results_0.1_cn_sampleID+new_treat.rds"
resDP <- read_rds(fn)%>%as.data.frame()%>%mutate(comb=paste(MCls, contrast, sep="_"))
res2 <- resDP%>%dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)
summ2 <- res2%>%group_by(comb)%>%summarise(ngene=n(),.groups="drop")


###
### spearman's correlation, option-2 
fn <- paste(outdir, "1.2_LFC_gene_ATAC_comb.rds", sep="")
df <- read_rds(fn)%>%drop_na(LFC_gene, LFC_ATAC) 
dfcorr <- df%>%group_by(MCls, contrast)%>%
    summarize(rr=cor(LFC_gene, LFC_ATAC, method="spearman"), .groups="drop")%>%as.data.frame()
dfcorr <- dfcorr%>%mutate(comb=paste(MCls, contrast, sep="_"))


plotDF <- dfcorr%>%inner_join(summ2, by="comb")

corr <- cor.test(plotDF$rr, plotDF$ngene, method="spearman")
rr <- round(as.numeric(corr$estimate), digits=3)
eq <- deparse(bquote(italic(rho)==.(rr)~",***"))

 
p1 <- ggplot(plotDF, aes(x=ngene, y=rr, color=MCls))+
   geom_point(size=1.5)+
   annotate("text", label=eq, x=1500, y=0.45, size=3, parse=T)+ 
   scale_color_manual(values=col1)+
   ggtitle("Option 1")+ 
   xlab("Number of DARs")+xlim(0, 5000)+
   ylab("SCC of LFC on gene expression and peaks")+ylim(-0.1, 0.5)+ 
   theme_bw()+
   theme(legend.position="none",
         legend.title=element_blank(),
         legend.text=element_text(size=8),
         legend.key.size=unit(0.4, "cm"),
         axis.title=element_text(size=9),
         axis.text=element_text(size=8),
         plot.title=element_text(hjust=0.5, size=10))


## p1 <- p1+theme(plot.title=element_blank())
## figfn <- paste(outdir, "Figure1.4.1_compare_scatter.png", sep="")
## ggsave(figfn, p1, device="png", width=600, height=380, units="px", dpi=120)

 

###
### option-3
fn <- paste(outdir, "1.3_LFC_gene_ATAC_comb.rds", sep="")
df <- read_rds(fn)%>%drop_na(LFC_gene, LFC_ATAC) 
dfcorr <- df%>%group_by(MCls, contrast)%>%
    summarize(rr=cor(LFC_gene, LFC_ATAC, method="spearman"), .groups="drop")%>%as.data.frame()
dfcorr <- dfcorr%>%mutate(comb=paste(MCls, contrast, sep="_"))


plotDF2 <- dfcorr%>%inner_join(summ2, by="comb")

corr <- cor.test(plotDF2$rr, plotDF2$ngene, method="spearman")
rr <- round(as.numeric(corr$estimate), digits=3)
eq <- deparse(bquote(italic(rho)==.(rr)~",***"))

p2 <- ggplot(plotDF2, aes(x=ngene, y=rr, color=MCls))+
   geom_point(size=1.5)+
   annotate("text", label=eq, x=1500, y=0.75, size=3, parse=T)+ 
   scale_color_manual(values=col1)+
   xlab("Number of DEGs")+xlim(0, 5000)+
   ylab("SCC of LFC on gene expression and ATAC")+ylim(-0.25, 0.8)+
   ggtitle("Option 2")+ 
   theme_bw()+
   theme(legend.title=element_blank(),
         legend.text=element_text(size=8),
         legend.key.size=unit(0.4, "cm"),
         axis.title=element_text(size=10),
         axis.text=element_text(size=8),
         plot.title=element_text(hjust=0.5, size=10))

###
###
comb_plots <- plot_grid(p1, p2, ncol=2, align="h", rel_widths=c(1, 1.5))
figfn <- paste(outdir, "Figure1.5_compare_scatter_DARs.png", sep="")
ggsave(figfn, comb_plots, device="png", width=780, height=380, units="px", dpi=120)










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
       enrich2 <- fisher.test(dmat, alternative="greater")       
       df2 <- df2%>%mutate(odds=enrich$estimate, pval=enrich2$p.value,
           lower=enrich$conf.int[1], upper=enrich$conf.int[2])
       
       ###
       DFcomb <- rbind(DFcomb, df2)
    }### END IF       
}

DFcomb <- DFcomb%>%mutate(log_odds=log(odds), log_lower=log(lower), log_upper=log(upper))

opfn <- paste(outdir, "2_enrich_Fisher.txt", sep="")
write.table(DFcomb, file=opfn, quote=F, row.names=F)


### next time skip the above enrich procedure and directly read data. 
fn <- "./3_compare_plots.outs/2_enrich_Fisher.txt"
DFcomb <- read.table(fn, header=T)
plotDF <- cvt%>%dplyr::select(comb)%>%left_join(DFcomb, by="comb")


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


rm(list=ls())

outdir <- "./3_compare_plots.outs/"


### DEG results
fn <- "./sc_multiome_data/2_Differential/1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn_sampleID+new_treat.rds"
resDG <- read_rds(fn)%>%as.data.frame()%>%mutate(comb=paste(MCls, contrast, sep="_"))
res2 <- resDG%>%dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5) 
DEG <- unique(res2$gene)
top10 <- res2%>%group_by(comb)%>%slice_max(order_by=abs(statistic), n=10)%>%ungroup()%>%pull(gene)%>%unique()
##top10 <- intersect(top10, DEG) ### 64 genes
###length(intersect(top10, DEG))



### gene annotation of peaks 
fn <- "./sc_multiome_data/2_Differential/2.2_compareRNAandATAC.outs/1_annot.signac_macs2_0.1_cn.rds"
peakAnno <- read_rds(fn)
peakAnno <- peakAnno%>%dplyr::select(peak=query_region, gene=gene_name, dtss=distance)%>%
    mutate(peak=gsub(":", "-", peak), dtss=abs(dtss))

peakAnno2 <- peakAnno%>%dplyr::filter(gene%in%top10, dtss<1e+5) 
geneSel <- intersect(peakAnno2$gene, top10)

## x <- peakAnno2%>%group_by(gene)%>%summarize(ny=n(),.groups="drop")


#### DP results
fn <- "./sc_multiome_data/2_Differential/1.2_DiffPeak.outs/2.0_DESeq.results_0.1_cn_sampleID+new_treat.rds"
resDP <- read_rds(fn)%>%as.data.frame()%>%mutate(comb=paste(MCls, contrast, sep="_"))
resDP <- resDP%>%dplyr::rename("peak"="gene")



######################
### heatmap of DEG ###
######################
 
df2 <- resDG%>%dplyr::filter(gene%in%geneSel)
dfmat <- df2%>%pivot_wider(id_cols=gene, names_from=comb, values_from=statistic) 
mat <- as.matrix(dfmat[,-1])
rownames(mat) <- dfmat$gene

## mat2 <- mat

rnz <- rowSums(is.na(mat))
mat2 <- as.matrix(mat[rnz==0,])
rownames(mat2) <- dfmat$gene[rnz==0]

###
### get colnames and re-order by treats
rn <- colnames(mat2)
x <- str_split(rn, "_", simplify=T)
cvt <- data.frame(comb=rn, MCls=paste(x[,1], x[,2], sep="_"), contrast=x[,3])
cvt <- cvt%>%arrange(contrast)


################################
### select the top 100 genes ###
################################

x <- mat2
vx <- diag(var(t(x)))
vx <- sort(vx, T)

geneSel2 <- names(vx[1:100])

## xx <- x%*%t(x)
## xx.eigen <- eigen(xx)
## dd <- xx.eigen$values
## U <- xx.eigen$vectors

### reorder column and select some variable gene
mat2 <- mat2[geneSel2, cvt$comb]



### color for heatmap value
y <- as.vector(mat2)
##y0 <- y[abs(y)<2] ### quantile(abs(y), probs=0.99)
mybreak <- c(min(y,na.rm=T), seq(-10, 10, length.out=98), max(y,na.rm=T))

quantile(abs(y), probs=c(0.9,0.95,0.99))
 
mycol <- colorRamp2(mybreak, colorRampPalette(rev(brewer.pal(n=7,name="RdBu")))(100))

 
###
### annotation columns
col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
  "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
   "6_Monocyte"="#984ea3", "7_dnT"="black")
col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
       "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3")


x <- str_split(colnames(mat2), "_", simplify=T)
df_col <- data.frame(celltype=paste(x[,1], x[,2], sep="_"), contrast=x[,3])
col_ha <- HeatmapAnnotation(df=df_col, col=list(celltype=col1, contrast=col2),
    annotation_name_gp=gpar(fontsize=10),
    annotation_legend_param=list(
  celltype=list(title_gp=gpar(fontsize=9), labels_gp=gpar(fontsize=9), title="celltype",
                grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm")),
  contrast=list(title_gp=gpar(fontsize=9), labels_gp=gpar(fontsize=9), title="contrast",
                grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm"))),
  show_legend=c(F,F),
  simple_anno_size=unit(0.3, "cm"))


###
### main heatmap

p1 <- Heatmap(mat2, col=mycol,
   cluster_rows=T, cluster_columns=F,
   show_row_names=T, row_names_gp=gpar(fontsize=4.8),
   show_column_names=T, column_names_gp=gpar(fontsize=6),
   show_row_dend=F, show_column_dend=F,
   ##column_names_rot=-45,   
   top_annotation=col_ha,
   heatmap_legend_param=list(title=bquote(~italic(Z)~"-score"),
      title_gp=gpar(fontsize=9),
      at=seq(-12, 12, by=4), 
      labels_gp=gpar(fontsize=9),
      grid_width=grid::unit(0.45, "cm"),
      legend_height=grid::unit(7.8, "cm")))
 
###
figfn <- paste(outdir, "Figure3.1_LFC_gene.heatmap.png", sep="")
## ggsave(figfn, plot=p2, device="png", width=680, height=700, units="px", dpi=300) ## ggsave 
png(figfn, width=700, height=1000,res=120)
set.seed(0)
p1 <- draw(p1)
dev.off()





### cluster gene and re-order gene
hmap <- Heatmap(mat2, cluster_rows=T, cluster_columns=F)
set.seed(0)
hmap <- draw(hmap)
cl <- row_order(hmap)

geneSel <- rownames(mat2)
DF_cl <- NULL
for (i in 1:length(cl)){
   cl_tmp <- data.frame(cluster=i, gene=geneSel[cl[[i]]])
   DF_cl <- rbind(DF_cl, cl_tmp)
}
opfn <- paste(outdir, "3.0_gene_reorder.xlsx", sep="")
write.xlsx(DF_cl, file=opfn, overwrite=T)



########################
#### Heatmap for DPs ###
########################

fn <- "3_compare_plots.outs/3.0_gene_reorder.xlsx"
geneSel <- read.xlsx(fn)$gene

resDP2 <- resDP%>%inner_join(peakAnno2, by="peak")
resDP2 <- resDP2%>%group_by(comb, gene)%>%
    slice_max(order_by=abs(statistic), n=1)%>%ungroup()%>%as.data.frame()
 
##df2 <- resDF%>%dplyr::filter(gene%in%top10)
dfmat <- resDP2%>%pivot_wider(id_cols=gene, names_from=comb, values_from=statistic) 
mat <- as.matrix(dfmat[,-1])
rownames(mat) <- dfmat$gene
 
mat2 <- mat[geneSel,]

## rnz <- rowSums(is.na(mat))
## mat2 <- as.matrix(mat[rnz==0,])
## rownames(mat2) <- dfmat$gene[rnz==0]


###
### re-order by treats
rn <- colnames(mat2)
x <- str_split(rn, "_", simplify=T)
cvt <- data.frame(comb=rn, MCls=paste(x[,1], x[,2], sep="_"), contrast=x[,3])
cvt <- cvt%>%arrange(contrast)

mat2 <- mat2[, cvt$comb]

rnz <- rowSums(is.na(mat2))
mat2 <- mat2[rnz==0,]

### color for heatmap value
y <- as.vector(mat2)
##y0 <- y[abs(y)<2] ### quantile(abs(y), probs=0.99)
mybreak <- c(min(y), seq(-10, 10, length.out=98), max(y)) 
mycol <- colorRamp2(mybreak, colorRampPalette(rev(brewer.pal(n=7,name="RdBu")))(100))

quantile(abs(y), probs=c(0.9, 0.95, 0.99))
 
###
### annotation columns
col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
  "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
   "6_Monocyte"="#984ea3", "7_dnT"="black")
col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
       "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3")


x <- str_split(colnames(mat2), "_", simplify=T)
df_col <- data.frame(celltype=paste(x[,1], x[,2], sep="_"), contrast=x[,3])
col_ha <- HeatmapAnnotation(df=df_col, col=list(celltype=col1, contrast=col2),
    annotation_name_gp=gpar(fontsize=10),
    annotation_legend_param=list(
  celltype=list(title_gp=gpar(fontsize=9), labels_gp=gpar(fontsize=9), title="celltype",
                grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm")),
  contrast=list(title_gp=gpar(fontsize=9), labels_gp=gpar(fontsize=9), title="contrast",
                grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm"))),
  simple_anno_size=unit(0.3, "cm"))


###
### main heatmap  
p2 <- Heatmap(mat2, col=mycol,
   cluster_rows=F, cluster_columns=F,
   show_row_names=T, row_names_gp=gpar(fontsize=4.8),
   show_column_names=T, column_names_gp=gpar(fontsize=6),
   show_row_dend=F,
   ##column_names_rot=-45,   
   top_annotation=col_ha,
   heatmap_legend_param=list(title=bquote(~italic(Z)~"-score"),
      title_gp=gpar(fontsize=9),
      at=seq(-12, 12, by=4), 
      labels_gp=gpar(fontsize=9),
      grid_width=grid::unit(0.45, "cm"),
      legend_height=grid::unit(7.8, "cm")))


   
###
figfn <- paste(outdir, "Figure3.2_LFC_peak.heatmap.png", sep="")
## ggsave(figfn, plot=p2, device="png", width=680, height=700, units="px", dpi=300) ## ggsave 
png(figfn, width=850, height=1000,res=120)
p2 <- draw(p2)
dev.off()
