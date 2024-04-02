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
library(annotables)

## library(MASS)
## library("rainbow", lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages")
## library("fds", lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages")
## library(fda, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages")
## library(CCA, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages")

## library("GGally")


###
### clean script, in the final we used union of DARs to calculate LFC of gene chromatin accessibility 


###################################################
### generate LFC on genes and LFC on peak 
###################################################


rm(list=ls())
outdir <- "./Plots_pub/3_compare_plots.outs/"
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)


###
### Read data

#### DEGs
fn <- "./sc_multiome_data/2_Differential/1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn_sampleID+new_treat.rds"
resDG <- read_rds(fn)%>%as.data.frame()%>%mutate(comb=paste(MCls, contrast, sep="_"))

gene0 <- unique(resDG$gene)
gene1 <- grch38%>%dplyr::filter(chr%in%as.character(1:22))%>%pull(symbol)%>%unique()
gene2 <- intersect(gene0, gene1)


DEG <- resDG%>%dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5, gene%in%gene2)%>%pull(gene)%>%unique()
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

 

###
### calcualte median value of LFC of DP ###
comb <- sort(unique(resDP$comb))
dfcomb <- map_dfr(comb, function(ii){
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
   cat(ii, nrow(df2), "\n") 
   df2
})


opfn <- paste(outdir, "1.2_LFC_gene_ATAC.comb.rds", sep="")
write_rds(dfcomb, file=opfn)


###
###

## fn <- "./3_compare_plots.outs/1.2_LFC_gene_ATAC_comb.rds"
## dfcomb <- read_rds(fn)
## dfcomb2 <- dfcomb%>%filter(gene%in%gene2)




#########################################################################
### Box plots of correlation between LFC on gene expression and peaks ###
#########################################################################

rm(list=ls())
outdir <- "./Plots_pub/3_compare_plots.outs/"
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)



####
#### box plot show correlation

fn <- paste(outdir, "1.2_LFC_gene_ATAC.comb.rds", sep="")
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
   stat_summary(fun=median, color="grey", geom="point", shape=23, size=2.5, stroke=0.9)+
   geom_jitter(width=0.2, size=0.5)+     
   scale_fill_manual(values=col1)+
   ## ggtitle("Between gene expression and chromatin")+ 
   ylab("SCC")+ylim(-0.1, 0.4)+
   theme_bw()+
   theme(legend.position="none",
         axis.title.x=element_blank(),
         axis.title.y=element_text(size=10),
         axis.text.x=element_blank(),
         axis.text.y=element_text(size=10),
         axis.ticks.x=element_blank(),
         plot.title=element_text(size=10))
         ## axis.title.y=element_text(size=9),
         ## axis.text.x=element_text(angle=30, hjust=1, size=8),
         ## axis.text.y=element_text(size=9))

###
### group by treat
p2 <- ggplot(dfcorr, aes(x=contrast, y=rr, fill=contrast))+
   ##geom_violin()+ 
   geom_boxplot(outlier.shape=NA, lwd=0.25)+
   stat_summary(fun=median, color="grey", geom="point", shape=23, size=2.5, stroke=0.9)+
   geom_jitter(width=0.2, size=0.5)+     
   scale_fill_manual(values=col2)+   ##ylim(-0.1, 0.5)+
   ylab("SCC")+ylim(-0.1, 0.4)+
   theme_bw()+
   theme(legend.position="none",
         axis.title.x=element_blank(),
         axis.title.y=element_text(size=10),
         axis.text.x=element_blank(),
         axis.text.y=element_text(size=10),
         axis.ticks.x=element_blank())
         
         ## axis.title.y=element_text(size=9),
         ## axis.text.x=element_text(angle=30, hjust=1, size=8),
         ## axis.text.y=element_text(size=9))

comb_plots <- plot_grid(p1, p2, nrow=2, align="v") ##, rel_heights=c(1, 0.9))
figfn <- paste(outdir, "Figure1.2_corr_gene_ATAC_box.png", sep="")
ggsave(figfn, comb_plots, device="png", width=340, height=580, units="px", dpi=120)





###########################################################################
### Scatter plots of LFC on gene expression and chromatin accessibility ###
###########################################################################


rm(list=ls())
outdir <- "./Plots_pub/3_compare_plots.outs/"
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)


###
### plot data

fn <- paste(outdir, "1.2_LFC_gene_ATAC.comb.rds", sep="")
df0 <- read_rds(fn)%>%drop_na(LFC_gene, LFC_ATAC) 


###
### First example
ii <- "6_Monocyte_vitD"
df2 <- df0%>%dplyr::filter(comb==ii)


###
### if it is DEG or not 
fn <- "./sc_multiome_data/2_Differential/1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn_sampleID+new_treat.rds"
resDG <- read_rds(fn)%>%as.data.frame()%>%
    mutate(comb=paste(MCls, contrast, sep="_"),
           is_sig=ifelse(p.adjusted<0.1&abs(estimate)>0.5, 1, 0),
           is_sig=ifelse(is.na(is_sig), 0, is_sig))



resDG2 <- resDG%>%filter(comb==ii)%>%dplyr::select(gene, is_sig)

### combine plot data with 
df2 <- df2%>%left_join(resDG2, by="gene")


###
### add annotation text
corr <- cor.test(df2$LFC_ATAC, df2$LFC_gene, method="spearman")
rr <- round(as.numeric(corr$estimate), digits=3)
pval <- corr$p.value

symb <- case_when(pval<=1e-03~"***",
          pval>1e-03&pval<=0.01~"**",
          pval>1e-01&pval<=0.05~"*",
          TRUE~"NS")
eq <- deparse(bquote(italic(rho)==.(rr)~.(symb)))


### modification 
df3 <- df2%>%mutate(LFC_ATAC2=ifelse(abs(LFC_ATAC)>5, 5*sign(LFC_ATAC), LFC_ATAC),
                    LFC_gene2=ifelse(abs(LFC_gene)>5, 5*sign(LFC_ATAC), LFC_gene),
                    gr=ifelse(abs(LFC_ATAC)>5|abs(LFC_gene)>5, "gr2", "gr1")
                    )%>%as.data.frame()

y0 <- 5
x0 <- -3.7
  
p0 <- ggplot(df3, aes(x=LFC_ATAC2, y=LFC_gene2))+
   geom_point(aes(color=factor(is_sig),
                  size=factor(gr), shape=factor(gr))   )+
   scale_color_manual(values=c("grey", "red"), labels=c("Not DEG", "DEG"),
      guide=guide_legend(override.aes=list(size=1.5)))+
   scale_size_manual(values=c("gr1"=0.8, "gr2"=1.5), guide="none")+
   scale_shape_manual(values=c("gr1"=16, "gr2"=21), guide="none")+ 
   annotate("text", label=eq, x=x0, y=y0, size=3.5, parse=T)+ 
   scale_x_continuous("LFC on gene chromatin accessibility", limits=c(-6,6))+
   scale_y_continuous("LFC on gene expression", limits=c(-6,6))+
   geom_vline(xintercept=0, linetype="dashed", size=0.4)+
   geom_hline(yintercept=0, linetype="dashed", size=0.4)+ 
   ## geom_smooth(method="lm", formula=y~x, size=0.5, se=F)+ 
   ggtitle("vitamin D in Monocyte")+
   theme_bw()+ 
   theme(legend.position=c(0.2, 0.75),
         legend.title=element_blank(),
         legend.text=element_text(size=10),
         legend.key.size=unit(0.4, "cm"),
         legend.background=element_blank(),
         legend.key=element_blank(),
         legend.box.background=element_blank(),
         axis.title=element_text(size=10),
         axis.text=element_text(size=10),
         plot.title=element_text(hjust=0.5, size=12))        
 
figfn <- paste(outdir, "Figure1.3_", ii, "_scatter_LFC.png", sep="")
ggsave(figfn, p0, device="png", width=380, height=350, units="px", dpi=120)
 

###
### another example

###
### plot data, combine the information about if it is DEG or not

ii <- "6_Monocyte_zinc"
df2 <- df0%>%dplyr::filter(comb==ii)
resDG2 <- resDG%>%filter(comb==ii)%>%dplyr::select(gene, is_sig)

df2 <- df2%>%left_join(resDG2, by="gene")


## df3 <- df3%>%mutate(gr2=ifelse(grepl("^MT", gene), 1, 0))

###
### add annotation text
corr <- cor.test(df2$LFC_ATAC, df2$LFC_gene, method="spearman")
rr <- round(as.numeric(corr$estimate), digits=3)
pval <- corr$p.value

symb <- case_when(pval<=1e-03~"***",
          pval>1e-03&pval<=0.01~"**",
          pval>1e-01&pval<=0.05~"*",
          TRUE~"NS")
eq <- deparse(bquote(italic(rho)==.(rr)~.(symb)))

y0 <- 4.5
x0 <- -4

 
### modification plot data
df3 <- df2%>%mutate(LFC_ATAC2=ifelse(abs(LFC_ATAC)>5, 5*sign(LFC_ATAC), LFC_ATAC),
                    LFC_gene2=ifelse(abs(LFC_gene)>5, 5*sign(LFC_ATAC), LFC_gene),
                    gr=ifelse(abs(LFC_ATAC)>5|abs(LFC_gene)>5, "gr2", "gr1")
                    )%>%as.data.frame()


###
### annotate genes 
anno2 <- df3%>%filter(LFC_gene>4, is_sig>0, grepl("^MT", gene))

###
### main figures
p1 <- ggplot(df3, aes(x=LFC_ATAC2, y=LFC_gene2))+
   geom_point(aes(color=factor(is_sig),
                  size=factor(gr), shape=factor(gr))   )+
   scale_color_manual(values=c("grey", "red"), labels=c("Not DEG", "DEG"),
      guide=guide_legend(override.aes=list(size=1.5)))+
   scale_size_manual(values=c("gr1"=0.8, "gr2"=1.5), guide="none")+
   scale_shape_manual(values=c("gr1"=16, "gr2"=21), guide="none")+        
   annotate("text", label=eq, x=x0, y=y0, size=3.5, parse=T)+
   geom_text_repel(data=anno2,
       aes(x=LFC_ATAC2, y=LFC_gene2, label=gene),
       max.overlaps=30, box.padding=0.3,
       point.padding=0.1, nudge_x=0.01, nudge_y=0.01,
       min.segment.length=0, segment.curvature=-1e-20, segment.size=unit(0.3, "mm"),
       arrow=arrow(length=unit(0.015, "npc")),
       color="maroon3", fontface="italic", size=2.8)+ 
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

figfn <- paste(outdir, "Figure1.3_", ii, "_scatter_LFC.png", sep="")
ggsave(figfn, p1, device="png", width=380, height=350, units="px", dpi=120)







####################################################################
#### Enrichment analysis to examine if DEGs are enriched in DARs ###
####################################################################

rm(list=ls()) 

outdir <- "./Plots_pub/3_compare_plots.outs/"


### DEG results
fn <- "./sc_multiome_data/2_Differential/1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn_sampleID+new_treat.rds"
resDG <- read_rds(fn)%>%as.data.frame()%>%mutate(comb=paste(MCls, contrast, sep="_"))

gene0 <- unique(resDG$gene)
gene1 <- grch38%>%dplyr::filter(chr%in%as.character(1:22))%>%pull(symbol)%>%unique()
gene2 <- intersect(gene0, gene1)

resDG <- resDG%>%filter(gene%in%gene2)
 
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

###
### next time skip the above enrich procedure and directly read data.

rm(list=ls()) 

outdir <- "./Plots_pub/3_compare_plots.outs/"

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

ylab2 <- plotDF$contrast
names(ylab2) <- plotDF$comb

p <- ggplot(plotDF, aes(x=log_odds, y=comb))+
   geom_errorbarh(aes(xmax=log_upper, xmin=log_lower, colour=MCls),
       size=0.5, height=0.2)+
   geom_point(aes(colour=MCls), shape=19, size=1.5)+
   scale_colour_manual(values=col1)+
   geom_vline(aes(xintercept=0), size=0.25, linetype="dashed")+ 
   xlab("log odds ratio")+xlim(-3, 6)+
   scale_y_discrete(labels=ylab2)+ 
   theme_bw()+
   theme(##plot.title=element_text(hjust=0.5, size=14),
         axis.title.y=element_blank(),
         axis.title.x=element_text(size=10),
         axis.text.x=element_text(size=10),
         axis.text.y=element_text(size=10),
         legend.position="none")
         ## legend.title=element_blank(),
         ## legend.text=element_text(size=8),
         ## legend.key.size=unit(0.4, "cm"))
         ##legend.position="none")
 
figfn <- paste(outdir, "Figure2_enrich_forest.png", sep="")
ggsave(figfn, p, device="png", width=380, height=580, units="px", dpi=120)




############################
#### supplementary  
############################



################
### FigS3_1
################

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

figfn <- paste(outdir, "FigS3_1_corr_heatmap.png", sep="")
###ggsave(figfn, p, width=520, height=520, units="px",dpi=120)
png(figfn, height=520, width=520, res=120)
print(p)
dev.off()    






###############
### FigS3_2 
###############


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
gene0 <- unique(res$gene)
gene1 <- grch38%>%dplyr::filter(chr%in%as.character(1:22))%>%pull(symbol)%>%unique()
gene2 <- intersect(gene0, gene1)

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
figfn <- paste(outdir, "FigS3_2_scatter_SCCvsDEGs_DARs.png", sep="")
ggsave(figfn, comb_plots, device="png", width=800, height=380, units="px", dpi=120)









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





