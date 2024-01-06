###
###
library(Matrix)
library(tidyverse)
library(data.table)
library(DESeq2)
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

outdir <- "./Plots_pub/2_diff_plots/"
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)
 
 
 
###
### script for supp plots in the manuscript  



#########################################################
### MA plots of differential gene expression analysis ### 
#########################################################

rm(list=ls())

outdir <- "./Plots_pub/2_diff_plots/"
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)

###
### data for MA plots
fn <- "./sc_multiome_data/2_Differential/1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn_sampleID+new_treat.rds"
res <- read_rds(fn)%>%as.data.frame()%>%
    mutate(comb=paste(MCls, contrast, sep="_"), color=ifelse(p.adjusted<0.1, T, F),
           color2=ifelse(is.na(color), F, color))%>% 
    drop_na(estimate)
    
MCls <- sort(unique(res$MCls))
treats <- sort(unique(res$contrast))

 
figfn <- paste(outdir, "FigS2_1_differential_gene.MA.png", sep="")
png(figfn, width=1200, height=1600, res=120)

par(mar=c(4,4,2,2),mgp=c(2,1,0))
x <- matrix(1:48, 8, 6, byrow=T)
layout(x)
 

for ( oneMCl in MCls){

##1
   for ( ii in treats){
      ##
       cat(oneMCl, ii, "\n")
       ##
       res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast==ii)%>%
           dplyr::select(baseMean, estimate, color=color2, p.value, p.adjusted)
        print(plotMA(res2[,1:3], colLine="NA", main=ii, cex.main=1.8, font.main=1, cex.axis=1, cex.lab=1.2))
   }       
   print(mtext(oneMCl, side=4, line=0.8, cex=1, font=1)) 
}

dev.off()
 
### End-plot 







###
### qq plots of differential gene expression analysis 
 

###
### setting colors
col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
  "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
   "6_Monocyte"="#984ea3", "7_dnT"="black")
col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
       "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3")


 
fn <- "./sc_multiome_data/2_Differential/1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn_sampleID+new_treat.rds"
res <- read_rds(fn)%>%as.data.frame()%>%mutate(comb=paste(MCls, contrast, sep="_"))%>%drop_na(estimate)


comb <- sort(unique(res$comb))
dfNew <- map_dfr(comb, function(ii){
  res2 <- res%>%dplyr::filter(comb==ii)
  ngene <- nrow(res2)
  res2 <- res2%>%
     arrange(p.value)%>%
     mutate(observed=-log10(p.value), expected=-log10(ppoints(ngene)))
  res2
})



###
p1 <- ggplot(dfNew, aes(x=expected, y=observed, color=contrast))+
   ggrastr::rasterise(geom_point(size=0.7),dpi=300)+    
   geom_abline(colour="black")+
   scale_color_manual(values=col2,
      guide=guide_legend(override.aes=list(size=3)))+
   facet_wrap(~MCls, scales="free", ncol=3)+
   xlab(bquote("Expected"~ -log[10]~"("~italic(plain(P))~")"))+
   ylab(bquote("Observed"~-log[10]~"("~italic(plain(P))~")"))+
   ###scale_y_continuous(expand=expansion(mult=c(0,0.3)))+
   theme_bw()+
   theme(legend.position=c(0.85, 0.18),
         legend.title=element_blank(),
         legend.background=element_blank(),
         legend.box.background=element_blank(),
         legend.key.size=grid::unit(0.6, "cm"),
         legend.text=element_text(size=12),
         axis.title=element_text(size=12),
         axis.text=element_text(size=12),
         strip.text.x=element_text(size=14))
###
###
figfn <- paste(outdir, "FigS2_2_differential_gene.qq.pdf", sep="")
ggsave(figfn, p1, width=8, height=8)




###########################################
### Differential accessibility analysis ### 
###########################################


rm(list=ls())

outdir <- "./Plots_pub/2_diff_plots/"
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)

###
### data for MA plots
fn <- "./sc_multiome_data/2_Differential/1.2_DiffPeak.outs/2.0_DESeq.results_0.1_cn_sampleID+new_treat.rds"

res <- read_rds(fn)%>%as.data.frame()%>%
    mutate(comb=paste(MCls, contrast, sep="_"), color=ifelse(p.adjusted<0.1, T, F),
           color2=ifelse(is.na(color), F, color))%>% 
    drop_na(estimate)
    
MCls <- sort(unique(res$MCls))
treats <- sort(unique(res$contrast))


figfn <- paste(outdir, "FigS2_6_differential_peak.MA.png", sep="")
png(figfn, width=1200, height=1600, res=120)

par(mar=c(4,4,2,2),mgp=c(2,1,0))
x <- matrix(1:48, 8, 6, byrow=T)
layout(x)
 

for ( oneMCl in MCls){

##1
   for ( ii in treats){
      ##
       cat(oneMCl, ii, "\n")
       ##
       res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast==ii)%>%
           dplyr::select(baseMean, estimate, color=color2, p.value, p.adjusted)
        print(plotMA(res2[,1:3], colLine="NA", main=ii, cex.main=1.8, font.main=1, cex.axis=1, cex.lab=1.2))
   }       
   print(mtext(oneMCl, side=4, line=0.8, cex=1, font=1)) 
}

dev.off()
 
### End-plot 







###
### qq plots of differential accessibility analysis

fn <- "./sc_multiome_data/2_Differential/1.2_DiffPeak.outs/2.0_DESeq.results_0.1_cn_sampleID+new_treat.rds"
res <- read_rds(fn)%>%as.data.frame()%>%mutate(comb=paste(MCls, contrast, sep="_"))%>%drop_na(estimate)


### 
comb <- sort(unique(res$comb))
dfNew <- map_dfr(comb, function(ii){
  res2 <- res%>%dplyr::filter(comb==ii)
  ngene <- nrow(res2)
  res2 <- res2%>%
     arrange(p.value)%>%
     mutate(observed=-log10(p.value), expected=-log10(ppoints(ngene)))
  res2
})



###
p2 <- ggplot(dfNew, aes(x=expected, y=observed, color=contrast))+
   ggrastr::rasterise(geom_point(size=0.7),dpi=300)+    
   geom_abline(colour="black")+
   scale_color_manual(values=col2,
      guide=guide_legend(override.aes=list(size=3)))+
   facet_wrap(~MCls, scales="free", ncol=3)+
   xlab(bquote("Expected"~ -log[10]~"("~italic(plain(P))~")"))+
   ylab(bquote("Observed"~-log[10]~"("~italic(plain(P))~")"))+
   ###scale_y_continuous(expand=expansion(mult=c(0,0.3)))+
   theme_bw()+
   theme(legend.position=c(0.85, 0.18),
         legend.title=element_blank(),
         legend.background=element_blank(),
         legend.box.background=element_blank(),
         legend.key.size=grid::unit(0.6, "cm"),
         legend.text=element_text(size=12),
         axis.title=element_text(size=12),
         axis.text=element_text(size=12),
         strip.text.x=element_text(size=14))
###
###
figfn <- paste(outdir, "FigS2_7_differential_peak.qq.pdf", sep="")
ggsave(figfn, p2, width=8, height=8)







