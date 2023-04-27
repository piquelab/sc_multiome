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
 
outdir <- "./3_summary.outs/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


################################################################################################
### 1. combine annotation, peaks, cell-type motifs, treatment motifs in different cell-types ###
################################################################################################

fn <- "./analysis2_th0.2/torus_output/Whole_Blood.est"
tmp <- read.table(fn)[2:42,]
tmp2 <- tmp%>%dplyr::select(V1)

###
###
fn <- "./torus_output/Whole_Blood.est"
est <- read.table(fn)[2:41,]
est2 <- est%>%full_join(tmp2, by="V1")

est2 <- est2%>%arrange(V1)%>%
   mutate(order=as.numeric(1:nrow(est2)),
   category=as.character(gsub(".1", "", V1)))%>%
   dplyr::rename("odds"="V2", "CI_lower"="V3", "CI_upper"="V4")
##
cvt <- str_split(est2$category, "_", simplify=T)

est2 <- est2%>%mutate(MCls=paste(cvt[,1], cvt[,2], sep="_"),
    treats=cvt[,3], treats=ifelse(treats=="", "peaking", treats),
    treat_val=ifelse(grepl("peaking", treats), "zzz_peak", treats),
    category2=fct_reorder(category, treat_val) )

###
col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4", "vitD"="seagreen4",
          "vitE"="salmon3", "zinc"="maroon3", "peaking"="#7570b3")
 
p2 <- ggplot(est2, aes(x=odds, y=category2, color=factor(treats)))+
    geom_errorbarh(aes(xmax=CI_upper, xmin=CI_lower), size=0.5, height=0.2)+
    geom_point(shape=19, size=0.5)+
    geom_vline(aes(xintercept=0), size=0.25, linetype="dashed")+
    ##scale_y_discrete(labels=ylab)+
    scale_x_continuous("log odds ratio", breaks=seq(-2, 3, by=1), limits=c(-2, 3))+
    scale_colour_manual(values=col2)+
    ggtitle("GTEx WBL")+
    theme_bw()+
    theme(legend.title=element_blank(),
          ##legend.key=element_blank(),
          legend.key.size=grid::unit(0.8, "lines"),
          legend.text=element_text(size=10),
          plot.title=element_text(hjust=0.5, size=12),
          axis.title.x=element_text(size=12),
          axis.title.y=element_blank(),
          axis.text=element_text(size=8))
          ##axis.title.x=element_text(size=10),
          ###strip.text=element_text(size=12),
          ##plot.title=element_text(hjust=0.5, size=12),
          ## axis.text.x=element_text(size=9, angle=-90, hjust=0, vjust=0.5),
          ###axis.text=element_text(size=9))
 
figfn <- paste(outdir, "Figure1.1_multiCondition_forest.png", sep="")
png(figfn, width=600, height=650, res=120)
p2
dev.off()




## figfn <- paste(outdir, "Figure1.1_multiCondition_forest.png", sep="")
## png(figfn, width=600, height=650, res=120)
## p2
## dev.off()




###################
### summary snp ###
###################
 
fn <- "/nfs/rprdata/julong/sc_multiome/analysis_multiomic_2023-03-27/3_motif/Response_motif/2_summary.xlsx"
summ <- read.xlsx(fn)
summ <- summ[,c(1:5,7)]


x <- fread("./torus_input/zzz_torus.annot.gz", header=T, data.table=F)
nsnp <- colSums(x[,-c(1,2)])
df <- data.frame(comb=names(nsnp), nsnp_0.1=nsnp)

###
x2 <- fread("./analysis2_th0.2/torus_input/zzz_torus.annot.gz", header=T, data.table=F)
nsnp <- colSums(x2[,-c(1,2)])
df2 <- data.frame(comb=names(nsnp), nsnp_0.2=nsnp)

###
df3 <- df%>%full_join(df2, by="comb")%>%mutate(comb=gsub("_d$", "", comb))
df3 <- df3%>%left_join(summ, by="comb")%>%
    mutate(nsnp_0.1=ifelse(is.na(nsnp_0.1), 0, nsnp_0.1),
           th_0.1=round(th_0.1, 3),
           th_0.2=round(th_0.2, 3))%>%
    arrange(treats)%>%
    dplyr::select(comb, th_0.1, nmotif_0.1, nsnp_0.1, th_0.2, nmotif_0.2, nsnp_0.2, treats)

###
opfn <- "./3_summary.outs/0_summ.xlsx"
write.xlsx(df3, opfn, overwrite=T)







