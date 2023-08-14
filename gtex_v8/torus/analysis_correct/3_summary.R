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


## fn <- "/nfs/rprdata/julong/sc_multiome/analysis_multiome_JW_2023-06-03/4_motif_plots.outs/response_motif_correct/2.2_response_motif_th0.1.txt"
## resp <- read.table(fn, header=T)
## df <- data.frame(comb=sort(unique(resp$comb)))


###
###
fn <- "./torus_output/Whole_Blood.est"
est <- read.table(fn)
est2 <- est[2:31,]

est2 <- est2%>%arrange(V1)%>%
   mutate(comb=as.character(gsub(".1$", "", V1)))%>%
   dplyr::rename("odds"="V2", "CI_lower"="V3", "CI_upper"="V4")


##
cvt <- str_split(est2$comb, "_", simplify=T)
est2 <- est2%>%mutate(MCls=paste(cvt[,1], cvt[,2], sep="_"),
    MCls_val=ifelse(grepl("peaking", MCls), "9_peak", MCls),
    comb2=fct_reorder(comb, MCls_val) )




###
### color cell type
col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
  "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
   "6_Monocyte"="#984ea3", "7_dnT"="black", "peaking"="#828282")

###
### color treats
col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4", "vitD"="seagreen4",
          "vitE"="salmon3", "zinc"="maroon3", "peaking"="#7570b3")
 
p <- ggplot(est2, aes(x=odds, y=comb, color=factor(MCls)))+
    geom_errorbarh(aes(xmax=CI_upper, xmin=CI_lower), size=0.5, height=0.2)+
    geom_point(shape=19, size=0.5)+
    geom_vline(aes(xintercept=0), size=0.25, linetype="dashed")+
    ##scale_y_discrete(labels=ylab)+
    scale_x_continuous("log odds ratio", breaks=seq(-2, 3, by=1), limits=c(-2, 3))+
    scale_colour_manual(values=col1)+
    ggtitle("GTEx WBL")+
    theme_bw()+
    theme(legend.title=element_blank(),
          ##legend.key=element_blank(),
          legend.key.size=grid::unit(0.8, "lines"),
          legend.text=element_text(size=10),
          plot.title=element_text(hjust=0.5, size=12),
          axis.title.x=element_text(size=12),
          axis.title.y=element_blank(),
          axis.text=element_text(size=9))
          ##axis.title.x=element_text(size=10),
          ###strip.text=element_text(size=12),
          ##plot.title=element_text(hjust=0.5, size=12),
          ## axis.text.x=element_text(size=9, angle=-90, hjust=0, vjust=0.5),
          ###axis.text=element_text(size=9))
 
figfn <- paste(outdir, "Figure1_multiCondition_forest.png", sep="")
ggsave(figfn, p, device="png", width=700, height=650, units="px", dpi=120)




###
### Forest union response motif

fn <- "./torus_output/Union_WBL.est"
est <- read.table(fn)
est2 <- est[2:3,]
names(est2)[2:4] <- c("odds", "CI_lower", "CI_upper")
est2$comb <- c("Peak", "Response motifs")


###
p2 <- ggplot(est2, aes(x=odds, y=comb, color=factor(comb)))+
    geom_errorbarh(aes(xmax=CI_upper, xmin=CI_lower), size=0.35,  height=0.1)+
    geom_point(shape=19, size=0.35)+
    geom_vline(aes(xintercept=0), size=0.25, linetype="dashed")+
    ##scale_y_discrete(labels=ylab)+
    scale_x_continuous("log odds ratio", breaks=seq(-2, 3, by=1), limits=c(-2, 3))+
    scale_colour_manual(values=c("Peak"="#7570b3", "Response motifs"="#e7298a"))+
    ggtitle("GTEx WBL")+
    theme_bw()+
    theme(legend.position="none",
          plot.title=element_text(hjust=0.5, size=12),
          axis.title.x=element_text(size=12),
          axis.title.y=element_blank(),
          axis.text=element_text(size=10))
          ##axis.title.x=element_text(size=10),
          ###strip.text=element_text(size=12),
          ##plot.title=element_text(hjust=0.5, size=12),
          ## axis.text.x=element_text(size=9, angle=-90, hjust=0, vjust=0.5),
          ###axis.text=element_text(size=9))
  
figfn <- paste(outdir, "Figure1.2_union_forest.png", sep="")
ggsave(figfn, p2, device="png", width=350, height=420, units="px", dpi=120)




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


#########################
#### sumary for final ###
#########################

fn <- "../../analysis_multiomic_2023-03-27/3_motif/Response_motif/2_summary.xlsx"
summ <- read.xlsx(fn)
summ <- summ[,c(1:3, 6, 7)]


x <- fread("./torus_input/zzz_torus.annot.gz", header=T, data.table=F)
nsnp <- colSums(x[,-c(1,2)])
df <- data.frame(comb=gsub("_d$", "", names(nsnp)), nsnp_0.1=nsnp)

df2 <- summ%>%full_join(df, by="comb")%>%arrange(comb)
df2 <- df2%>%mutate(th_0.1=round(th_0.1, 2))
opfn <- "./3_summary.outs/0.1_summ.xlsx"
write.xlsx(df2, file=opfn, overwrite=T)
