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


## fn <- "../../../analysis_multiomic_2023-03-27/3_motif/Response_motif/2_summary.xlsx"
## summ <- read.xlsx(fn)
## summ <- summ[,c(1, 6, 7)]
## tmp <- data.frame(comb="peaking", MCls="peaking", treats="peaking")
## summ <- rbind(summ, tmp)

comb <- read.table("annot_ls.txt")$V1
cvt <- str_split(comb, "_", simplify=T)

summ <- data.frame("comb"=comb, MCls=paste(cvt[,1], cvt[,2], sep="_"), treats=cvt[,3])
## summ2 <- data.frame("comb"=c("Peak", "Response"),
##                     MCls=c("88_Peak", "88_Response"), treats=c("88_Peak", "88_Response"))
## summ <- rbind(summ, summ2)


###
est2 <- map_dfr(comb, function(ii){
   ##
   fn <- paste("./torus_output/", ii, "_WBL.est", sep="")
   x <- read.table(fn)
   ##kk <- ifelse(grepl("peaking", ii), 2, 3)
   kk <- 3
   x2 <- x[kk,]
   x2$comb <- ii
   print(ii) 
   x2
})

### Union response
## fn <- "../torus_output/Union_WBL.est"
## x <- read.table(fn)
## x2 <- x[2:3,]
## x2$comb <- c("Peak", "Response")
## est2 <- rbind(est2, x2)
 
est2 <- est2%>%dplyr::rename("odds"="V2", "CI_lower"="V3", "CI_upper"="V4")
plotDF <- summ%>%full_join(est2, by="comb")%>%arrange(MCls)



## est2 <- est2%>%mutate(MCls=paste(cvt[,1], cvt[,2], sep="_"),
##     treats=cvt[,3], treats=ifelse(treats=="", "peaking", treats),
##     treat_val=ifelse(grepl("peaking", treats), "zzz_peak", treats),
##     category2=fct_reorder(category, treat_val) )

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
 
p2 <- ggplot(plotDF, aes(x=odds, y=comb, color=factor(MCls)))+
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
 
figfn <- paste(outdir, "Figure1.3_multiCondition_forest.png", sep="")
png(figfn, width=700, height=650, res=120)
p2
dev.off()






## figfn <- paste(outdir, "Figure1.1_multiCondition_forest.png", sep="")
## png(figfn, width=600, height=650, res=120)
## p2
## dev.off()


############################
#### compare old results ###
############################









###################
### summary snp ###
###################

fn <- "./torus_input/zzz_torus.annot.gz"
x <- fread(fn, header=T, data.table=F)
nsnp <- colSums(x[,-1])
df1 <- data.frame(comb=names(nsnp), nsnp=nsnp)

##
###
fn <- "./analysis2_each/torus_input/zzz_torus.annot.gz"
x2 <- fread(fn, header=T, data.table=F)
nsnp <- colSums(x2[,-1])
df2 <- data.frame(comb=names(nsnp), nsnp2=nsnp)
 
df3 <- df1%>%full_join(df2, by="comb")%>%arrange(comb)%>%
    mutate(nsnp=ifelse(is.na(nsnp), 0, nsnp))

opfn <- "./3_summary.outs/0_summ.xlsx"
write.xlsx(df3, opfn, overwrite=T)







