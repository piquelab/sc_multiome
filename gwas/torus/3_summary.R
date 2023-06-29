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
 
outdir <- "./3_summary.outs/no_shrinkage/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)

###
### column names
fn <- "./analysis2_th0.2/torus_output/Asthma.est"
tmp <- read.table(fn)[-1,]
tmp2 <- tmp%>%dplyr::select(V1)


### forest plots for each trait
traits <- read.table("traits_ls.txt")$V1
for (ii in traits){
 
cat(ii, "\n")    
###    
fn <- paste("./torus_output/", ii, ".est", sep="")
est <- read.table(fn)
est2 <- est[-1,]
est2 <- est2%>%full_join(tmp2, by="V1")


###
est2 <- est2%>%mutate(category=as.character(gsub("\\.1", "", V1)))%>%
   dplyr::rename("comb"="V1", "odds"="V2", "CI_lower"="V3", "CI_upper"="V4")

est2 <- est2%>%mutate(se2=(abs((odds-CI_lower))/1.96)^2,
    shrink=1/(1+se2), bhat=odds*shrink, ss=sqrt(shrink*se2),
    b_lower=bhat-1.96*ss, b_upper=bhat+1.96*ss)

    
##
cvt <- str_split(est2$category, "_", simplify=T)
 
est2 <- est2%>%mutate(MCls=paste(cvt[,1], cvt[,2], sep="_"),
    treats=cvt[,3], treats=ifelse(treats=="", "peaking", treats),
    treat_val=ifelse(grepl("peaking", treats), "zzz_peaking", treats),
    category2=fct_reorder(category, treat_val) )

###
col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4", "vitD"="seagreen4",
          "vitE"="salmon3", "zinc"="maroon3", "peaking"="#7570b3")
est3 <- est2%>%mutate(bhat=ifelse(b_lower<(-17)|b_upper>12, NA, bhat))
    
p2 <- ggplot(est2, aes(x=odds, y=category2, color=factor(treats)))+
    geom_errorbarh(aes(xmax=CI_upper, xmin=CI_lower), size=0.5, height=0.2)+
    geom_point(shape=19, size=0.5)+
    geom_vline(aes(xintercept=0), size=0.25, linetype="dashed")+
    ##scale_y_discrete(labels=ylab)+
    scale_x_continuous("log odds ratio", breaks=seq(-16,8, by=4), limits=c(-17,8))+
   ##breaks=seq(-4, 4, by=2), limits=c(-4.5,4.5))+
    scale_colour_manual(values=col2)+
    ggtitle(ii)+
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

 
figfn <- paste(outdir, "Figure1_", ii, "_forest.png", sep="")
png(figfn, width=600, height=650, res=120)
print(p2)
dev.off()
###
}    

###
###


####################
### summary SNPs ###
####################

fn <- "/nfs/rprdata/julong/sc_multiome/analysis_multiomic_2023-03-27/3_motif/Response_motif/2_summary.xlsx"
summ <- read.xlsx(fn)
summ <- summ[,c(1:5,7)]



x <- fread("./torus_input/zzz_torus.annot.gz", header=T, data.table=F)
nsnp <- colSums(x[,-c(1,2)])
df <- data.frame(comb=names(nsnp), nsnp_0.1=nsnp)

###
x2 <- fread("./analysis2_th0.2/torus_input/zzz_torus.annot.gz", header=T, data.table=F)
nsnp <- colSums(x2[,-c(1,42)])
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


#####################
### check results ###
#####################

traits <- read.table("traits_ls.txt")$V1

outdir <- "./3_summary.outs/check/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


for (trait in traits){

###    
fn <- paste("./torus_output/", trait, ".est", sep="")
x <- read.table(fn)
x <- x[2:41,]
x <- x%>%mutate(category=gsub("\\.1", "", V1))%>%
    dplyr::select(odds=V2, CI_lower=V3, CI_upper=V4, category)

fn <- paste("./analysis3_each/torus_output/", trait, "_comb.est", sep="")
x2 <- read.table(fn, header=T)
names(x2) <- c(paste0("ind_", names(x2)[1:3]), "category")

xcomb <- x%>%inner_join(x2, by="category")
xcomb <- xcomb%>%
    mutate(treats=gsub(".*_", "", category),
           treat2=ifelse(grepl("peaking", treats), "zzz_peaking", treats))%>%arrange(treat2)    

xcomb <- xcomb[,c(4, 1:3, 5:8)] 
    
opfn <- paste(outdir, trait, "_check.xlsx", sep="")
write.xlsx(xcomb, file=opfn, overwrite=T)

cat(trait, "\n")
}



