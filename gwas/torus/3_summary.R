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



### color setting
col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
  "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
   "6_Monocyte"="#984ea3", "7_dnT"="black")
col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4", "vitD"="seagreen4",
          "vitE"="salmon3", "zinc"="maroon3", "peaking"="#7570b3")



###
### column names
x <- fread("./torus_input/zzz_torus.annot.gz", header=T, data.table=F)
nsnp <- colSums(x[,-1])
df <- data.frame(comb=names(nsnp), nsnp=nsnp)

tmp2 <- df%>%dplyr::select(comb)


### forest plots for each trait
traits <- read.table("traits_of_interest.txt")$V1
traits <- sort(traits)

for (ii in traits){
 
cat(ii, "\n")    
###    
fn <- paste("./torus_output/", ii, ".est", sep="")
est <- read.table(fn)
est2 <- est[-1,]
## est2 <- est2%>%full_join(tmp2, by="V1")


###
est2 <- est2%>%mutate(comb=as.character(gsub("\\.1", "", V1)))%>%
   dplyr::select(comb, "odds"="V2", "CI_lower"="V3", "CI_upper"="V4")

est2 <- est2%>%mutate(se2=(abs((odds-CI_lower))/1.96)^2,
    shrink=1/(1+se2), bhat=odds*shrink, ss=sqrt(shrink*se2),
    b_lower=bhat-1.96*ss, b_upper=bhat+1.96*ss)

    
##
x <- str_split(est2$comb, "_", simplify=T)
 
est2 <- est2%>%mutate(MCls=paste(x[,1], x[,2], sep="_"),
    treats=x[,3], treats=ifelse(treats=="", "peaking", treats),
    MCls_val=ifelse(grepl("peaking", MCls), "8_peaking", MCls),
    comb2=fct_reorder(comb, MCls_val) )

## est3 <- est2%>%mutate(bhat=ifelse(b_lower<(-12)|b_upper>12, NA, bhat))
est3 <- est2%>%mutate(odds=ifelse(CI_lower<(-10)|CI_upper>10, NA, odds),
                      CI_lower=ifelse(is.na(odds), NA, CI_lower),
                      CI_upper=ifelse(is.na(odds), NA, CI_upper))

    
p2 <- ggplot(est3, aes(x=odds, y=comb2, color=factor(MCls)))+
    geom_errorbarh(aes(xmax=CI_upper, xmin=CI_lower), linewidth=0.5, height=0.2)+
    geom_point(shape=19, size=0.5)+
    geom_vline(aes(xintercept=0), size=0.25, linetype="dashed")+
    ##scale_y_discrete(labels=ylab)+
    scale_x_continuous("log odds ratio", breaks=seq(-10,10, by=2), limits=c(-11,11))+
   ##breaks=seq(-4, 4, by=2), limits=c(-4.5,4.5))+
    scale_colour_manual(values=col1)+
    ggtitle(ii)+
    theme_bw()+
    theme(legend.title=element_blank(),
          ##legend.key=element_blank(),
          legend.key.size=grid::unit(0.8, "lines"),
          legend.text=element_text(size=10),
          legend.position="none",
          plot.title=element_text(hjust=0.5, size=10),
          axis.title.x=element_text(size=12),
          axis.title.y=element_blank(),
          axis.text.x=element_text(size=10),
          axis.text.y=element_text(size=8))
          ##axis.title.x=element_text(size=10),
          ###strip.text=element_text(size=12),
          ##plot.title=element_text(hjust=0.5, size=12),
          ## axis.text.x=element_text(size=9, angle=-90, hjust=0, vjust=0.5),
          ###axis.text=element_text(size=9))

figfn <- paste(outdir, "Figure1_", ii, "_forest.png", sep="")
ggsave(figfn, plot=p2, device="png", width=550, height=650, units="px", dpi=120)     
## png(figfn, width=580, height=650, res=120)
## print(p2)
## dev.off()

###
}    
 
    
###
###


####################
### summary SNPs ###
####################

outdir <- "./3_summary.outs/"
     
fn <- "/nfs/rprdata/julong/sc_multiome/analysis_multiomic_2023-03-27/3_motif/Response_motif/2_summary.xlsx"
summ <- read.xlsx(fn)
summ <- summ[,c(1:3,6:7)]


x <- fread("./torus_input/zzz_torus.annot.gz", header=T, data.table=F)
nsnp <- colSums(x[,-c(1,2)])
df <- data.frame(comb=names(nsnp), nsnp_0.1=nsnp)

###
## x2 <- fread("./analysis2_th0.2/torus_input/zzz_torus.annot.gz", header=T, data.table=F)
## nsnp <- colSums(x2[,-c(1,42)])
## df2 <- data.frame(comb=names(nsnp), nsnp_0.2=nsnp)

### 
## df3 <- df%>%full_join(df2, by="comb")%>%mutate(comb=gsub("_d$", "", comb))
df3 <- df%>%mutate(comb=gsub("_d$", "", comb))%>%left_join(summ, by="comb")%>%
    mutate(nsnp_0.1=ifelse(is.na(nsnp_0.1), 0, nsnp_0.1),
           th_0.1=round(th_0.1, 3))%>%
           ## th_0.2=round(th_0.2, 3))%>%
    arrange(comb)

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



