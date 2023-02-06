###
library(tidyverse)
## library(parallel)
library(data.table)
## library(purrr)
library(GenomicRanges)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(Signac)
library(SeuratWrappers)
library(SeuratObject)
## library(cicero, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
## library(monocle3)
###
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(viridis)
theme_set(theme_grey())

outdir <- "2_demux.outs"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)



###############################################################
### read results of demuxlet
###############################################################
basefolder <- "/wsu/home/groups/piquelab/scGxE/counts_cellranger_arc/demuxlet-RNA/demuxlet"

list.files(basefolder, "^scGxE_1-*")

expNames <- unique(gsub("\\..*", "", list.files(basefolder, "^scGxE_1-*")))
expNames

folders <- paste(basefolder, "/", expNames, ".out.best", sep="")
folders

names(folders) <- expNames
folders


###
demux <- map_dfr(expNames, function(ii){
  fn <- folders[ii]
  dd <- fread(fn)
  dd <- dd%>%mutate(BARCODE=gsub("-1","",BARCODE),
      NEW_BARCODE=paste(ii, "_", BARCODE, sep=""),
      EXP=ii, treats=gsub("scGxE_1-[1-8]", "", EXP))
  dd
})

head(demux)

demux <- demux%>%mutate(NEW_BARCODE=gsub("-1", "", NEW_BARCODE))
table(demux$EXP)

library(dplyr)
demux <- demux %>% mutate(treats = case_when(
                     endsWith(EXP, "1") ~ "caffeine",
                     endsWith(EXP, "2") ~ "vitE",
                     endsWith(EXP, "3") ~ "zinc",
                     endsWith(EXP, "4") ~ "vitA",
                     endsWith(EXP, "5") ~ "vitD",
                     endsWith(EXP, "6") ~ "nicotine",
                     endsWith(EXP, "7") ~ "etOH",
                     endsWith(EXP, "8") ~ "water"
                     ))

head(demux)

### output
opfn <- "./2_demux.outs/1_demux.rds" ## 115,812 barcodes
write_rds(demux, opfn)







#########################################
### subset SNG cells of seurat object ###
#########################################
rm(list=ls())

atac <- read_rds("./1_Merge.outs/1_seurat.merge.rds")

colnames(atac@meta.data)
head(rownames(atac@meta.data))
head(Cells(atac))

atac <- RenameCells(atac, new.names=gsub("-1","",Cells(atac)))
head(Cells(atac))

demux <- read_rds("./2_demux.outs/1_demux.rds")%>%
   dplyr::select(NEW_BARCODE, SNG.BEST.GUESS, DROPLET.TYPE)
head(demux)

colnames(atac@meta.data)

meta <- atac@meta.data%>%mutate("NEW_BARCODE"=gsub("-1", "", barcode))
colnames(meta)
head(rownames(meta))

meta <- meta%>%left_join(demux, by="NEW_BARCODE")
colnames(meta)
rownames(meta) <- meta$NEW_BARCODE

atac <- AddMetaData(atac, meta)
colnames(atac@meta.data)
head(rownames(atac@meta.data))
table(atac@meta.data$DROPLET.TYPE)
table(atac@meta.data$SNG.BEST.GUESS)

##
## atac2 <- subset(atac, subset=DROPLET.TYPE=="SNG")
## colnames(atac2@meta.data)
## table(atac2@meta.data$DROPLET.TYPE)
## table(atac2@meta.data$SNG.BEST.GUESS)



##
atac2 <- subset(atac, subset= (SNG.BEST.GUESS=="AL-002" | SNG.BEST.GUESS=="AL-018" | SNG.BEST.GUESS=="AL-026" | SNG.BEST.GUESS=="AL-030" | SNG.BEST.GUESS=="AL-041" | SNG.BEST.GUESS=="AL-044" | SNG.BEST.GUESS=="AL-045" | SNG.BEST.GUESS=="AL-056" | SNG.BEST.GUESS=="AL-058" | SNG.BEST.GUESS=="AL-062" | SNG.BEST.GUESS=="AL-070" | SNG.BEST.GUESS=="AL-105"))
table(atac2@meta.data$DROPLET.TYPE)
table(atac2@meta.data$SNG.BEST.GUESS)

atac2 <- subset(atac2, subset=DROPLET.TYPE=="SNG")
colnames(atac2@meta.data)
table(atac2@meta.data$DROPLET.TYPE)
table(atac2@meta.data$SNG.BEST.GUESS)


opfn <- "./2_demux.outs/2_seurat.merge.SNG.rds"
write_rds(atac2, opfn)




####################################
### summary results - demux file ###
####################################
demux <- read_rds("./2_demux.outs/1_demux.rds")

colnames(demux)

head(demux$EXP)

demux <- demux%>%mutate(BATCH=gsub("-.*","",EXP))                
## demux <- demux%>%dplyr::select("SNG.BEST.GUESS", "EXP", "BATCH", "DROPLET.TYPE")

table(demux$SNG.BEST.GUESS)

demux2 <- subset(demux, subset= (SNG.BEST.GUESS=="AL-002" | SNG.BEST.GUESS=="AL-018" | SNG.BEST.GUESS=="AL-026" | SNG.BEST.GUESS=="AL-030" | SNG.BEST.GUESS=="AL-041" | SNG.BEST.GUESS=="AL-044" | SNG.BEST.GUESS=="AL-045" | SNG.BEST.GUESS=="AL-056" | SNG.BEST.GUESS=="AL-058" | SNG.BEST.GUESS=="AL-062" | SNG.BEST.GUESS=="AL-070" | SNG.BEST.GUESS=="AL-105"))

table(demux2$SNG.BEST.GUESS)

table(demux2$DROPLET.TYPE)

##(1)
dd <- demux2%>%
   group_by(EXP, DROPLET.TYPE)%>%
   summarise(ncell=n(),.groups="drop")        

fig1 <- ggplot(dd)+
   geom_bar(stat="identity",
            position=position_fill(reverse=F),
            aes(x=EXP, y=ncell, fill=DROPLET.TYPE))+
   theme_bw()+
   theme(legend.title=element_blank(),
         axis.title=element_blank(), 
         axis.text.x=element_text(angle=90, hjust=1, size=8))   
                
png("./2_demux.outs/Figure1.1_droplet.png", width=600, height=500, res=120)
fig1
dev.off()


### (2)
demux3 <- demux2%>%dplyr::filter(DROPLET.TYPE=="SNG")

dd <- demux3%>%
   group_by(SNG.BEST.GUESS)%>%
   summarise(ncell=n(),.groups="drop")

fig2 <- ggplot(dd)+
   geom_bar(stat="identity", aes(x=SNG.BEST.GUESS, y=ncell), fill="#1c9099")+
   ggtitle("Number of cells")+
   theme_bw()+
   theme(legend.title=element_blank(),
         axis.title=element_blank(), 
         axis.text.x=element_text(angle=90, hjust=1, size=8),
         plot.title=element_text(hjust=0.5))

png("./2_demux.outs/Figure1.2_IND.png", width=600, height=500, res=120)
fig2
dev.off()


### (3)
dd <- demux3%>%
   group_by(EXP)%>%
   summarise(ncell=n(),.groups="drop")
fig3 <- ggplot(dd)+
   geom_bar(stat="identity", aes(x=EXP, y=ncell), fill="#1c9099")+
   ggtitle("Number of cells")+
   theme_bw()+
   theme(legend.title=element_blank(),
         axis.title=element_blank(), 
         axis.text.x=element_text(angle=90, hjust=1, size=8),
         plot.title=element_text(hjust=0.5))
png("./2_demux.outs/Figure1.3_EXP.png", width=500, height=500, res=120)
fig3
dev.off() 


### (4)
fig4.1 <- ggplot(demux2, aes(x=NUM.SNPS))+
        geom_histogram(fill="grey70", color="grey30", position="identity")+
        xlab("NUM.SNPs")+theme_bw()
 
fig4.2 <- ggplot(demux2, aes(x=NUM.READS))+
        geom_histogram(fill="grey70", color="grey30", position="identity")+
        xlab("NUM.Reads")+theme_bw() 
        
png("./2_demux.outs/Figure1.4_hist.png", width=600, height=400, res=120)
plot_grid(fig4.1, fig4.2, ncol=2)
dev.off()  







######################
### summary reads  ###
#####################
fn <- "./2_demux.outs/2_seurat.merge.SNG.rds"
atac<- read_rds(fn)

## counts <- GetAssayData(sc, slot="counts")
## nCount_ATAC <- colSums(counts)
## nFeature_ATAC <- colSums(counts>0)
###

meta <- atac@meta.data

colnames(demux)

demux_treats <- demux %>%
   dplyr::select(NEW_BARCODE, treats, EXP)

colnames(demux_treats)

head(demux_treats)

meta <- meta%>%left_join(demux_treats, by="NEW_BARCODE")

colnames(meta)

head(meta$NEW_BARCODE)

table(meta$treats)

## meta2 <- meta%>%mutate(treats=gsub(".*-ATAC-|_.*", "", NEW_BARCODE),
##                        treat2=gsub("-EtOH", "", treats),
##                        EXP=gsub("_.*","",NEW_BARCODE))

## dd2 <- meta2%>%group_by(SNG.BEST.GUESS, treat2)%>%
##    summarise(ncell=n(), reads=mean(nCount_ATAC), ngene=mean(nFeature_ATAC),.groups="drop")



## xx2 <- meta2%>%group_by(EXP)%>%
##    summarise(ncell=n(), reads=mean(nCount_ATAC), ngene=mean(nFeature_ATAC),.groups="drop")

## tmp <- dd2%>%
##            group_by(treat2)%>%
##            summarise(nind=n(), ncell=median(ncell), reads=median(reads),ngene=median(ngene),.groups="drop")


## col1 <- c("CTRL"="#828282",
##    "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
##    "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
## lab1 <- c("CTRL"="CTRL",
##    "LPS"="LPS", "LPS-DEX"="LPS+DEX",
##    "PHA"="PHA", "PHA-DEX"="PHA+DEX")


## fig1 <- ggplot(dd2, aes(x=treat2, y=ncell, fill=treat2))+
##    geom_violin()+xlab("")+ylab("")+
##    ggtitle("#Cells per individual")+
##    scale_fill_manual(values=col1)+
##    scale_x_discrete(labels=lab1)+
##    theme_bw()+
##    theme(legend.position="none",
##    plot.title=element_text(hjust=0.5),
##    axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5))

## ###
## fig2 <- ggplot(dd2,aes(x=treat2, y=reads, fill=treat2))+
##    geom_violin()+xlab("")+ylab("")+
##    ggtitle("#UMIs per cell")+
##    scale_fill_manual(values=col1)+
##    scale_x_discrete(labels=lab1)+
##    theme_bw()+
##    theme(legend.position="none",
##          plot.title=element_text(hjust=0.5),
##          axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5))

## ###
## fig3 <- ggplot(dd2,aes(x=treat2, y=ngene, fill=treat2))+
##     geom_violin()+xlab("")+ylab("")+
##     ggtitle("#Genes per cell")+
##     scale_fill_manual(values=col1)+
##     scale_x_discrete(labels=lab1)+
##     theme_bw()+
##     theme(legend.position="none",
##           plot.title=element_text(hjust=0.5),
##           axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5))

## ## png("./2_kb2_output/Figure2.5_violin.png", width=800, height=500, res=120)
## pdf("./2_demux.outs/Figure2.1_violin.pdf", width=8, height=5)
## print(plot_grid(fig1, fig2, fig3, labels="AUTO", label_fontface="plain", label_x=0.1,  ncol=3))
## dev.off()


dd2 <- meta%>%group_by(SNG.BEST.GUESS, treats)%>%
   summarise(ncell=n(), reads=mean(nCount_ATAC), ngene=mean(nFeature_ATAC),.groups="drop")



xx2 <- meta%>%group_by(EXP)%>%
   summarise(ncell=n(), reads=mean(nCount_ATAC), ngene=mean(nFeature_ATAC),.groups="drop")

tmp <- dd2%>%
           group_by(treats)%>%
           summarise(nind=n(), ncell=median(ncell), reads=median(reads),ngene=median(ngene),.groups="drop")


col1 <- c("caffeine"="#828282", "etOH"="#fb9a99",
          "nicotine"="#e31a1c", "vitA"="#ff3300",
          "vitD"="#ff6600", "vitE"="#ff9900",
          "water"="#a6cee3", "zinc"="#1f78b4")


## lab1 <- c("CTRL"="CTRL",
##    "LPS"="LPS", "LPS-DEX"="LPS+DEX",
##    "PHA"="PHA", "PHA-DEX"="PHA+DEX")


fig1 <- ggplot(dd2, aes(x=treats, y=ncell, fill=treats))+
   geom_violin()+xlab("")+ylab("")+
   ggtitle("#Cells per individual")+
   scale_fill_manual(values=col1)+
   theme_bw()+
   theme(legend.position="none",
   plot.title=element_text(hjust=0.5),
   axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5))

###
fig2 <- ggplot(dd2,aes(x=treats, y=reads, fill=treats))+
   geom_violin()+xlab("")+ylab("")+
   ggtitle("#UMIs per cell")+
   scale_fill_manual(values=col1)+
   theme_bw()+
   theme(legend.position="none",
         plot.title=element_text(hjust=0.5),
         axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5))

###
fig3 <- ggplot(dd2,aes(x=treats, y=ngene, fill=treats))+
    geom_violin()+xlab("")+ylab("")+
    ggtitle("#Genes per cell")+
    scale_fill_manual(values=col1)+
    theme_bw()+
    theme(legend.position="none",
          plot.title=element_text(hjust=0.5),
          axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5))

## png("./2_kb2_output/Figure2.5_violin.png", width=800, height=500, res=120)
pdf("./2_demux.outs/Figure2.1_violin.pdf", width=8, height=5)
print(plot_grid(fig1, fig2, fig3, labels="AUTO", label_fontface="plain", label_x=0.1,  ncol=3))
dev.off()





########################################################################
########################################################################
## RNA ##
########################################################################
########################################################################
rm(list=ls())

fn <- "./1_Merge.outs/1_seurat.merge.combined.rds"

combined <- read_rds(fn)

combined

colnames(combined@meta.data)

head(Cells(combined))

head(combined@meta.data$NEW_BARCODE)



demux <- read_rds("./2_demux.outs/1_demux.rds")

colnames(demux)




## png("./2_demux.outs/scatter_ATAC_ori.png")
## p <- plot(combined@meta.data$nFeature_ATAC, combined@meta.data$nCount_ATAC, xlab="ATAC_features", ylab="ATAC_counts")
## abline(lm(combined@meta.data$nCount_ATAC~combined@meta.data$nFeature_ATAC), col="red") # regression line (y~x)
## lines(lowess(combined@meta.data$nFeature_ATAC,combined@meta.data$nCount_ATAC), col="blue") # lowess line (x,y)
## print(p)
## dev.off()







#########################################
### subset SNG cells of seurat object ###
#########################################
combined <- RenameCells(combined, new.names=gsub("-1","",Cells(combined)))
colnames(combined@meta.data)
head(Cells(combined))


demux3 <- demux %>%
   dplyr::select(NEW_BARCODE, SNG.BEST.GUESS, DROPLET.TYPE, EXP, treats)
head(demux3)

meta <- combined@meta.data%>%mutate("NEW_BARCODE"=gsub("-1", "", barcode))
head(meta$NEW_BARCODE)
head(rownames(meta))

meta <- meta%>%left_join(demux3, by="NEW_BARCODE")
colnames(meta)

rownames(meta) <- meta$NEW_BARCODE

combined <- AddMetaData(combined, meta)
colnames(combined@meta.data)
head(rownames(combined@meta.data))
table(combined@meta.data$DROPLET.TYPE)
table(combined@meta.data$SNG.BEST.GUESS)

##
## combined2 <- subset(combined, subset=DROPLET.TYPE=="SNG")
## colnames(combined2@meta.data)
## table(combined2@meta.data$DROPLET.TYPE)
## table(combined2@meta.data$SNG.BEST.GUESS)

opfn <- "./2_demux.outs/2_seurat.merge.demux.combined.rds"
write_rds(combined, opfn)

combined <- readRDS("./2_demux.outs/2_seurat.merge.demux.combined.rds")

combined

colnames(combined@meta.data)

DefaultAssay(combined) <- "RNA"

combined[["percent.mt.RNA"]] <- PercentageFeatureSet(combined, pattern = "^MT-")




##--------------------------------------------
combined2 <- subset(combined, subset= (SNG.BEST.GUESS=="AL-002" | SNG.BEST.GUESS=="AL-018" | SNG.BEST.GUESS=="AL-026" | SNG.BEST.GUESS=="AL-030" | SNG.BEST.GUESS=="AL-041" | SNG.BEST.GUESS=="AL-044" | SNG.BEST.GUESS=="AL-045" | SNG.BEST.GUESS=="AL-056" | SNG.BEST.GUESS=="AL-058" | SNG.BEST.GUESS=="AL-062" | SNG.BEST.GUESS=="AL-070" | SNG.BEST.GUESS=="AL-105"))
table(combined2@meta.data$DROPLET.TYPE)
table(combined2@meta.data$SNG.BEST.GUESS)

combined2 <- subset(combined2, subset=DROPLET.TYPE=="SNG")
colnames(combined2@meta.data)
table(combined2@meta.data$DROPLET.TYPE)
table(combined2@meta.data$SNG.BEST.GUESS)


opfn <- "./2_demux.outs/2_seurat.merge.SNG.combined.rds"
write_rds(combined2, opfn)

combined2 <- read_rds("./2_demux.outs/2_seurat.merge.SNG.combined.rds")

colnames(combined2@meta.data)


#-----------------------------------------------
colnames(combined@meta.data)

table(combined@meta.data$DROPLET.TYPE)

#combined@meta.data$droplet.type.num < NULL

number <- c("SNG"="0", "DBL"="1", "AMB"="2")

combined@meta.data$droplet.type.num <- number[as.character(combined@meta.data$DROPLET.TYPE)]

table(combined@meta.data$droplet.type.num)

str(combined@meta.data$droplet.type.num)




#############################################################################
# filter out low quality cells
#############################################################################
## png("./2_demux.outs/scatter_RNA.png")
## p <- plot(combined@meta.data$nFeature_RNA, combined@meta.data$nCount_RNA, xlab="RNA_features", ylab="RNA_counts", col = combined@meta.data$droplet.type.num)
## abline(lm(combined@meta.data$nCount_RNA~combined@meta.data$nFeature_RNA), col="red") # regression line (y~x)
## lines(lowess(combined@meta.data$nFeature_RNA,combined@meta.data$nCount_RNA), col="blue") # lowess line (x,y)
## print(p)
## dev.off()

## png("./2_demux.outs/scatter_ATAC.png")
## p <- plot(pbmc.combined@meta.data$nFeature_ATAC, pbmc.combined@meta.data$nCount_ATAC, xlab="ATAC_features", ylab="ATAC_counts")
## abline(lm(pbmc.combined@meta.data$nCount_ATAC~pbmc.combined@meta.data$nFeature_ATAC), col="red") # regression line (y~x)
## lines(lowess(pbmc.combined@meta.data$nFeature_ATAC,pbmc.combined@meta.data$nCount_ATAC), col="blue") # lowess line (x,y)
## print(p)
## dev.off()

## png("./2_demux.outs/scatter_peaks.png")
## p <- plot(pbmc.combined@meta.data$nFeature_peaks, pbmc.combined@meta.data$nCount_peaks, xlab="peaks_features", ylab="peaks_counts")
## abline(lm(pbmc.combined@meta.data$nCount_peaks~pbmc.combined@meta.data$nFeature_peaks), col="red") # regression line (y~x)
## lines(lowess(pbmc.combined@meta.data$nFeature_peaks,pbmc.combined@meta.data$nCount_peaks), col="blue") # lowess line (x,y)
## print(p)
## dev.off()

## png("./2_demux.outs/scatter_macs2.png")
## p <- plot(pbmc.combined@meta.data$nFeature_macs2, pbmc.combined@meta.data$nCount_macs2, xlab="macs2_features", ylab="macs2_counts")
## abline(lm(pbmc.combined@meta.data$nCount_macs2~pbmc.combined@meta.data$nFeature_macs2), col="red") # regression line (y~x)
## lines(lowess(pbmc.combined@meta.data$nFeature_macs2,pbmc.combined@meta.data$nCount_macs2), col="blue") # lowess line (x,y)
## print(p)
## dev.off()

## png("./2_demux.outs/tss_enrichment.png")
## p <- plot(pbmc.combined@meta.data$TSS.percentile, pbmc.combined@meta.data$TSS.enrichment, xlab="tss_percentile", ylab="tss_enrichment")
## abline(lm(pbmc.combined@meta.data$TSS.enrichment~pbmc.combined@meta.data$TSS.percentile), col="red") # regression line (y~x)
## lines(lowess(pbmc.combined@meta.data$TSS.percentile,pbmc.combined@meta.data$TSS.enrichment), col="blue") # lowess line (x,y)
## print(p)
## dev.off()

## png("./2_demux.outs/nucleosome_signal.png")
## p <- plot(pbmc.combined@meta.data$nucleosome_percentile, pbmc.combined@meta.data$nucleosome_signal, xlab="nucleosome_percentile", ylab="nucleosome_signal")
## abline(lm(pbmc.combined@meta.data$nucleosome_signal~pbmc.combined@meta.data$nucleosome_percentile), col="red") # regression line (y~x)
## lines(lowess(pbmc.combined@meta.data$nucleosome_percentile,pbmc.combined@meta.data$nucleosome_signal), col="blue") # lowess line (x,y)
## print(p)
## dev.off()



#----------------------------------------------
colnames(combined@meta.data)

png("./2_demux.outs/scatter_RNA_ggplot_log.png")
p <- ggplot(combined@meta.data, aes(x=log10(nFeature_RNA), y=log10(nCount_RNA), shape=DROPLET.TYPE, color=DROPLET.TYPE)) +
      geom_point()
print(p)
dev.off()

png("./2_demux.outs/scatter_ATAC_ggplot_log.png")
p <- ggplot(combined@meta.data, aes(x=log10(nFeature_ATAC), y=log10(nCount_ATAC), shape=DROPLET.TYPE, color=DROPLET.TYPE)) +
      geom_point()
print(p)
dev.off()

## png("./2_demux.outs/scatter_RNA_ggplot.png")
## p <- ggplot(combined@meta.data, aes(x=nFeature_RNA, y=nCount_RNA, shape=DROPLET.TYPE, color=DROPLET.TYPE)) +
##       geom_point()
## print(p)
## dev.off()

## png("./2_demux.outs/scatter_RNA_ggplot.png")
## p <- ggplot(combined@meta.data, aes(x=nFeature_RNA, y=nCount_RNA, shape=DROPLET.TYPE, color=DROPLET.TYPE)) +
##       geom_point()
## print(p)
## dev.off()

png("./2_demux.outs/scatter_nucleosome_ggplot.png")
p <- ggplot(combined@meta.data, aes(x=nucleosome_percentile, y=nucleosome_signal, shape=DROPLET.TYPE, color=DROPLET.TYPE)) +
      geom_point()
print(p)
dev.off()

png("./2_demux.outs/scatter_TSS_ggplot.png")
p <- ggplot(combined@meta.data, aes(x=TSS.percentile, y=TSS.enrichment, shape=DROPLET.TYPE, color=DROPLET.TYPE)) +
      geom_point()
print(p)
dev.off()


#--------
png("./2_demux.outs/nucleosome_ggplot2_violin.png")
p <- ggplot(combined@meta.data, aes(x=DROPLET.TYPE, y=nucleosome_signal)) +
      geom_violin()
print(p)
dev.off()


png("./2_demux.outs/TSS_ggplot2_violin.png")
p <- ggplot(combined@meta.data, aes(x=DROPLET.TYPE, y=TSS.enrichment)) +
      geom_violin()
print(p)
dev.off()

png("./2_demux.outs/mt_RNA_ggplot2_violin.png")
p <- ggplot(combined@meta.data, aes(x=DROPLET.TYPE, y=percent.mt.RNA)) +
      geom_violin()
print(p)
dev.off()

max(combined$percent.mt.RNA)
     
