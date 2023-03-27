##
library(Matrix)
library(tidyverse)
## library(parallel)
## library(data.table)
## library(future)
## library(purrr)
## library(furrr)
## library("BiocParallel")
## library(Rcpp)
## library(reshape)
library(qqman)
library(qvalue)
##
library(DESeq2)
library(biobroom)
library(ashr)
library(GenomicRanges)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(Signac)
library(clusterProfiler)
library(org.Hs.eg.db)
library(annotables)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(EnsDb.Hsapiens.v75)
library(ChIPseeker,lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
###
library(ggplot2)
library(cowplot,lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(grid)
library(gridExtra)
library(ggExtra)
library(RColorBrewer)
library(ggsci)
library(RColorBrewer)
library(viridis)
theme_set(theme_grey())

rm(list=ls())

outdir <- "./3.2_Example.outs/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


###
### read atac data
atac <- read_rds("../3_motif/1_motif.outs/1_scATAC.motif.rds")
## x <- atac@meta.data
## x$treat <- gsub(".*-ATAC-|_.*", "", x$NEW_BARCODE)
## atac <- AddMetaData(atac,x)

###
### peak annotation
fn <- "../2_Differential/2.2_compareRNAandATAC.outs/2_annot.ChIPseeker.rds"
anno <- read_rds(fn)%>%as.data.frame()%>%
       mutate(peak=paste(gsub("chr","",seqnames), start, end, sep="-"))



### Fun, getPeaks
getPeaks <- function(atac, resMotif, resDP){
   ###
   ### motif data   
   motif <- Motifs(atac)
   x <- motif@data
   peaks <- rownames(x)
    
   motifDf <- resMotif%>%distinct(gene, .keep_all=T)%>%
       dplyr::select(ID=gene, motif.name=motif)%>%as.data.frame()
    
   ## motif.ID <- ConvertMotifID(motif, motif.name)
   peaks_motif <- map_dfr(1:nrow(motifDf), function(i){ 
      ###
      ID <- motifDf[i,1]
      motif.name <- motifDf[i,2] 
       
   ### peak containing motifs
      peak.sel <- peaks[x[,ID]==1] ## peak containing motifs
      resDP2 <- resDP%>%dplyr::filter(gene%in%peak.sel)
      ##DEX
      tmp <- resMotif%>%dplyr::filter(gene==ID, grepl("DEX", contrast), qval<0.5, abs(beta)>1.41)
      peak.DEX <- NULL
      if ( nrow(tmp)>0){ 
         Bmotif <- max(tmp$beta)
         peak.DEX <- resDP%>%
            dplyr::filter(gene%in%peak.sel,
            sign(estimate)==sign(Bmotif), grepl("DEX", contrast))%>%dplyr::pull(gene)%>%unique()
      }
      ### immune stimuli
      tmp <- resMotif%>%dplyr::filter(gene==ID, !grepl("DEX", contrast), qval<0.5, abs(beta)>1.41)
      peak.stimuli <- NULL 
      if ( nrow(tmp)>0){ 
         Bmotif <- max(tmp$beta)
         peak.stimuli <- resDP%>%
             dplyr::filter(gene%in%peak.sel,
             sign(estimate)==sign(Bmotif), !grepl("DEX", contrast))%>%dplyr::pull(gene)%>%unique()
      }
   ##    
       peak.final <- union(peak.DEX, peak.stimuli)
       ##
       peakDf <- data.frame(peaks=peak.final, motif.ID=ID, motif.name=as.character(motif.name))
       peakDf
    })
    ###
    peaks_motif
}

    

    

####
#### fun-1, get gene id
getGene <- function(anno, peakSel){
   ##
   anno2 <- anno%>%dplyr::filter(peak%in%peakSel, abs(distanceToTSS)<100000) ##100000
   ### gene list in a set of peaks that contain specific motif and are differentially accessible
   gene <- unique(anno2$geneId)
   gene 
}    


###
### fun-2, average gene expression
avePathway <- function(X){

   ### filtering more missing value
   ii <- apply(!is.na(X), 1, sum)
   X <- X[ii>20,]
   bti <- colnames(X)
   cvt <- str_split(bti, "_", simplify=T)
   cvt <- data.frame(rn=bti, MCls=cvt[,1],treats=cvt[,2],sampleID=cvt[,3],Batch=cvt[,4])%>%
       mutate(comb=paste(MCls, Batch, sep="_"))
   comb <- unique(cvt$comb)
   for (ii in comb){
      bti0 <- cvt%>%dplyr::filter(comb==ii)%>%dplyr::pull(rn)
      x <- X[,bti0]
      x.mean <- apply(x, 1, mean, na.rm=T)
      x.scale <- sweep(x, 1, x.mean, "-")
      X[,bti0] <- x.scale
   }
   pathway <- apply(X, 2, mean, na.rm=T)
}

###
### fun-3
getData <- function(gene, datatype="bulk"){

   ### bulk 
   if ( datatype=="bulk"){
      fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew_output/Filter2/YtX_sel.comb.RData"
      load(fn)
      rn <- gsub("\\..*", "", rownames(YtX_sel))
      rownames(YtX_sel) <- rn
      X <- YtX_sel[rn%in%gene,]+1
      counts <- colSums(YtX_sel)
      X <- sweep(X, 2, counts, "/")
      X <- X*1e+06
      X <- log2(X)
     
      bti <- colnames(X)
      cvt <- str_split(bti, "_", simplify=T)
      cvt <- data.frame(rn=bti, MCls=cvt[,1], treats=gsub("-EtOH", "", cvt[,2]),
                        sampleID=cvt[,3], Batch=cvt[,4])
      cvt$y <- avePathway(X)
   ##
   } ### end bulk
    
    
   #### variability    
   if ( datatype=="NB.phi"){
      fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/10_RNA.Variance_output/tmp10/1.2_Sel.PhxNew.RData"
      load(fn)
      rn <- gsub("\\..*", "", rownames(PhxNew2))
      rownames(PhxNew2) <- rn

      X <- PhxNew2[rn%in%gene,]
      X <- log2(X)
      bti <- colnames(PhxNew2)
      cvt <- str_split(bti, "_", simplify=T)
      cvt <- data.frame(rn=bti, MCls=cvt[,1], treats=gsub("-EtOH", "", cvt[,2]), sampleID=cvt[,3], Batch=cvt[,4])
      cvt$y <- avePathway(X)
   } ### end variability

    
   MCls <- sort(unique(cvt$MCls))
   ### 
   cvt <- map_dfr(MCls, function(oneMCl){
      cvt2 <- cvt%>%dplyr::filter(MCls==oneMCl)
      df0 <- cvt2%>%dplyr::filter(treats=="CTRL")%>%
         dplyr::select(sampleID, y)%>%
         dplyr::rename("y0"="y")
      cvt2 <- cvt2%>%left_join(df0, by=c("sampleID"))%>%mutate(y2=y-y0)
      cvt2
   })
   ###    
   cvt
}
    


## lab1 <- c("CTRL"="CTRL",
lab1 <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX",
   "PHA"="PHA", "PHA-DEX"="PHA+DEX")
## col1 <- c("CTRL"="#828282",
col1 <- c("LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
   "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")


## motif.name <- "RELA"
res <- read_rds("../2_Differential/1.2_DiffPeak.outs/2.0_DESeq.results.rds")%>%
    as.data.frame()

motifList <- c("NR3C1", "NR3C2", "RELA", "FOS::JUN", "NFKB1", "NFKB2")


##################
### expression ###
##################

## NR3C2
###
###
for (i in 1:2){
   ##
   motif <- motifList[i]
   cat(i, motif, "\n") 
 ### differential peaks 
   res2 <- res%>%
      dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5, estimate>0, MCls=="Tcell", grepl("DEX",contrast))
    
   gene <- getGene(atac, anno, motif.name=motif, resDP=res2)
   cvt <- getData(gene=gene, datatype="bulk")
###
###
   cvt2 <- cvt%>%dplyr::filter(treats!="CTRL", MCls=="Tcell", Batch=="SCAIP6")

   p <- ggplot(cvt2,aes(x=factor(treats), y=y, fill=treats))+
      geom_boxplot(outlier.size=0.8)+
      ylab("Expression value")+
      scale_y_continuous(expand=expansion(mult=c(0, 0.2)))+
      scale_fill_manual("", values=col1, labels=lab1)+
      scale_x_discrete("", labels=lab1)+
      ggtitle(bquote(~italic(.(motif))~"(T cell)"))+
      theme_bw()+
      theme(axis.text.x=element_text(angle=-90, hjust=0, size=10),
          axis.title.x=element_blank(),
          axis.title.y=element_text(size=12),
          plot.title=element_text(hjust=0.5, size=12),
          legend.position="none")

###
   figfn <- paste("./3.2_Example.outs/Figure", i, "_", motif, ".boxplot.png", sep="")
   png(figfn, width=480, height=420, res=120)
   print(p)
   dev.off()
###
}##

###
###
for (i in 3:5){
   motif <- motifList[i]
   cat(i, motif, "\n")
 ### differential peaks 
   res2 <- res%>%
      dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5, estimate<0, MCls=="Tcell", grepl("DEX",contrast))
    
   gene <- getGene(atac, anno, motif.name=motif, resDP=res2)
   cvt <- getData(gene=gene, datatype="bulk")
###
###
   cvt2 <- cvt%>%dplyr::filter(treats!="CTRL", MCls=="Tcell", Batch=="SCAIP6")

   p <- ggplot(cvt2,aes(x=factor(treats), y=y, fill=treats))+
      geom_boxplot(outlier.size=0.8)+
      ylab("Expression value")+
      scale_y_continuous(expand=expansion(mult=c(0, 0.2)))+
      scale_fill_manual("", values=col1, labels=lab1)+
      scale_x_discrete("", labels=lab1)+
      ggtitle(bquote(~italic(.(motif))~"(T cell)"))+
      theme_bw()+
      theme(axis.text.x=element_text(angle=-90, hjust=0, size=10),
          axis.title.x=element_blank(),
          axis.title.y=element_text(size=12),
          plot.title=element_text(hjust=0.5, size=12),
          legend.position="none")

###
   figfn <- paste("./3.2_Example.outs/Figure", i, "_", motif, ".boxplot.png", sep="")
   png(figfn, width=480, height=420, res=120)
   print(p)
   dev.off()
}



###################
### variability ###
###################

for (i in 1:2){
   ##
   motif <- motifList[i]
   cat(i, motif, "\n") 
 ### differential peaks 
   res2 <- res%>%
      dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5, estimate>0, MCls=="Tcell", grepl("DEX",contrast))
    
   gene <- getGene(atac, anno, motif.name=motif, resDP=res2)
   cvt <- getData(gene=gene, datatype="NB.phi")
###
###
   cvt2 <- cvt%>%dplyr::filter(treats!="CTRL", MCls=="Tcell", Batch=="SCAIP6")

   p <- ggplot(cvt2,aes(x=factor(treats), y=y, fill=treats))+
      geom_boxplot(outlier.size=0.8)+
      ylab("Variability")+
      scale_y_continuous(expand=expansion(mult=c(0, 0.2)))+
      scale_fill_manual("", values=col1, labels=lab1)+
      scale_x_discrete("", labels=lab1)+
      ggtitle(bquote(~italic(.(motif))~"(T cell)"))+
      theme_bw()+
      theme(axis.text.x=element_text(angle=-90, hjust=0, size=10),
          axis.title.x=element_blank(),
          axis.title.y=element_text(size=12),
          plot.title=element_text(hjust=0.5, size=12),
          legend.position="none")

###
   figfn <- paste("./3.2_Example.outs/Figure", i, "_", motif, "_va.boxplot.png", sep="")
   png(figfn, width=480, height=420, res=120)
   print(p)
   dev.off()
###
}##



###########################################################
### compare motif activties and transcriptional changes ###
###########################################################

fn <- "../2_Differential/1.2_DiffPeak.outs/2.0_DESeq.results.rds"
resDP <- read_rds(fn)%>%drop_na(p.value)%>%
    mutate(comb=paste(MCls, contrast, sep="_"))%>%
    dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)%>%as.data.frame()


fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew_output/Filter2/1_DESeq.results.rds"
resDG <- read_rds(fn)%>%drop_na(p.value)%>%dplyr::filter(Batch2=="SCAIP6")%>%
    mutate(comb=paste(MCls, contrast, sep="_"))%>%as.data.frame()


resMotif <- read_rds("./2_motif.activities.outs/3_motif.diff.results.rds")%>%
    mutate(comb=paste(MCls, contrast, sep="_"))%>%as.data.frame()

topmotif <- resMotif%>%dplyr::filter(qval<0.1, abs(beta)>1.41)%>%dplyr::pull(gene)%>%unique()

comb <- sort(unique(resMotif$comb))

resMotif2 <- resMotif%>%dplyr::filter(gene%in%topmotif)

### motif >> peaks



###
###

peak_motif <- getPeaks(atac, resMotif2, resDP)

comb <- sort(unique(resMotif2$comb))

plotDf <- map_dfr(comb, function(ii){
   ###
   res2 <- resMotif2%>%dplyr::filter(comb==ii)
   LFC.RNA <- sapply(res2$motif, function(mm){
       ##
       peakSel <- peak_motif%>%dplyr::filter(motif.name==mm)%>%pull(peaks)
       gene2 <- getGene(anno, peakSel)
       LFC <- resDG%>%dplyr::filter(comb==ii, gene%in%gene2)%>%dplyr::pull(estimate)
       median(LFC)
   })
   res2$LFC.RNA <- LFC.RNA
   res2
})



## df2 <- plotdf%>%group_by(MCls)%>%mutate(min.x=min(beta.x), max.x=max(beta.x))%>%ungroup()%>%
##    group_by(contrast)%>%mutate(min.y=min(beta.y), max.y=max(beta.y))%>%ungroup() 
## anno_df2 <- df2%>%
##     group_by(contrast, MCls)%>%
##     nest()%>%
##     mutate(corr=map(data, ~cor.test((.x)$beta.x, (.x)$beta.y, method="pearson")),
##           eq=map(corr,feq),
##           r2=map_dbl(corr,~(.x)$estimate),
##           xpos=map_dbl(data,~xFun(.x, a=0.7)),
##           ypos=map_dbl(data,~yFun(.x, a=1)))%>%
##    dplyr::select(-data,-corr)
library(ggrastr)
##
plotDf <- plotDf%>%dplyr::rename("beta.motif"="beta", "beta.RNA"="LFC.RNA")
plotDf2 <- plotDf%>%dplyr::filter(MCls=="Tcell")
cor.test(plotDf2$beta.motif, plotDf2$beta.RNA)

plotDf2 <- plotDf2%>%mutate(gr2=ifelse(motif%in%c("NR3C1", "NFKB1"), contrast, 0)) 
     
annoDf <- plotDf2%>%dplyr::filter(motif%in%c("NR3C1", "NFKB1"))%>%
    mutate(motif2=paste("italic(", motif, ")", sep="")) 

## eq <- bquote(italic(R)==~"0.814,"~"***")
  
p <- ggplot(plotDf2)+
    rasterise(geom_point(aes(x=beta.motif, y=beta.RNA, colour=gr2), size=0.8),
                   dpi=300)+
    annotate("text", x=-1.5, y=0.23, label="italic(R)==0.630", parse=T, color="blue")+
    ggrepel::geom_text_repel(data=annoDf,
        aes(x=beta.motif, y=beta.RNA, label=motif2, colour=gr2), size=3, parse=T,
        box.padding=0.8, max.overlaps=Inf)+
    scale_colour_manual(values=c("0"="grey50",
                                 "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
                                 "PHA"="#a6cee3", "PHA-DEX"="#1f78b4"))+ 
    ## facet_grid(contrast~MCls, scales="fixed",
    ##    labeller=labeller(
    ##       contrast=as_labeller(c("LPS"="LPS", "LPS-DEX"="LPS+DEX",
    ##                              "PHA"="PHA", "PHA-DEX"="PHA+DEX")) ))+
    ggtitle("T-cell")+
    scale_x_continuous("LFC on motif activities", expand=expansion(mult=0.1))+
    scale_y_continuous("LFC on motif regulated genes transcription", expand=expansion(mult=0.1))+
    theme_bw()+
    theme(legend.position="none",
          axis.title=element_text(size=12),
          plot.title=element_text(hjust=0.5))

p2 <- p+geom_smooth(data=plotDf2, aes(x=beta.motif, y=beta.RNA),
   method="lm", formula=y~x, color="blue", size=0.5, se=F)

## slides
figfn <- paste(outdir, "Figure0.1_motif.png", sep="")
png(filename=figfn, width=380, height=380, res=100)  
print(p2)
dev.off()


### poster                  
figfn <- paste(outdir, "Figure0.1.1_motif.png", sep="")
png(filename=figfn, width=300, height=450, res=100)  
print(p2)
dev.off()


### nearby gene expression
## gene0 <- res5$geneId[1]
## symbol0 <- res5$SYMBOL[1]
## cvt <- getData(gene=gene0, datatype="RNA")
## ###
## ###
## cvt2 <- cvt%>%dplyr::filter(treats!="CTRL", MCls=="Tcell")
## p <- ggplot(cvt2,aes(x=factor(treats), y=yscale2, fill=treats))+
##    geom_boxplot(outlier.size=0.8)+
##    ylab(bquote(~log[2]~"(Expression)"))+
##    scale_y_continuous(expand=expansion(mult=c(0, 0.2)))+
##    scale_fill_manual("", values=col1, labels=lab1)+
##    scale_x_discrete("", labels=lab1)+
##    ggtitle(bquote(~italic(.(symbol0))~"(T cell)"))+
##    theme_bw()+
##    theme(axis.text.x=element_text(angle=-90, hjust=0, size=10),
##          axis.title.x=element_blank(),
##          axis.title.y=element_text(size=12),
##          plot.title=element_text(hjust=0.5, size=12),
##          legend.position="none")

## ###
## figfn <- paste("./5_Example.outs/", motif.name, "/Figure", i, ".3_", symbol0, ".boxplot.png", sep="")
## png(figfn, width=480, height=420, res=120)
## print(p)
## dev.off()

## ### gene variability
## cvt <- getData(gene=gene0, datatype="NB.phi")
## cvt2 <- cvt%>%dplyr::filter(treats!="CTRL", MCls=="Tcell")
## p <- ggplot(cvt2,aes(x=factor(treats), y=yscale2, fill=treats))+
##    geom_boxplot(outlier.size=0.8)+
##    ylab(bquote(~log[2]~"(Variability)"))+
##    scale_y_continuous(expand=expansion(mult=c(0, 0.2)))+
##    scale_fill_manual("", values=col1, labels=lab1)+
##    scale_x_discrete("", labels=lab1)+
##    ggtitle(bquote(~italic(.(symbol0))~"(T cell)"))+
##    theme_bw()+
##    theme(axis.text.x=element_text(angle=-90, hjust=0, size=10),
##          axis.title.x=element_blank(),
##          axis.title.y=element_text(size=12),
##          plot.title=element_text(hjust=0.5, size=12),
##          legend.position="none")

## ###
## figfn <- paste("./5_Example.outs/", motif.name, "/Figure", i, ".4_", symbol0, ".va.boxplot.png", sep="")
## png(figfn, width=480, height=420, res=120)
## print(p)
## dev.off()



## geneID <- bitr(geneID=unique(resDE$gene), fromType="ENSEMBL", toType="SYMBOL",
##    OrgDb=org.Hs.eg.db)
## geneID <- geneID%>%dplyr::rename("gene"="ENSEMBL")

## resDE <- resDE%>%dplyr::left_join(geneID, by="gene")
