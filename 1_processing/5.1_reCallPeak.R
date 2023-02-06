##
library(tidyverse)
## library(parallel)
library(data.table)
## library(purrr)
library(GenomicRanges)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(SeuratObject)
library(Signac)
library(SeuratWrappers)
library(cicero, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(monocle3)
library(EnsDb.Hsapiens.v75)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
###
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(viridis)
theme_set(theme_grey())

###
###
outdir <- "./5.1_reCallPeak.outs/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)
rm(list=ls())


#####################
### re-call peaks ###
#####################

###
# atac <- read_rds("./4.2_Integrate.outs/3_scATAC.annot.rds")
###

pbmc.combined <- read_rds("./3_Clustering.outs/1_seurat.cluster.combined.mitofilt.rds")
pbmc.combined
colnames(pbmc.combined@meta.data)


### define cell type based on cluster
#Julong
## cell <- c("0"="Tcell", "1"="Tcell", "2"="NKcell", "3"="Bcell",
##           "4"="Tcell", "5"="Tcell", "6"="Monocyte", "7"="Tcell", "8"="Tcell",
##           "9"="Tcell", "10"="Tcell", "11"="DC", "12"="DC")

#0.15
## cell <- c("0"="CD4Naive", "1"="TCM", "2"="NKcell", "3"="Bcell",
##           "4"="CD8Naive", "5"="TEM", "6"="Monocyte", "7"="CTL", "8"="dnT",
##           "9"="Treg", "10"="Platelet", "11"="DC")

#0.1
cell <- c("0"="CD4Naive", "1"="TCM", "2"="NKcell", "3"="TEM",
           "4"="Bcell", "5"="CD8Naive", "6"="Monocyte", "7"="dnT", "8"="MAIT",
           "9"="Platelet", "10"="DC")

#0.05
#cell <- c("0"="Tcell", "1"="Tcell", "2"="Tcell", "3"="NKcell",
#          "4"="Bcell", "5"="Monocyte", "6"="Platelet", "7"="DC")

pbmc.combined$MCls <- cell[as.character(pbmc.combined$wsnn_res.0.1)]

#--atac--
## atac <- pbmc.combined[["ATAC"]]
## atac
## colnames(atac@meta.data)

## atac <- CreateSeuratObject(atac, assay = "ATAC")
## atac
## colnames(atac@meta.data)
#--atac--

DefaultAssay(pbmc.combined) <- "ATAC"

peaks <- CallPeaks(pbmc.combined, group.by="MCls",
   macs2.path="/wsu/home/groups/piquelab/apps/el7/anaconda3python/envs/macs2/bin/macs2")
###
head(peaks)
length(peaks)

opfn <- "./5.1_reCallPeak.outs/1.1_CallPeak.macs.mitofilt.rds"
write_rds(peaks, opfn)

peaks <- read_rds("./5.1_reCallPeak.outs/1.1_CallPeak.macs.mitofilt.rds")
head(peaks)
length(peaks)

peaks.2 <- keepStandardChromosomes(peaks, pruning.mode = "coarse")

length(peaks.2)



#--------------trial------------
## metafn <- paste0("/wsu/home/groups/piquelab/scGxE/counts_cellranger_arc/scGxE_1-1/outs/filtered_feature_bc_matrix/barcodes.tsv.gz")
## meta <- fread(metafn, header=F)
## head(meta)


#--------------function---------
##
## frag <- Fragments(atac)
## counts <- FeatureMatrix(fragments=frag, features=peaks)
## opfn <- "./5.1_reCallPeak.outs/1.2_counts.rds"
## write_rds(counts, opfn)

####
#### readATAC funtion for read atac data 
## readATAC <- function(run, peaks){
## ### load metadata
## ### extract filtered barcodes from cell ranger
##    metafn <- paste(run, "outs/filtered_feature_bc_matrix/barcodes.tsv.gz", sep="")
##    meta <- fread(metafn, header=F)
                                                                         
##    fragfn <- paste(run, "outs/atac_fragments.tsv.gz", sep="")
## ### create fragment object
##    frags <- CreateFragmentObject(path=fragfn, cells=meta$V1)

## ### quantify peaks
##    counts <- FeatureMatrix(fragments=frags, features=peaks, cells=meta$V1)

## ### create a seurat object
##    assay <- CreateChromatinAssay(counts,
##       sep=c(":", "-"), genome="hg38", fragments=frags)
##    atac <- CreateSeuratObject(assay, assay="ATAC")
##    atac 
## }###


readATAC <- function(run, peaks){
   metafn <- paste(run, "outs/filtered_feature_bc_matrix/barcodes.tsv.gz", sep="")
   meta <- fread(metafn, header=F)                                                                         
   fragfn <- paste(run, "outs/atac_fragments.tsv.gz", sep="")
   frags <- CreateFragmentObject(path=fragfn, cells=meta$V1)

   counts <- FeatureMatrix(fragments=frags, features=peaks, cells=meta$V1)

   assay <- CreateChromatinAssay(counts,
      sep=c(":", "-"), fragments=frags)

   atac <- CreateSeuratObject(assay, assay="ATAC")
   atac 
}



##############################
### merge multiple objects ###
##############################

### folders
basefolder <- "/wsu/home/groups/piquelab/scGxE/counts_cellranger_arc/"
expNames <- dir(basefolder,"^scGxE_1-*")
folders <- paste0(basefolder, expNames, "/", sep="")
ind <- dir.exists(folders)
folders <- folders[ind]
expNames <- expNames[ind]
names(folders) <- expNames

folders
expNames

#peaks <- read_rds("./5.1_reCallPeak.outs/1.1_CallPeaks.macs.mitofilt.rds")
expNames <- names(folders)
expNames

atac_ls <- lapply(expNames, function(ii){
   cat(ii,"\n")
   atac <- readATAC(run=folders[ii], peaks=peaks.2)
   atac <- RenameCells(atac, add.cell.id=ii)
   atac <- RenameCells(atac, new.names=gsub("-1","",Cells(atac)))
   atac
})
###

#this is macs2
combined <- merge(atac_ls[[1]], atac_ls[-1], project="sc-atac")
combined

write_rds(combined, file="./5.1_reCallPeak.outs/2_scATAC.merge.mitofilt.ndup.rds")

#combined <- read_rds("./5.1_reCallPeak.outs/2_scATAC.merge.mitofilt.ndup.rds")
#combined




### Add meta data
## combined <- read_rds("./5.1_reCallPeak.outs/2_seurat.merge.rds")
###add the gene information to the object
annotations <- GetGRangesFromEnsDb(ensdb=EnsDb.Hsapiens.v86)
#seqlevelsStyle(annotations) <- "1000"
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"
Annotation(combined) <- annotations

#this is RNA+ATAC
pbmc.combined <- read_rds("./3_Clustering.outs/1_seurat.cluster.combined.mitofilt.rds")

pbmc.combined
pbmc.combined[["ATAC"]]

x <- pbmc.combined@meta.data
colnames(x)

combined2 <- subset(combined,cells=Cells(pbmc.combined))
combined2

meta <- combined2@meta.data
colnames(meta)
head(meta)

colnames(pbmc.combined@meta.data)

#metaNew <- cbind(meta,x[,28:37])

metaNew2 <- cbind(meta,x[,25:61])
colnames(metaNew2)

#--skipping this for now--
#metaNew2$SCT.weight < NULL
#metaNew2$ATAC.weight < NULL
metaNew2 <- subset(metaNew2, select = -SCT.weight)
metaNew2 <- subset(metaNew2, select = -ATAC.weight)
#--skipping this for now--


combined <- AddMetaData(combined2,metaNew2)
combined
colnames(combined@meta.data)
head(combined)
head(combined@assays$ATAC@counts)

opfn <- "./5.1_reCallPeak.outs/3_scATAC.annot.mitofilt.ndup.rds"
write_rds(combined, opfn)



###
### bed files
#this is macs2
#atac <- read_rds("./5.1_reCallPeak.outs/3_scATAC.annot.mitofilt.rds")
atac <- combined
atac

peak <- granges(atac)
chr <- seqnames(peak)
range <- ranges(peak)
bed <- data.frame(chr=chr, start=start(range), end=end(range))

write.table(bed, file="./5.1_reCallPeak.outs/peak.mitofilt.ndup.bed", quote=F, sep="\t", row.names=F, col.names=F)




###############
### summary ###
###############
#atac <- read_rds("./5.1_reCallPeak.outs/3_scATAC.annot.mitofilt.rds")
###

meta <- atac@meta.data

colnames(meta)

head(meta)

table(meta$SNG.BEST.GUESS)

table(meta$treats)

meta2 <- meta

## meta2 <- meta%>%mutate(treats=gsub(".*-ATAC-|_.*", "", NEW_BARCODE),
##                        treat2=gsub("-EtOH", "", treats),
##                        EXP=gsub("_.*","",NEW_BARCODE))


dd2 <- meta2%>%group_by(SNG.BEST.GUESS, treats)%>%
   summarise(ncell=n(), reads=mean(nCount_ATAC), ngene=mean(nFeature_ATAC),.groups="drop")


 
xx2 <- meta2%>%group_by(EXP)%>%
   summarise(ncell=n(), reads=mean(nCount_ATAC), ngene=mean(nFeature_ATAC),.groups="drop")

tmp <- dd2%>%
           group_by(treats)%>%
           summarise(nind=n(), ncell=median(ncell), reads=median(reads),ngene=median(ngene),.groups="drop")


## col1 <- c("CTRL"="#828282",
##    "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
##    "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
## lab1 <- c("CTRL"="CTRL",
##    "LPS"="LPS", "LPS-DEX"="LPS+DEX",
##    "PHA"="PHA", "PHA-DEX"="PHA+DEX")

colnames(atac@meta.data)

## col1 <- c("caffeine"="#828282", "etOH"="#fb9a99",
##                     "nicotine"="#e31a1c", "vitA"="#ff3300",
##                     "vitD"="#ff6600", "vitE"="#ff9900",
##                     "water"="#a6cee3", "zinc"="#1f78b4")


col1 <- c("caffeine"="red", "etOH"="grey",
                    "nicotine"="tan", "vitA"="tan4",
                    "vitD"="seagreen4", "vitE"="salmon3",
                    "water"="grey", "zinc"="maroon3")


fig1 <- ggplot(dd2, aes(x=treats, y=ncell, fill=treats))+
   geom_violin()+xlab("")+ylab("")+
   ggtitle("#Cells per individual")+
   scale_fill_manual(values=col1)+
   ylim(0,3000)+
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
    ggtitle("#Features per cell")+
    scale_fill_manual(values=col1)+
    theme_bw()+
    theme(legend.position="none",
          plot.title=element_text(hjust=0.5),
          axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5))

## png("./2_kb2_output/Figure2.5_violin.png", width=800, height=500, res=120)
pdf("./5.1_reCallPeak.outs/Figure2.1_violin.mitofilt.pdf", width=8, height=5)
print(plot_grid(fig1, fig2, fig3, labels="AUTO", label_fontface="plain", label_x=0.1,  ncol=3))
dev.off()

## combined <- combined%>%NucleosomeSignal()%>%TSSEnrichment(fast=F)
## ## combined <- NucleosomeSignal(combined)
## #combined <- TSSEnrichment(combined, fast=F)
## combined$pct_reads_in_peaks <- combined$peak_region_fragments/combined$passed_filters*100
## combined$blacklist_ratio <- combined$blacklist_region_fragments/combined$peak_region_fragments

## opfn <- "./1_Merge.outs/1_seurat.merge.rds"
## write_rds(combined, file=opfn)  



####
###
## atac <- read_rds("./5.1_reCallPeak.outs/3_scATAC.annot.rds")

## ###

## MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
## ##
## atac2 <- subset(atac, subset=MCls!="DC")

## FindAllMarkers(
