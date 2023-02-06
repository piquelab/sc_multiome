#module load R/test_4.1
###
#setwd("/nfs/rprdata/julong/sc-atac/analyses.2021-02-05/1_seurat")
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
library(EnsDb.Hsapiens.v75)
library(EnsDb.Hsapiens.v86)
## library(cicero, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
## library(monocle3)
###
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(viridis)
theme_set(theme_grey())

outdir <- "./1_Merge.outs/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)

rm(list=ls())

###########################
### 0, generate folders ###
###########################

basefolder <- "/wsu/home/groups/piquelab/scGxE/counts_cellranger_arc/"

expNames <- dir(basefolder,"^scGxE_1-[1-8]")
expNames

folders <- paste0(basefolder, expNames, "/", sep="")
folders

ind <- dir.exists(folders)
ind

folders <- folders[ind]
folders

expNames <- expNames[ind]
expNames

names(folders) <- expNames
folders







##################################################################
### 1, creating a common peak set from filtered peak.bed files ###
##################################################################

expNames <- names(folders)
expNames

peaks <- lapply(expNames, function(ii){
   fn <- paste(folders[ii], "outs/atac_peaks.bed", sep="")
   bed <- read.table(fn, col.names=c("chr", "start", "end"))  
   gr <- makeGRangesFromDataFrame(bed)
   gr  
})

head(peaks)

x2 <- unlist(GRangesList(peaks))
x2

combined.peaks <- GenomicRanges::reduce(x=x2)
###

peakwidths <- GenomicRanges::width(combined.peaks)

combined.peaks <- combined.peaks[peakwidths<10000&peakwidths>20]

write_rds(combined.peaks, file="./1_Merge.outs/0_combined.peaks.rds")

combined.peaks





#####################################
### 3, combined seurat object
#####################################

### readATAC, quantify peaks and create ATAC object
readATAC <- function(run, peaks){
    metafn <- paste(run, "outs/filtered_feature_bc_matrix/barcodes.tsv.gz", sep="")
    meta <- fread(metafn, header=F)
    fragfn <- paste(run, "outs/atac_fragments.tsv.gz", sep="")
    frags <- CreateFragmentObject(path=fragfn, cells=meta$V1)
    counts <- FeatureMatrix(fragments=frags, features=peaks, cells=meta$V1)
    assay <- CreateChromatinAssay(counts,
                                  sep=c(":", "-"), genome="hg38", fragments=frags)
    atac <- CreateSeuratObject(assay, assay="ATAC")
    atac
}






####################################
### merge multiple objects
####################################
combined.peaks <- read_rds("./1_Merge.outs/0_combined.peaks.rds")

expNames <- names(folders)
expNames

atac_ls <- lapply(expNames, function(ii){
   cat(ii,"\n")
   atac <- readATAC(run=folders[ii], peaks=combined.peaks)
   atac <- RenameCells(atac, add.cell.id=ii)
   atac <- RenameCells(atac, new.names=gsub("-1","",Cells(atac)))
   atac
})
###

atac_ls

combined <- merge(atac_ls[[1]], atac_ls[-1], project="sc-atac")
combined[["ATAC"]]
write_rds(combined, file="./1_Merge.outs/1_seurat.merge.rds")

combined <- read_rds("./1_Merge.outs/1_seurat.merge.rds")

head(combined@meta.data)

combined@meta.data[rownames(combined@meta.data) %like% "AAACAGCCAAACAACA", ]





####################################
### Add meta data
####################################
combined <- read_rds("./1_Merge.outs/1_seurat.merge.rds")


### read meta data across 10 experiments
meta <- map_dfr(expNames, function(ii){
   run <- folders[ii]
   fn <- paste(run, "outs/per_barcode_metrics.csv", sep="")
   meta <- fread(fn)%>%mutate(barcode=paste(ii, barcode, sep="_"))
   meta <- meta%>%mutate(barcode=gsub("-1", "", barcode))
})

colnames(meta)
head(rownames(meta))
#head(meta)
#head(meta$barcode)

colnames(combined@meta.data)
head(rownames(combined@meta.data))
#head(combined@meta.data)

x <- combined@meta.data
x <- x%>%mutate(barcode=rownames(x))%>%left_join(meta)

colnames(x)
head(rownames(x))
#head(x)
#head(x$barcode)

rownames(x) <- x$barcode
colnames(x)
head(rownames(x))

combined <- AddMetaData(combined, x)
colnames(combined@meta.data)
head(rownames(combined@meta.data))
#head(combined@meta.data)

combined@meta.data[rownames(combined@meta.data) %like% "scGxE_1-4_AAACAGCCAAACAACA", ]





###add the gene information to the object
annotations <- GetGRangesFromEnsDb(ensdb=EnsDb.Hsapiens.v86)
#seqlevelsStyle(annotations) <- "1000"
genome(annotations) <- "hg38"
Annotation(combined) <- annotations



combined <- combined%>%NucleosomeSignal()%>%TSSEnrichment(fast=F)
## combined <- NucleosomeSignal(combined)
#combined <- TSSEnrichment(combined, fast=F)





## single.csv <- fread("/nfs/rprdata/julong/sc-atac/count.SCAIP.2021-01-14/SCAIP6A-ATAC-CTRL/outs/singlecell.csv")
## colnames(single.csv)
dot.csv <- fread("/wsu/home/groups/piquelab/scGxE/counts_cellranger_arc/scGxE_1-1/outs/per_barcode_metrics.csv")
colnames(dot.csv)
head(dot.csv$barcode)



combined$pct_reads_in_peaks <- combined$atac_peak_region_fragments/combined$atac_fragments*100
colnames(combined@meta.data)
head(combined@meta.data$pct_reads_in_peaks)

#combined$blacklist_ratio <- combined$blacklist_region_fragments/combined$peak_region_fragments


head(rownames(combined@meta.data))

opfn <- "./1_Merge.outs/1_seurat.merge.rds"
write_rds(combined, file=opfn)            

###
###
 

combined <- read_rds("./1_Merge.outs/1_seurat.merge.rds")

combined

colnames(combined@meta.data)





## ###############
## ### Summary ###
## ###############
## atac <- read_rds("./1_Merge.outs/1_seurat.merge.rds")
## atac$high.tss <- ifelse(atac$TSS.enrichment>2, "High", "Low")

## fig0 <- TSSPlot(atac, group.by="high.tss")+
##         NoLegend()+
##         theme(plot.title=element_text(hjust=0.5))        
## figfn <- "./1_Merge.outs/Figure1.tss.png"
## png(figfn, width=500, height=400, res=120)
## print(fig0)
## dev.off() 

## ##
## #atac$nucleosome_group <- ifelse(atac$nucleosome_signal>4, "NS > 4", "NS < 4")
## #fig1 <- FragmentHistogram(object=atac, group.by = "nucleosome_group")
## #figfn <- "./outs/Figure2.fragment.png"
## #png(figfn, width=500, height=400, res=120)
## #print(fig1)
## #dev.off()

## ###
## fig2 <- VlnPlot(object=atac,
##                 features=c("pct_reads_in_peaks", "peak_region_fragments", "TSS.enrichment", 
##                            "blacklist_ratio", "nucleosome_signal"),
##                 pt.size=0,
##                 ncol=3)&
##         theme_bw()+
##         theme(axis.title=element_blank(),
##               axis.text.x=element_blank(),
##               legend.position="none",
##               plot.title=element_text(hjust=0.5,size=10))
## figfn <- "./1_Merge.outs/Figure3.vlnplot.png"
## png(figfn, width=700, height=500, res=120)
## print(fig2)
## dev.off()





##########################################################################
##                      RNA                                             ##
##########################################################################
pbmc.combined <- readRDS(file = "/wsu/home/groups/piquelab/scGxE/merging_objects/pbmc_combined.rds", refhook = NULL)

pbmc.combined

pbmc.combined[["RNA"]]

pbmc.combined[["macs2"]]

colnames(pbmc.combined@meta.data)

opfn <- "./1_Merge.outs/1_seurat.merge.RNA.rds"
write_rds(pbmc.combined[["RNA"]], file=opfn)            



#----------------------------------------------------------
rna <- read_rds("./1_Merge.outs/1_seurat.merge.RNA.rds")

rna

#combined[["RNA"]] <- rna

Cells(rna) %in% Cells(combined)
head(Cells(rna))
head(Cells(combined))
table(Cells(combined))

rna <- RenameCells(rna, new.names=gsub("-1","",Cells(rna)))
rna <- RenameCells(rna, new.names=gsub("11","scGxE_1",Cells(rna)))
rna <- RenameCells(rna, new.names=gsub("12","scGxE_1-2",Cells(rna)))
rna <- RenameCells(rna, new.names=gsub("13","scGxE_1-3",Cells(rna)))
rna <- RenameCells(rna, new.names=gsub("14","scGxE_1-4",Cells(rna)))
rna <- RenameCells(rna, new.names=gsub("15","scGxE_1-5",Cells(rna)))
rna <- RenameCells(rna, new.names=gsub("16","scGxE_1-6",Cells(rna)))
rna <- RenameCells(rna, new.names=gsub("17","scGxE_1-7",Cells(rna)))
rna <- RenameCells(rna, new.names=gsub("18","scGxE_1-8",Cells(rna)))

opfn <- "./1_Merge.outs/1_seurat.merge.RNA.rds"
write_rds(rna, file=opfn)            

rna

combined[["RNA"]] <- rna
colnames(combined@meta.data)

#------------------------------------------------------------
pbmc.combined <- RenameCells(pbmc.combined, new.names=gsub("-1","",Cells(pbmc.combined)))
pbmc.combined <- RenameCells(pbmc.combined, new.names=gsub("11","scGxE_1",Cells(pbmc.combined)))
pbmc.combined <- RenameCells(pbmc.combined, new.names=gsub("12","scGxE_1-2",Cells(pbmc.combined)))
pbmc.combined <- RenameCells(pbmc.combined, new.names=gsub("13","scGxE_1-3",Cells(pbmc.combined)))
pbmc.combined <- RenameCells(pbmc.combined, new.names=gsub("14","scGxE_1-4",Cells(pbmc.combined)))
pbmc.combined <- RenameCells(pbmc.combined, new.names=gsub("15","scGxE_1-5",Cells(pbmc.combined)))
pbmc.combined <- RenameCells(pbmc.combined, new.names=gsub("16","scGxE_1-6",Cells(pbmc.combined)))
pbmc.combined <- RenameCells(pbmc.combined, new.names=gsub("17","scGxE_1-7",Cells(pbmc.combined)))
pbmc.combined <- RenameCells(pbmc.combined, new.names=gsub("18","scGxE_1-8",Cells(pbmc.combined)))



x <- pbmc.combined@meta.data

colnames(x)
head(rownames(x))
head(x)

x2 <- x %>% dplyr::select(nucleosome_signal, nucleosome_percentile, TSS.enrichment, TSS.percentile)

head(x2)

combined <- AddMetaData(combined, x2) 

colnames(combined@meta.data)
         
opfn <- "./1_Merge.outs/1_seurat.merge.combined.rds"
write_rds(combined, file=opfn)            



#############################################################################
# summary
#############################################################################



