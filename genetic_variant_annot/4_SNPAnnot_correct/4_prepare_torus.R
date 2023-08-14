##
library(tidyverse)
library(data.table)
library(Matrix)
library(openxlsx)
###
options(scipen=18)


### Here we use newly generated response motif files to prepare annotation files
### ../../analysis_multiome_JW_2023-06-03/4_motif_plots.outs/response_motif_correct/2.2_response_motif_th0.1.txt
### cp the file in this folder

###
### read all snps
fn <- "../gtex_v8_snpinfor.txt.gz"
allsnp <- fread(fn, header=F, data.table=F)
allsnp <- allsnp%>%mutate(chr_pos=paste(V1, V2, sep="_"))



###
### SNP annotation in peak  
fn <- "zzz_peaking_torus.annot.gz"
anno <- fread(fn, header=T, data.table=F)

snpList <- anno%>%filter(peaking_d==1)%>%pull(SNP)%>%unique()


## fn <- "./Response_motif/2.2_response_motif_th0.1.txt"
fn <- "2.2_response_motif_th0.1.txt"
resMotif <- read.table(fn, header=T)
summ <- resMotif%>%group_by(comb)%>%summarise(ny=n())%>%filter(ny>0)
comb2 <- sort(unique(summ$comb)) ## 35 conditions

 
### annotation
mat <- lapply(comb2, function(ii){
  ###
  motifList <- resMotif%>%filter(comb==ii)%>%pull(motif_ID)%>%unique()  
  chr_pos2 <- lapply(motifList, function(motif){
      ###
      ##cat(motif, "\n")
      fn <- paste("../annot_jaspar2022/allsnp_", motif, ".bed.gz", sep="")   
      x <- try(fread(fn, header=F, data.table=F, stringsAsFactors=F), silent=T)
      if ( class(x)!="try-error"){ 
         x <- x%>%mutate(chr_pos=paste(V1, V2, sep="_"))
         xx <- x$chr_pos
      }else{
         xx <- NA
      }   
      xx
   })
   ###snpmotif <- unlist(snpmotif[!is.na(snpmotif)])
   chr_pos2 <- unique(unlist(chr_pos2[!is.na(chr_pos2)]))
   snpmotif <- allsnp%>%filter(chr_pos%in%chr_pos2)%>%pull(V3)
   snp2 <- intersect(snpmotif, snpList)

   cat(ii, length(snp2), "\n") 
   ###
   ### 
   anno2 <- anno%>%mutate(anno_d=ifelse(SNP%in%snp2, 1, 0))
   dd <- anno2%>%pull(anno_d)
   dd
})

mat <- do.call(cbind, mat)
colnames(mat) <- paste(comb2, "_d", sep="")
###
###
anno2 <- cbind(anno, mat)


opfn <- "zzz_th0.1_multiConditions_torus.annot.gz"
fwrite(anno2, opfn, quote=F, sep=" ")



##################################################
### New annotation-peaks and peak bound by TFs ###
##################################################

### response motif
fn <- "2.2_response_motif_th0.1.txt"
resMotif <- read.table(fn, header=T)%>%pull(motif_ID)%>%unique()


### all snps
fn <- "../gtex_v8_snpinfor.txt.gz"
allsnp <- fread(fn, header=F, data.table=F)
allsnp <- allsnp%>%mutate(chr_pos=paste(V1, V2, sep="_"))


### annotation torus
fn <- "zzz_peaking_torus.annot.gz"
anno <- fread(fn, header=T, data.table=F)
snpList <- anno%>%filter(peaking_d==1)%>%pull(SNP)%>%unique() 

###
fn <- "../annot_jaspar2022/zzz_allmotif.bed.gz"
x <- fread(fn, header=F, data.table=F)
x <- x%>%mutate(chr_pos=paste(V1, V2, sep="_"))
x2 <- x%>%filter(V5%in%resMotif)

chr_pos2 <- unique(x2$chr_pos)
snpmotif <- allsnp%>%filter(chr_pos%in%chr_pos2)%>%pull(V3)%>%unique()

### motif binding snp in peaks
snp2  <- intersect(snpList, snpmotif)

anno2 <- anno%>%mutate(peaking_d=ifelse(SNP%in%snp2, 2, peaking_d))

opfn <- "zzz_th0.1_union_torus.annot.gz"
fwrite(anno2, opfn, quote=F, sep=" ")

## x <- fread(opfn, header=T, data.table=F)










