##
library(tidyverse)
library(data.table)
library(Matrix)
library(openxlsx)
###
options(scipen=200)

outdir <- "./4_SNPAnnot.outs/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)



###
### read all snps
fn <- "gtex_v8_snpinfor.txt.gz"
allsnp <- fread(fn, header=F, data.table=F)
allsnp <- allsnp%>%mutate(chr_pos=paste(V1, V2, sep="_"))


###
### cell type active peak bed files
fns <- list.files("Peak_bed", pattern="_peak.bed.gz")
MCls <- gsub("SNP_|_peak.bed.gz", "", fns)%>%unique()%>%sort()

chr_pos <- lapply(MCls, function(ii){
   ### 
   fn <- paste("./Peak_bed/SNP_", ii, "_peak.bed.gz", sep="")
   peak <- fread(fn, header=F, data.table=F)
   peak2 <- peak%>%dplyr::filter(V4>0)%>%mutate(chr_pos=paste(V1, V2, sep="_"))
   chr_pos <- peak2$chr_pos
   chr_pos
})
chr_pos2 <- unique(unlist(chr_pos))
snpList <- allsnp%>%dplyr::filter(chr_pos%in%chr_pos2)%>%pull(V3)

### torus format
anno <- data.frame(SNP=allsnp$V3)%>%mutate(peaking_d=ifelse(SNP%in%snpList, 1, 0))

opfn <- gzfile(paste(outdir, "zzz_peaking_torus.annot.gz", sep=""))
write.table(anno, opfn, quote=F, row.names=F, col.names=T)


###
### response motifs
fn <- "./4_SNPAnnot.outs/zzz_peaking_torus.annot.gz"
anno <- fread(fn, header=T, data.table=F)

snpList <- anno%>%filter(peaking_d==1)%>%pull(SNP)%>%unique()
 
## fn <- "./Response_motif/2.2_response_motif_th0.1.txt"
fn <- "./Response_motif/2.3_response_motif_th0.2.txt"
resMotif <- read.table(fn, header=T)
summ <- resMotif%>%group_by(comb)%>%summarise(ny=n())%>%filter(ny>0)
comb2 <- sort(unique(summ$comb)) ## 49 conditions
comb2 <- comb2[!grepl("8_MAIT|water",comb2)] ## 40 conditions


### annotation
mat <- lapply(comb2, function(ii){
  ###
  motifList <- resMotif%>%filter(comb==ii)%>%pull(motif_ID)%>%unique()  
  chr_pos2 <- lapply(motifList, function(motif){
      ###
      ##cat(motif, "\n")
      fn <- paste("./annot_jaspar2022/allsnp_", motif, ".bed.gz", sep="")   
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

##
opfn <- gzfile(paste(outdir, "zzz_th0.2_multiConditions_torus.annot.gz", sep=""))
write.table(anno2, opfn, quote=F, row.names=F, col.names=T)



##################################################
### New annotation-peaks and peak bound by TFs ###
##################################################

### response motif
fn <- "./Response_motif/2.3_response_motif_th0.2.txt"
resMotif <- read.table(fn, header=T)%>%pull(motif_ID)%>%unique()


### all snps
fn <- "gtex_v8_snpinfor.txt.gz"
allsnp <- fread(fn, header=F, data.table=F)
allsnp <- allsnp%>%mutate(chr_pos=paste(V1, V2, sep="_"))


### annotation torus
fn <- "./4_SNPAnnot.outs/zzz_peaking_torus.annot.gz"
anno <- fread(fn, header=T, data.table=F)
snpList <- anno%>%filter(peaking_d==1)%>%pull(SNP)%>%unique() 

###
fn <- "./annot_jaspar2022/zzz_allmotif.bed.gz"
x <- fread(fn, header=F, data.table=F)
x <- x%>%mutate(chr_pos=paste(V1, V2, sep="_"))
x2 <- x%>%filter(V5%in%resMotif)

chr_pos2 <- unique(x2$chr_pos)
snpmotif <- allsnp%>%filter(chr_pos%in%chr_pos2)%>%pull(V3)%>%unique()

### motif binding snp in peaks
snp2  <- intersect(snpList, snpmotif)

anno2 <- anno%>%mutate(peaking_d=ifelse(SNP%in%snp2, 2, peaking_d))

opfn <- gzfile(paste(outdir, "zzz_th0.2_union_torus.annot.gz", sep=""))
write.table(anno2, opfn, quote=F, row.names=F, col.names=T)




##########################################
### summary cell type peak information ###
##########################################

fns <- list.files("Peak_bed", pattern="_peak.bed.gz")
MCls <- gsub("SNP_|_peak.bed.gz", "", fns)%>%unique()%>%sort()

 
summDF <- map_dfr(MCls, function(ii){
   ###
   cat(ii, "\n") 
   fn <- paste("./Peak_bed/SNP_", ii, "_peak.bed.gz", sep="")
   peak <- fread(fn, header=F, data.table=F)
   peak2 <- peak%>%dplyr::filter(V4>0)%>%mutate(chr_pos=paste(V1, V2, sep="_"))
   chr_pos2 <- peak2$chr_pos
   snps <- allsnp%>%dplyr::filter(chr_pos%in%chr_pos2)%>%pull(V3) 
   ###
   fn <- paste("./Peak_bed/", ii, "_peak.bed", sep="")
   peak <- read.table(fn, header=F)
   ###
   df2 <- data.frame(MCls=ii, npeak=nrow(peak), nsnps=length(snps)) 
})

opfn <- paste(outdir, "1.1_MCl_peaks.xlsx", sep="")
write.xlsx(summDF, file=opfn, overwrite=T)


################################
### summary response motifs  ###
################################
 
##fn <- "./Response_motif/2.2_response_motif_each.txt"
fn <- "./Response_motif/2_response_motif.txt"
resMotif <- read.table(fn, header=T)
summ <- resMotif%>%group_by(comb)%>%summarise(ny=n())%>%
    dplyr::filter(!grepl("8_MAIT", comb))%>%arrange(comb)
 
df2   <- data.frame(comb=gsub("_d$", "", colnames(mat)), nsnps=colSums(mat))

summ <- summ%>%left_join(df2, by="comb")%>%as.data.frame()
 
opfn <- paste(outdir, "2_response_motif.xlsx", sep="")
write.xlsx(summ, file=opfn, overwrite=T)


##
x <- read.xlsx("./4_SNPAnnot.outs/2_response_motif.xlsx")
names(x) <- c("comb", "R1_nmotif", "R1_nsnps")
##
x2 <- read.xlsx("./4_SNPAnnot.outs/2.2_response_motif.xlsx")
names(x2) <- c("comb", "R2_nmotif", "R2_nsnps")

x <- x%>%full_join(x2, by="comb")
x <- x%>%arrange(comb)

opfn <- "./4_SNPAnnot.outs/2.3_comb.xlsx"
write.xlsx(x, file=opfn, overwrite=T)






