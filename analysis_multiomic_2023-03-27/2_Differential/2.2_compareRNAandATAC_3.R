###
library(Matrix)
library(tidyverse)
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
### annotation required package
library(clusterProfiler)
library(org.Hs.eg.db)
library(annotables)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(EnsDb.Hsapiens.v86)
library(ChIPseeker)
#library(ChIPseeker, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
###
library(ggplot2)
library(cowplot, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(grid)
library(gridExtra)
library(ggExtra)
library(RColorBrewer)
library(ggsci)
library(RColorBrewer)
library(viridis)
theme_set(theme_grey())
##
library(ComplexHeatmap)
library(circlize)

outdir <- "./2.2_compareRNAandATAC_3.outs/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)


###########################
### annotation database ###
###########################
###
### annotation database
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
edb <- EnsDb.Hsapiens.v75
seqlevelsStyle(edb) <- "UCSC"



###########################
### atac object ###
###########################
## atac <- read_rds("../1_processing/5.1_reCallPeak.outs/3_scATAC.annot.mitofilt.ndup.rds")

## anno <- Annotation(atac)
## head(anno)


## annotations <- GetGRangesFromEnsDb(ensdb=EnsDb.Hsapiens.v86)
## #seqlevelsStyle(annotations) <- "1000"
## seqlevelsStyle(annotations) <- "UCSC"
## genome(annotations) <- "hg38"
## Annotation(atac) <- annotations


## anno <- Annotation(atac)
## head(anno)








###########################################
### peak annotation from signac - peakAnno
###########################################
## granges(atac)


## peakAnno <- ClosestFeature(atac, regions=granges(atac))
## opfn <- "./2.2_compareRNAandATAC.outs/1_annot.signac_macs2_0.1_cn.rds"
## write_rds(peakAnno, opfn)


peakAnno <- read_rds("./2.2_compareRNAandATAC.outs/1_annot.signac_macs2_0.1_cn.rds")
peakAnno <- peakAnno %>% mutate(query_region=gsub(":", "-", query_region))
colnames(peakAnno)[which(names(peakAnno) == "query_region")] <- "peak"
head(peakAnno)
nrow(peakAnno)

peakAnno.select <- peakAnno%>%dplyr::select(peak, gene_id, gene_name, symbol, entrezid, distance)
head(peakAnno.select)
nrow(peakAnno.select)

peakAnno.filt <- peakAnno.select%>%dplyr::filter(distance<100000)%>% dplyr::filter(!is.na(symbol))
head(peakAnno.filt)
nrow(peakAnno.filt)






## ################################################
## ### peak annotation from ChIPseeker - edb
## ################################################
## #x <- as.data.frame(granges(atac))%>%
## #   mutate(seqnames=paste("chr",seqnames,sep=""))

## x <- as.data.frame(granges(atac))
## head(x)

## peak <- makeGRangesFromDataFrame(x)
## head(peak)

## edb <- annotatePeak(peak, tssRegion=c(-3000,3000), TxDb=edb,
##    addFlankGeneInfo=T, flankDistance=100000, annoDb="org.Hs.eg.db")

## opfn <- "./2.2_compareRNAandATAC.outs/2_annot.ChIPseeker_edb_macs2_0.1_cn.rds"
## write_rds(edb, opfn)


## edb <- read_rds("./2.2_compareRNAandATAC.outs/2_annot.ChIPseeker_edb_macs2_0.1_cn.rds") %>% as.data.frame() %>% mutate(peak = paste(seqnames, start, end, sep="-"))
## head(edb)
## nrow(edb)
## colnames(edb)

## edb.select <- edb %>% dplyr::select(peak, geneChr, geneId, GENENAME, SYMBOL, ENTREZID, distanceToTSS)
## head(edb.select)
## nrow(edb.select)
## #head(edb.select %>% dplyr::filter(distanceToTSS > 10000)) 


## edb.filt <- edb.select %>% dplyr::filter(distanceToTSS < 100000) %>% dplyr::filter(!is.na(SYMBOL))
## head(edb.filt)
## nrow(edb.filt)









## ################################################
## ### peak annotation from ChIPseeker - txdb
## ################################################
## opfn <- "./2.2_compareRNAandATAC.outs/2_annot.ChIPseeker_macs2_0.1_cn.rds"
## write_rds(peakAnno, opfn)

## txdb <- read_rds("./2.2_compareRNAandATAC.outs/2_annot.ChIPseeker_macs2_0.1_cn.rds") %>% as.data.frame() %>% mutate(peak = paste(seqnames, start, end, sep="-"))
## head(txdb)
## nrow(txdb)
## colnames(txdb)

## txdb.select <- txdb %>% dplyr::select(peak, geneChr, gene_id, gene_name, symbol, entrezid, distanceToTSS)
## head(txdb.select)
## #head(txdb.select %>% dplyr::filter(distanceToTSS > 10000)) 


## txdb.filt <- txdb.select %>% dplyr::filter(distanceToTSS < 100000) %>% dplyr::filter(!is.na(symbol))
## head(txdb.filt)
## nrow(txdb.filt)






## ####################################################
## ### peak annotattion using linkpeakstogenes - link
## ####################################################
## ##link <- read_rds("../1_processing/6.1_linkPeakstoGenes.outs/links_dist100000.rds")
link <- read_rds("../1_processing/6.1_linkPeakstoGenes.outs/links.rds")
## head(link)
## nrow(link)
## max(link$pvalue)


link.filt <- link %>% arrange(pvalue)
link.filt <- link.filt %>% distinct(peak, .keep_all = TRUE) %>% dplyr::filter(!is.na(gene))

head(link.filt)
nrow(link.filt)     

end






## ###################################
## ### compare signac and ChIPseeker
## ###################################
## ### compare the closest genes
## ## x1 <- read_rds("./2.2_compareRNAandATAC.outs/1_annot.signac.rds")
## ## x1 <- x1%>%dplyr::rename("peak_region"="query_region")%>%dplyr::select(peak_region, gene_id)  

## ## x2 <- read_rds("./2.2_compareRNAandATAC.outs/2_annot.ChIPseeker.rds")%>%
## ##     as.data.frame()%>%
## ##     mutate(chr=gsub("chr","",seqnames),
## ##            peak_region=paste(chr,start,end,sep="-"))%>%
## ##     dplyr::select(peak_region, geneId)

## ## x <- x1%>%left_join(x2, by="peak_region")


## peakAnno.txdb <- left_join(peakAnno, txdb, by="peak")

## head(peakAnno.txdb)

## ## ### compare the genes around regions
## ## x2 <- read_rds("./2.2_compareRNAandATAC.outs/2_annot.ChIPseeker.rds")%>%
## ##    as.data.frame()%>%
## ##    mutate(chr=gsub("chr","",seqnames),
## ##           peak_region=paste(chr,start,end,sep="-"))%>%
## ##    dplyr::select(peak_region, flank_geneIds)
## ## x <- x1%>%left_join(x2, by="peak_region")

## ## geneList <- str_split(x$flank_geneIds, ";")
## ## ii <- map_dbl(1:length(geneList), function(i){
## ##   gene0 <- as.character(x[i, "gene_id"])
## ##   x <- ifelse(gene0%in%geneList[[i]], 1, 0)
## ##   x
## ## })


#################################################
### resDP ###
#################################################
resDP <- read_rds("../2_Differential/3_treatEffect_2_2.outs/resDP_control.rds") %>% as.data.frame()

#resDP <- resDP %>% drop_na(p.adjusted)
##%>% drop_na(p.value)
#head(resDP)
#nrow(resDP)


table(resDP$contrast)
table(resDP$MCls)


#nrow(resDP %>% dplyr::filter(contrast=="caffeine"))
#nrow(resDP %>% dplyr::filter(contrast=="nicotine"))
#nrow(resDP %>% dplyr::filter(contrast=="zinc"))
#nrow(resDP %>% dplyr::filter(contrast=="vitA"))
#nrow(resDP %>% dplyr::filter(contrast=="vitD"))
#nrow(resDP %>% dplyr::filter(contrast=="vitE"))


sig_resDP <- resDP %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)

head(sig_resDP)
nrow(sig_resDP)

peakAll <- unique(sig_resDP$peak)

end




#################################################
### resDP.peakAnno ###
#################################################
resDP.peakAnno <- left_join(resDP, peakAnno.filt, by="peak")

resDP.peakAnno <- resDP.peakAnno %>% mutate(contrast.gene=paste0(comb, ".", symbol)) %>% drop_na(symbol)

head(resDP.peakAnno)
nrow(resDP.peakAnno)
colnames(resDP.peakAnno)


write_rds(resDP.peakAnno, "./2.2_compareRNAandATAC_3.outs/resDP.peakAnno.rds")

#resDP.peakAnno.2 <- readRDS("./2.2_compareRNAandATAC_3.outs/resDP.peakAnno.rds")
#nrow(resDP.peakAnno.2)

sig_resDP.peakAnno <- resDP.peakAnno %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)

head(sig_resDP.peakAnno)
nrow(sig_resDP.peakAnno)

#sig_resDP.peakAnno.2 <- resDP.peakAnno.2 %>% dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)
##head(sig_resDP.peakAnno.2)
#nrow(sig_resDP.peakAnno.2)


#unique(sig_resDP.peakAnno$comb)



#################################################
### resDP.peakAnno.sum ###
#################################################
comb <- unique(resDP.peakAnno$comb) ## this can be significant data
comb

resDP.pa.sum <- map_dfr(comb, function(ii){
    res.temp <- resDP.peakAnno %>% dplyr::filter(comb==ii)
    res.temp <- res.temp %>% arrange(p.adjusted)
    res.temp <- res.temp %>% distinct(symbol, .keep_all = TRUE)
})

head(resDP.pa.sum)

nrow(resDP.pa.sum)

#comb <- unique(resDP.pa.sum$comb) ## this can be significant data
#comb


sig_resDP.pa.sum <- resDP.pa.sum %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)
head(sig_resDP.pa.sum)

nrow(sig_resDP.pa.sum)

#comb <- unique(sig_resDP.pa.sum$comb)
#comb







## #################################################
## ### resDP.txdb ###
## #################################################
## resDP.txdb <- left_join(resDP, txdb.filt, by="peak") 

## resDP.txdb <- resDP.txdb %>% mutate(contrast.gene=paste0(comb, ".", symbol)) %>% drop_na(symbol)

## head(resDP.txdb)
## nrow(resDP.txdb)


## #################################################
## ### resDP.edb ###
## #################################################
## resDP.edb <- left_join(resDP, edb.filt, by="peak")


## resDP.edb <- resDP.edb %>% mutate(contrast.gene=paste0(comb, ".", SYMBOL)) %>% drop_na(SYMBOL)

## head(resDP.edb)
## nrow(resDP.edb)


## #################################################
## ### resDP.link ###
## #################################################
resDP.link <- left_join(resDP, link.filt, by="peak")

resDP.link <- resDP.link %>% mutate(contrast.gene=paste0(comb, ".", gene)) %>% drop_na(gene)

head(resDP.link)
nrow(resDP.link)


sig_resDP.link <- resDP.link %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)
#head(sig_resDP.link)

nrow(sig_resDP.link)




end









## #################################################
## ### ##----mat data for significant----###
## #################################################
## comb <- unique(resDP.pa.sum$comb) ## this can be significant data
## comb

## sig.genes <- unique(sig_resDP.pa.sum$symbol)
## head(sig.genes)
## length(sig.genes)

## df <- resDP.pa.sum %>% dplyr::select(estimate, statistic)
## head(df)

## ##for (i in c("estimate", "statistic")){
## ##for (i in colnames(df)){
##     ##print(i)
##     ##print(noquote(i))
##     ##print(noquote(paste0("res2$", i)))
##     mat <- map_dfc(comb, function(ii){
##         res2 <- resDP.pa.sum %>% dplyr::filter(comb==ii)
##         b <- res2$estimate
##         print(head(b))
##         names(b) <- res2$symbol
##         b[sig.genes]
##     })
##     nrow(mat)
##     ##
##     mat <- as.matrix(mat)
##     rownames(mat) <- sig.genes 
##     colnames(mat) <- comb

##     head(mat)

##     mat[is.na(mat)] = 0
##     max(mat)
##     min(mat)
##     mat2 <- mat
##     comb.2 <- sort(comb)
##     mat3 <- mat2[, comb.2]

## head(mat3)

##     ##----DP sig plot Fig1.1----
##     ##breaks and color scale
##     y <- as.numeric(mat2)
##     max(y)
##     min(y)
##     scale <- quantile(abs(y),probs=0.99)
##     scale
##     ##y0 <- y[abs(y)< quantile(abs(y),probs=0.99)] #99% percent quantile(abs(y),probs=0.99)
##     ##mybreaks <- c(min(y),quantile(y0,probs=seq(0,1,length.out=98)),max(y))
##     mybreaks <- seq(-scale, scale, length.out=100)
##     names(mybreaks) <- NULL
##     mybreaks
##     mycol <- colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100)
##     col1 <- c("caffeine"="red", "nicotine"="tan",
##               "vitA"="tan4", "vitD"="seagreen4",
##               "vitE"="salmon3", "water"="grey", "zinc"="maroon3")
##     col2 <- c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
##               "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00", "1-TCM"="pink",
##               "3-TEM"="blue", "5-CD8Naive"="green", "7-dnT"="black")
##     tmp_colors <- list(treatment=col1, celltype=col2) #brewer.pal(4,"Set1")

##     ##
##     x <- str_split(comb, "_", simplify=T)
##     tmp_column <- data.frame(treatment=x[,1], celltype=x[,2])
##     rownames(tmp_column) <- comb
##     fig1 <- pheatmap(mat2,  col=mycol, breaks=mybreaks, border_color="NA",
##                      cluster_rows=T, cluster_cols=F,
##                      annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend=T,
##                      show_colnames=T,
##                      show_rownames=F,
##                      treeheight_row = 0,
##                      na_col="white",
##                      fontsize_col=7,
##                      fontsize_row=3)

##     figfn <- paste0("./2.2_compareRNAandATAC_3.outs/Figure1.1_heatmap.beta.DP.sig.cell_", "estimate", "_control_quantile_0.99.png")
##     png(figfn, width=1800, height=2100,res=225)
##     print(fig1)
##     dev.off()

##     ##
##     x <- str_split(comb.2, "_", simplify=T)
##     tmp_column <- data.frame(treatment=x[,1], celltype=x[,2])
##     rownames(tmp_column) <- comb.2
##     fig1 <- pheatmap(mat3,  col=mycol, breaks=mybreaks, border_color="NA",
##                      cluster_rows=T, cluster_cols=F,
##                      annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend=T,
##                      show_colnames=T,
##                      show_rownames=F,
##                      treeheight_row = 0,
##                      na_col="white",
##                      fontsize_col=7,
##                      fontsize_row=3)
##     figfn <- paste0("./2.2_compareRNAandATAC_3.outs/Figure1.1_heatmap.beta.DP.sig.treat_", "estimate", "_control_quantile_0.99.png")
##     png(figfn, width=1800, height=2100,res=225)
##     print(fig1)
##     dev.off()
## ##}








## #################################################
## #----mat data for for top 10 significant only----
## #################################################
## top10 <- map_dfr(comb, function(ii){
##         res2 <- sig_resDP.pa.sum %>% dplyr::filter(comb==ii)
##         res2 <- res2 %>% arrange(p.adjusted)
##         res2 <- res2 [1:10, ]
##         print(ii)
##         print(res2$symbol)
## })

## head(top10)
## nrow(top10)

## write.csv(top10, "./2.2_compareRNAandATAC_3.outs/DP.top10_control.csv")


## print("Number of top 10 significant")
## nrow(top10)

## ##
## sig.top10.genes <- unique(top10$symbol)
## head(sig.top10.genes)
## length(sig.top10.genes)

## ##
## mat <- map_dfc(comb, function(ii){
##     res2 <- resDP.pa.sum %>% dplyr::filter(comb==ii)
##     b <- res2$statistic
##     names(b) <- res2$symbol
##     b[sig.top10.genes]
##     ##print(ii)
##     ##print(b)
## })
## nrow(mat)

## ##
## mat <- as.matrix(mat)
## rownames(mat) <- sig.top10.genes
## colnames(mat) <- comb
## head(mat)
## mat[is.na(mat)] = 0
## max(mat)
## min(mat)
## mat2 <- mat
## comb.2 <- sort(comb)
## mat3 <- mat2[, comb.2]

## mat3.hc <- hclust(mat3, method = "ward.D2")
## print(rownames(mat3))

## #----top10 DP plot Fig1.1----
## #breaks and color scale
## y <- as.numeric(mat2)
## max(y)
## abs(min(y))
## scale <- quantile(abs(y),probs=0.99)
## scale
## ##y0 <- y[abs(y)< quantile(abs(y),probs=0.99)] #99% percent quantile(abs(y),probs=0.99)
## ##mybreaks <- c(min(y),quantile(y0,probs=seq(0,1,length.out=98)),max(y))
## mybreaks <- seq(-scale, scale, length.out=100)
## names(mybreaks) <- NULL
## mybreaks

## mycol <- colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100)
## col1 <- c("caffeine"="red", "nicotine"="tan",
##           "vitA"="tan4", "vitD"="seagreen4",
##           "vitE"="salmon3", "water"="grey", "zinc"="maroon3")
## col2 <- c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
##           "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00", "1-TCM"="pink",
##           "3-TEM"="blue", "5-CD8Naive"="green", "7-dnT"="black")
## tmp_colors <- list(treatment=col1, celltype=col2) #brewer.pal(4,"Set1")

## ##
## x <- str_split(comb, "_", simplify=T)
## tmp_column <- data.frame(treatment=x[,1], celltype=x[,2])
## rownames(tmp_column) <- comb
## fig1 <- pheatmap(mat2,  col=mycol, breaks=mybreaks, border_color="NA",
##                  cluster_rows=T, cluster_cols=F, treeheight_row = 0,
##                  annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend=T,
##                  show_colnames=T,
##                  show_rownames=T,
##                  na_col="white",
##                  fontsize_col=7,
##                  fontsize_row=2)

## ##figfn <- paste0("./2.2_compareRNAandATAC_3.outs/Figure1.2_heatmap.beta.DP.top10.cell_", "estimate", "_control_quantile_0.99.png")
## figfn <- paste0("./2.2_compareRNAandATAC_3.outs/Figure1.2_heatmap.beta.DP.top10.cell_", "zscore", "_control_quantile_0.99.png")
## png(figfn, width=2500, height=3000,res=225)
## print(fig1)
## dev.off()

## ##
## x <- str_split(comb.2, "_", simplify=T)
## tmp_column <- data.frame(treatment=x[,1], celltype=x[,2])
## rownames(tmp_column) <- comb.2

## fig1 <- pheatmap(mat3,  col=mycol, breaks=mybreaks, border_color="NA",
##                  cluster_rows=T, cluster_cols=F,treeheight_row = 0,
##                  annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend=T,
##                  show_colnames=T,
##                  show_rownames=T,
##                  na_col="white",
##                  fontsize_col=7,
##                  fontsize_row=2)

## ##figfn <- paste0("./2.2_compareRNAandATAC_3.outs/Figure1.2_heatmap.beta.DP.top10.treat_", "estimate", "_control_quantile_0.99.png")
## figfn <- paste0("./2.2_compareRNAandATAC_3.outs/Figure1.2_heatmap.beta.DP.top10.treat_", "zscore", "_control_quantile_0.99_highres.png")
## png(figfn, width=2500, height=3500,res=225)
## print(fig1)
## dev.off()


## end























#################################################
### resDE ###
#################################################
resDE <- read_rds("../2_Differential/3_treatEffect_2_2.outs/resDE_control.mitofilt.rds") %>% as.data.frame()

#resDE <- resDE %>% drop_na(p.adjusted) %>% mutate(contrast.gene=paste0(comb, ".", gene))
#%>% drop_na(p.value)
#head(resDE)
#nrow(resDE)

sig_resDE <- resDE %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)
head(sig_resDE)
nrow(sig_resDE)

geneAll <- unique(sig_resDE$gene)






#################################################
### CORRELATION
#################################################
MCls <- unique(resDE$MCls)
MCls
comb.2 <- MCls

treats <- unique(resDE$contrast)

treats <- c("caffeine", "nicotine", "vitA", "vitD", "vitE", "zinc")
treats


#----all----
logFC <- "NA"
include <- "all"


#----Significance.0----
logFC <- 0
conditions <- c("both", "sc", "GXE")
include <- "union"


## alt_func <- function(j){
##         if (j=="caffeine"){
##             #scorr_caff <- c(scorr_caff, "NA")
##             scorr_caff_z <- c(scorr_caff_z, 0)
##         } else if (j=="nicotine"){
##             #scorr_nic <- c(scorr_nic, "NA")
##             scorr_nic_z <- c(scorr_nic_z, 0)
##         } else if (j=="vitA"){
##             #scorr_vitA <- c(scorr_vitA, "NA")
##             scorr_vitA_z <- c(scorr_vitA_z, 0)
##         } else if (j=="vitD"){
##             #scorr_vitD <- c(scorr_vitD, "NA")
##             scorr_vitD_z <- c(scorr_vitD_z, 0)
##         } else if (j=="vitE"){
##             #scorr_vitE <- c(scorr_vitE, "NA")
##             scorr_vitE_z <- c(scorr_vitE_z, 0)
##         } else if (j=="water"){
##             #scorr_water <- c(scorr_water, "NA")
##             scorr_water_z <- c(scorr_water, 0)
##         } else if (j=="zinc"){
##             #scorr_zinc <- c(scorr_zinc, "NA")
##             scorr_zinc_z <- c(scorr_zinc_z, 0)
##         }
## }


scorr_caff <- c()
scorr_nic <- c()
scorr_vitA <- c()
scorr_vitD <- c()
scorr_vitE <- c()
scorr_zinc <- c()
scorr_caff_z <- c()
scorr_nic_z <- c()
scorr_vitA_z <- c()
scorr_vitD_z <- c()
scorr_vitE_z <- c()
scorr_water_z <- c()
scorr_zinc_z <- c()

head(resDP.link)

head(resDP.peakAnno)


for (j in treats){
    resDE2 <- resDE %>% dplyr::filter(contrast==j)
    ##
    resDP2 <- resDP.peakAnno %>% dplyr::filter(contrast==j)
    ##
    ##resDP2 <- resDP.link %>% dplyr::filter(contrast==j)
    ##
    resDP2.sum <- resDP2 %>% arrange(p.adjusted)
    resDP2.sum <- resDP2.sum %>% distinct(symbol, .keep_all = TRUE) #for peakAnno
    ##resDP2.sum <- resDP2.sum %>% distinct(gene, .keep_all = TRUE) #for link
    ##
    resDE.DP <- left_join(resDE2, resDP2.sum, by="contrast.gene")
    resDE.DP <- resDE.DP %>% dplyr::filter(!is.na(estimate.y))
    ##
    resDE.DP <- resDE.DP %>% mutate(Significance=ifelse(Significance.y %in% c("Up", "Down") &  Significance.x %in% c("Up", "Down"), "both", ifelse(Significance.y %in% c("Up", "Down") & (Significance.x %in% c("Not Significant", "NA") | is.na(Significance.x)), "DP", ifelse((Significance.y %in% c("Not Significant", "NA") | is.na(Significance.y)) & Significance.x %in% c("Up", "Down"), "DG", "NA"))))
    ##
    resDE.DP <- resDE.DP %>% mutate(Significance.25=ifelse(Significance.25.y %in% c("Up", "Down") &  Significance.25.x %in% c("Up", "Down"), "both", ifelse(Significance.25.y %in% c("Up", "Down") & (Significance.25.x %in% c("Not Significant", "NA") | is.na(Significance.25.x)), "DP", ifelse((Significance.25.y %in% c("Not Significant", "NA") | is.na(Significance.25.y)) & Significance.25.x %in% c("Up", "Down"), "DG", "NA"))))
    ##
    resDE.DP <- resDE.DP %>% mutate(Significance.0=ifelse(Significance.0.y %in% c("Up", "Down") &  Significance.0.x %in% c("Up", "Down"), "both", ifelse(Significance.0.y %in% c("Up", "Down") & (Significance.0.x %in% c("Not Significant", "NA") | is.na(Significance.0.x)), "DP", ifelse((Significance.0.y %in% c("Not Significant", "NA") | is.na(Significance.0.y)) & Significance.0.x %in% c("Up", "Down"), "DG", "NA"))))
    ####################
    ## scorr_values ##
    ####################
    for (i in comb.2){
        tryCatch({
        res1 <- resDE.DP %>% dplyr::filter(contrast.x==j, MCls.x==i)       #                          %>% dplyr::filter(Significance.0 %in% conditions)
        #scorr <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$estimate)
        scorr_z <- as.numeric(cor.test(res1$statistic.x, res1$statistic.y, method = "spearman")$estimate)
        #print(scorr)
        print(scorr_z)
        if (is.numeric(scorr_z) & j=="caffeine"){
            #scorr_caff <- c(scorr_caff, scorr)
            scorr_caff_z <- c(scorr_caff_z, scorr_z)
        } else if (is.numeric(scorr_z) & j=="nicotine"){
            #scorr_nic <- c(scorr_nic, scorr)
            scorr_nic_z <- c(scorr_nic_z, scorr_z)
        } else if (is.numeric(scorr_z) & j=="vitA"){
            #scorr_vitA <- c(scorr_vitA, scorr)
            scorr_vitA_z <- c(scorr_vitA_z, scorr_z)
        } else if (is.numeric(scorr_z) & j=="vitD"){
            #scorr_vitD <- c(scorr_vitD, scorr)
            scorr_vitD_z <- c(scorr_vitD_z, scorr_z)
        } else if (is.numeric(scorr_z) & j=="vitE"){
            #scorr_vitE <- c(scorr_vitE, scorr)
            scorr_vitE_z <- c(scorr_vitE_z, scorr_z)
        } else if (is.numeric(scorr_z) & j=="water"){
            #scorr_water <- c(scorr_water, scorr)
            scorr_water_z <- c(scorr_water, scorr_z)
        } else if (is.numeric(scorr_z) & j=="zinc"){
            #scorr_zinc <- c(scorr_zinc, scorr)
            scorr_zinc_z <- c(scorr_zinc_z, scorr_z)
        }
        ##
        },  error=function(j){
            if (j=="caffeine"){
            #scorr_caff <- c(scorr_caff, "NA")
            scorr_caff_z <- c(scorr_caff_z, 0)
        } else if (j=="nicotine"){
            #scorr_nic <- c(scorr_nic, "NA")
            scorr_nic_z <- c(scorr_nic_z, 0)
        } else if (j=="vitA"){
            #scorr_vitA <- c(scorr_vitA, "NA")
            scorr_vitA_z <- c(scorr_vitA_z, 0)
        } else if (j=="vitD"){
            #scorr_vitD <- c(scorr_vitD, "NA")
            scorr_vitD_z <- c(scorr_vitD_z, 0)
        } else if (j=="vitE"){
            #scorr_vitE <- c(scorr_vitE, "NA")
            scorr_vitE_z <- c(scorr_vitE_z, 0)
        } else if (j=="water"){
            #scorr_water <- c(scorr_water, "NA")
            scorr_water_z <- c(scorr_water, 0)
        } else if (j=="zinc"){
            #scorr_zinc <- c(scorr_zinc, "NA")
            scorr_zinc_z <- c(scorr_zinc_z, 0)
        }
        })    
    }
}


#cat("ERROR", "\n")

scorr_caff_z <- c(scorr_caff_z, 0)
scorr_nic_z <- c(scorr_nic_z, 0)
scorr_zinc_z <- c(scorr_zinc_z, 0)
scorr_vitD_z <- c(scorr_vitD_z, 0)


scorr_caff_z
scorr_nic_z
scorr_vitA_z
scorr_vitD_z
scorr_vitE_z
scorr_water_z
scorr_zinc_z






library(ggsci)
library(circlize)
library(reshape2)
library("ggpubr")


scorr_data <- data.frame(MCls, scorr_caff_z, scorr_nic_z,
                         scorr_vitA_z, scorr_vitD_z,
                         scorr_vitE_z, scorr_zinc_z)
is.num <- sapply(scorr_data, is.numeric)
scorr_data[is.num] <- lapply(scorr_data[is.num], round, 2)
##
sig.melt <- melt(scorr_data)
##
sig.p <- ggplot(sig.melt, aes(x = variable, y = MCls)) +
    geom_tile(aes(fill=value), color="black")+
    scale_fill_gradient(low = "white", high = "white")+
    geom_text(aes(label = value), color = "black", size = 5) +
    ggtitle(paste0("logFC_", logFC, "_", include, "_zscorevszscore_corr_DG_vs_DP")) +
    theme(axis.text.x = element_text(angle = 90),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position="none")
##
figure <- ggarrange(sig.p + font("x.text", size = 14))
png(paste0("./2.2_compareRNAandATAC_3.outs/logFC_", logFC, "_", include, "_zscorevszscore_corr_link.png"), width=500, height=500, pointsize=16, res=125)
print(figure)
dev.off()


end

## checking if gene_name and symbol are same
for (i in (sig_resDP.peakAnno$gene_name==sig_resDP.peakAnno$symbol)){
    if(i==FALSE){
    print(sig_resDP.peakAnno$gene_name)
    print(sig_resDP.peakAnno$symbol)
    } 
}








#################################################
### SCATTER PLOT 
#################################################
library(ggplot2)
library(ggpmisc)

comb

head(sig_resDP.peakAnno)

#head(sig_resDP.link)

sig_comb.filt <- map_dfr(comb, function(ii){
    resDP2 <- sig_resDP.peakAnno%>%dplyr::filter(comb==ii) #peakAnno
    ##resDP2 <- sig_resDP.link%>%dplyr::filter(comb==ii) #link
    resDE2 <- sig_resDE%>%dplyr::filter(comb==ii)
    ##DP <- unique(resDP2$peak)
    DEG <- unique(resDE2$gene)
    ##x <- peakAnno2%>%dplyr::filter(peak%in%DP, gene_name%in%DEG)
    x <- resDP2 %>% dplyr::filter(gene_name%in%DEG)
    if(nrow(x)==0){
        res1 <- NULL
        ##next
    }else{
        ##res1 <- unique(x$gene_name)
        ##res1 <- x
        res1 <- x%>%mutate(gene=gene_name)
        res1 <- res1%>%left_join(resDE2, by="gene")
    }
    ##res <- list(DEG=res1)
    res1
})

head(sig_comb.filt)
nrow(sig_comb.filt)
length(unique(sig_comb.filt$gene))
table(sig_comb.filt$comb.x)
max(sig_comb.filt$estimate.y)
min(sig_comb.filt$estimate.y)
max(sig_comb.filt$estimate.x)
min(sig_comb.filt$estimate.x)






##----plot loop----
#install.packages("ggpmisc")
library(ggpmisc)

typeofDP <- "filt"
comb.2 <- unique(sig_comb.filt$contrast.x)
comb.2

comb.3 <- unique(sig_comb.filt$MCls.x)
comb.3


##----general----
## for(i in comb.2){
##     res1 <- sig_comb.filt %>% dplyr::filter(contrast.x==i)
##     ##res1 <- sig_comb.filt %>% dplyr::filter(MCls.x==i)
##     png(paste0("./2.2_compareRNAandATAC_2.outs/Figure1.1_significant_scatter_macs2_0.1_cn_", typeofDP, "_", i, ".png"))
##     p <- ggplot(res1, aes(x=estimate.y, y=estimate.x, color=contrast.x)) +
##         geom_point()+
##         xlab("sig_genes") + ylab("sig_peaks")+
##         xlim(-4, 10) + ylim(-4, 5)+
##         stat_poly_line() +
##         stat_poly_eq() +
##         ggtitle("fold.enrichment")
##     print(p)
##     dev.off()
## }



##----celltype----
## figfn <- paste0("./2.2_compareRNAandATAC_3.outs/Figure1.1_significant_scatter_macs2_0.1_cn_peakAnno_control_celltype.png")
## png(figfn, width=6000, height=8000, pointsize=33, res=250)
## par(mar=c(4,4,2,2),mgp=c(2,1,0))
## x <- matrix(1:8, 2, 4, byrow=T)
## layout(x)
## ##for (i in c(1)){
res1 <- sig_comb.filt %>% dplyr::filter(MCls.x==comb.3[1])
p1 <- ggplot(res1, aes(x=estimate.y, y=estimate.x, color=MCls.x), show.legend = FALSE) +
    geom_point()+
    xlab("sig_genes") + ylab("sig_peaks")+
    xlim(-4, 10) + ylim(-4, 5)+
    ##stat_poly_line() +
    ##stat_poly_eq() +
    ggtitle(comb.3[1])
p1 + theme(legend.position = "none")
scorr_est_1 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$estimate)
scorr_pval_1 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$p.value)
##print(p)
##
res1 <- sig_comb.filt %>% dplyr::filter(MCls.x==comb.3[2])
p2 <- ggplot(res1, aes(x=estimate.y, y=estimate.x, color=MCls.x), show.legend = FALSE) +
    geom_point()+
    xlab("sig_genes") + ylab("sig_peaks")+
    xlim(-4, 10) + ylim(-4, 5)+
    ##stat_poly_line() +
    ##stat_poly_eq() +
    ggtitle(comb.3[2])
p2 + theme(legend.position = "none")
scorr_est_2 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$estimate)
scorr_pval_2 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$p.value)
##print(p)
##
res1 <- sig_comb.filt %>% dplyr::filter(MCls.x==comb.3[3])
p3 <- ggplot(res1, aes(x=estimate.y, y=estimate.x, color=MCls.x), show.legend = FALSE) +
    geom_point()+
    xlab("sig_genes") + ylab("sig_peaks")+
    xlim(-4, 10) + ylim(-4, 5)+
    ##stat_poly_line() +
    ##stat_poly_eq() +
    ggtitle(comb.3[3])
p3 + theme(legend.position = "none")
scorr_est_3 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$estimate)
scorr_pval_3 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$p.value)
##print(p)
##
res1 <- sig_comb.filt %>% dplyr::filter(MCls.x==comb.3[4])
p4 <- ggplot(res1, aes(x=estimate.y, y=estimate.x, color=MCls.x), show.legend = FALSE) +
    geom_point()+
    xlab("sig_genes") + ylab("sig_peaks")+
    xlim(-4, 10) + ylim(-4, 5)+
    ##stat_poly_line() +
    ##stat_poly_eq() +
    ggtitle(comb.3[4])
p4 + theme(legend.position = "none")
scorr_est_4 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$estimate)
scorr_pval_4 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$p.value)
##print(p)
##
res1 <- sig_comb.filt %>% dplyr::filter(MCls.x==comb.3[5])
p5 <- ggplot(res1, aes(x=estimate.y, y=estimate.x, color=MCls.x), show.legend = FALSE) +
    geom_point()+
    xlab("sig_genes") + ylab("sig_peaks")+
    xlim(-4, 10) + ylim(-4, 5)+
    ##stat_poly_line() +
    ##stat_poly_eq() +
    ggtitle(comb.3[5])
p5 + theme(legend.position = "none")
scorr_est_5 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$estimate)
scorr_pval_5 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$p.value)
##print(p)
##
res1 <- sig_comb.filt %>% dplyr::filter(MCls.x==comb.3[6])
p6 <- ggplot(res1, aes(x=estimate.y, y=estimate.x, color=MCls.x), show.legend = FALSE) +
    geom_point()+
    xlab("sig_genes") + ylab("sig_peaks")+
    xlim(-4, 10) + ylim(-4, 5)+
    ##stat_poly_line() +
    ##stat_poly_eq() +
    ggtitle(comb.3[6])
p6 + theme(legend.position = "none")
scorr_est_6 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$estimate)
scorr_pval_6 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$p.value)
##print(p)
##
res1 <- sig_comb.filt %>% dplyr::filter(MCls.x==comb.3[7])
p7 <- ggplot(res1, aes(x=estimate.y, y=estimate.x, color=MCls.x), show.legend = FALSE) +
    geom_point()+
    xlab("sig_genes") + ylab("sig_peaks")+
    xlim(-4, 10) + ylim(-4, 5)+
    ##stat_poly_line() +
    ##stat_poly_eq() +
    ggtitle(comb.3[7])
p7 + theme(legend.position = "none")
scorr_est_7 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$estimate)
scorr_pval_7 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$p.value)
##print(p)
##
res1 <- sig_comb.filt %>% dplyr::filter(MCls.x==comb.3[8])
p8 <- ggplot(res1, aes(x=estimate.y, y=estimate.x, color=MCls.x), show.legend = FALSE) +
    geom_point()+
    xlab("sig_genes") + ylab("sig_peaks")+
    xlim(-4, 10) + ylim(-4, 5)+
    ##stat_poly_line() +
    ##stat_poly_eq() +
    ggtitle(comb.3[8])
p8 + theme(legend.position = "none")
scorr_est_8 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$estimate)
scorr_pval_8 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$p.value)
##print(p)
##
##print(mtext(1, side=4, line=0.5, cex=1, col="blue"))
##}
##dev.off()



print(scorr_est_1)
print(scorr_est_2)
print(scorr_est_3)
print(scorr_est_4)
print(scorr_est_5)
print(scorr_est_6)
print(scorr_est_7)
print(scorr_est_8)
print(scorr_pval_1)
print(scorr_pval_2)
print(scorr_pval_3)
print(scorr_pval_4)
print(scorr_pval_5)
print(scorr_pval_6)
print(scorr_pval_7)
print(scorr_pval_8)



library("ggpubr")
figure <- ggarrange(p1, p2, p3, p4, p5, p6, p7, p8 +
                                             font("x.text", size = 18),
                    ncol = 4, nrow = 2)
png(paste0("./2.2_compareRNAandATAC_3.outs/Figure1.1_significant_scatter_macs2_0.1_cn_peakAnno_control_celltype_noline.png"), width=4000, height=2000, pointsize=20, res=175)
print(figure)
dev.off()




##----treatment----
## figfn <- paste0("./2.2_compareRNAandATAC_3.outs/Figure1.1_significant_scatter_macs2_0.1_cn_peakAnno_control_treat.png")
## png(figfn, width=6000, height=8000, pointsize=33, res=250)
## par(mar=c(4,4,2,2),mgp=c(2,1,0))
## x <- matrix(1:5, 1, 5, byrow=T)
## layout(x)
##print(layout(x))
##for (i in comb.2){
res1 <- sig_comb.filt %>% dplyr::filter(contrast.x=="caffeine")
p1 <- ggplot(res1, aes(x=estimate.y, y=estimate.x, color=contrast.x), show.legend = FALSE) +
    geom_point()+
    xlab("sig_genes") + ylab("sig_peaks")+
    xlim(-4, 10) + ylim(-4, 5)+
    ##stat_poly_line() +
    ##stat_poly_eq() +
    ggtitle("caffeine")
p1 + theme(legend.position = "none")
p1 + guides(color = FALSE)
scorr_est_1 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$estimate)
scorr_pval_1 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$p.value)
##print(p)
##
res1 <- sig_comb.filt %>% dplyr::filter(contrast.x=="nicotine")
p2 <- ggplot(res1, aes(x=estimate.y, y=estimate.x, color=contrast.x), show.legend = FALSE) +
    geom_point()+
    xlab("sig_genes") + ylab("sig_peaks")+
    xlim(-4, 10) + ylim(-4, 5)+
    ##stat_poly_line() +
    ##stat_poly_eq() +
    ggtitle("nicotine")
p2 + theme(legend.position = "none")
p2 + guides(color = FALSE)
scorr_est_2 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$estimate)
scorr_pval_2 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$p.value)
##print(p)
##
res1 <- sig_comb.filt %>% dplyr::filter(contrast.x=="vitA")
p3 <- ggplot(res1, aes(x=estimate.y, y=estimate.x, color=contrast.x), show.legend = FALSE) +
    geom_point()+
    xlab("sig_genes") + ylab("sig_peaks")+
    xlim(-4, 10) + ylim(-4, 5)+
    ##stat_poly_line() +
    ##stat_poly_eq() +
    ggtitle("vitA")
p3 + theme(legend.position = "none")
p3 + guides(color = FALSE)
scorr_est_3 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$estimate)
scorr_pval_3 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$p.value)
##print(p)
##
res1 <- sig_comb.filt %>% dplyr::filter(contrast.x=="vitD")
p4 <- ggplot(res1, aes(x=estimate.y, y=estimate.x, color=contrast.x), show.legend = FALSE) +
    geom_point()+
    xlab("sig_genes") + ylab("sig_peaks")+
    xlim(-4, 10) + ylim(-4, 5)+
    ##stat_poly_line() +
    ##stat_poly_eq() +
    ggtitle("vitD")
p4 + theme(legend.position = "none")
p4 + guides(color = FALSE)
scorr_est_4 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$estimate)
scorr_pval_4 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$p.value)
##print(p)
##
res1 <- sig_comb.filt %>% dplyr::filter(contrast.x=="zinc")
p5 <- ggplot(res1, aes(x=estimate.y, y=estimate.x, color=contrast.x), show.legend = FALSE) +
    geom_point()+
    xlab("sig_genes") + ylab("sig_peaks")+
    xlim(-4, 10) + ylim(-4, 5)+
    ##stat_poly_line() +
    ##stat_poly_eq() +
    ggtitle("zinc")
p5 + theme(legend.position = "none")
p5 + guides(color = FALSE)
scorr_est_5 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$estimate)
scorr_pval_5 <- as.numeric(cor.test(res1$estimate.x, res1$estimate.y, method = "spearman")$p.value)
##print(p)
##
##print(mtext(1, side=4, line=0.5, cex=1, col="blue"))
#}
##dev.off()


print(scorr_est_1)
print(scorr_est_2)
print(scorr_est_3)
print(scorr_est_4)
print(scorr_est_5)
print(scorr_pval_1)
print(scorr_pval_2)
print(scorr_pval_3)
print(scorr_pval_4)
print(scorr_pval_5)


figure <- ggarrange(p1, p2, p3, p4, p5 +
                                    font("x.text", size = 16),
                    ncol = 3, nrow = 2)
png(paste0("./2.2_compareRNAandATAC_3.outs/Figure1.1_significant_scatter_macs2_0.1_cn_peakAnno_control_treat_noline.png"), width=3000, height=2000, pointsize=18, res=175)
print(figure)
dev.off()




##----celltype and treatments----
figfn <- paste0("./2.2_compareRNAandATAC_3.outs/Figure1.1_significant_scatter_macs2_0.1_cn_peakAnno_control_allcomb.png")
png(figfn, width=6000, height=8000, pointsize=33, res=250)
par(mar=c(4,4,2,2),mgp=c(2,1,0))
x <- matrix(1:40, 8, 5, byrow=T)
layout(x)
for (i in comb.3){
    res1 <- sig_comb.filt %>% dplyr::filter(MCls.x==i, contrast.x=="caffeine")
    p <- ggplot(res1, aes(x=estimate.y, y=estimate.x, color=contrast.x)) +
        geom_point()+
        xlab("sig_genes") + ylab("sig_peaks")+
        xlim(-4, 10) + ylim(-4, 5)+
        stat_poly_line() +
        stat_poly_eq() +
        ggtitle("fold.enrichment")
    print(p)
    ##
    res1 <- sig_comb.filt %>% dplyr::filter(MCls.x==i, contrast.x=="nicotine")
    p <- ggplot(res1, aes(x=estimate.y, y=estimate.x, color=contrast.x)) +
        geom_point()+
        xlab("sig_genes") + ylab("sig_peaks")+
        xlim(-4, 10) + ylim(-4, 5)+
        stat_poly_line() +
        stat_poly_eq() +
        ggtitle("fold.enrichment")
    print(p)
    ##
    res1 <- sig_comb.filt %>% dplyr::filter(MCls.x==i, contrast.x=="vitA")
    p <- ggplot(res1, aes(x=estimate.y, y=estimate.x, color=contrast.x)) +
        geom_point()+
        xlab("sig_genes") + ylab("sig_peaks")+
        xlim(-4, 10) + ylim(-4, 5)+
        stat_poly_line() +
        stat_poly_eq() +
        ggtitle("fold.enrichment")
    print(p)
    ##
    res1 <- sig_comb.filt %>% dplyr::filter(MCls.x==i, contrast.x=="vitD")
    p <- ggplot(res1, aes(x=estimate.y, y=estimate.x, color=contrast.x)) +
        geom_point()+
        xlab("sig_genes") + ylab("sig_peaks")+
        xlim(-4, 10) + ylim(-4, 5)+
        stat_poly_line() +
        stat_poly_eq() +
        ggtitle("fold.enrichment")
    print(p)
    ##
    res1 <- sig_comb.filt %>% dplyr::filter(MCls.x==i, contrast.x=="zinc")
    p <- ggplot(res1, aes(x=estimate.y, y=estimate.x, color=contrast.x)) +
        geom_point()+
        xlab("sig_genes") + ylab("sig_peaks")+
        xlim(-4, 10) + ylim(-4, 5)+
        stat_poly_line() +
        stat_poly_eq() +
        ggtitle("fold.enrichment")
    print(p)
    ##
    print(mtext(i, side=4, line=0.5, cex=1, col="blue"))
}
dev.off()
    


end


## ##----plot----
## lm_eqn <- function(df){
##     m <- lm(y ~ x, df);
##     eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
##                      list(a = format(unname(coef(m)[1]), digits = 2),
##                           b = format(unname(coef(m)[2]), digits = 2),
##                           r2 = format(summary(m)$r.squared, digits = 3)))
##     as.character(as.expression(eq));
## }


## png("./2.2_compareRNAandATAC_3.outs/Figure1.1_significant_scatter_macs2_0.1_cn_peakAnno_control_scatter.png")

## library("ggpubr")
## p <- ggscatter(sig_comb.filt, x = "estimate.y", y = "estimate.x",
##           #color=sig_comb.filt$comb.x,
##           add = "reg.line", conf.int = TRUE,
##           cor.coef = TRUE, cor.method = "spearman",
##           xlab = "sig_genes", ylab = "sig_peaks")


## p <- ggplot(sig_comb.filt, aes(x=estimate.y, y=estimate.x, color=comb.x)) +
##     geom_point()+
##     xlab("sig_genes") + ylab("sig_peaks")+
##     xlim(-4, 10) + ylim(-4, 5)+
##     ##stat_poly_line() +
##     ##stat_poly_eq() +
##     ##geom_text(label = lm_eqn(df), parse = TRUE)+
##     ggtitle("fold.enrichment")

## print(p)
## dev.off()

## scorr <- as.numeric(cor.test(sig_comb.filt$estimate.y, sig_comb.filt$estimate.x, method = "spearman")$estimate)
## print(scorr)

## scorr.pval <- as.numeric(cor.test(sig_comb.filt$estimate.y, sig_comb.filt$estimate.x, method = "spearman")$p.value)
## print(scorr.pval)






#################################################
### QQ PLOT 
#################################################
head(resDE)
nrow(resDE)

head(resDP.peakAnno)
nrow(resDP.peakAnno)

head(resDP.pa.sum)
nrow(resDP.pa.sum)

comb <- sort(unique(resDE$comb))
comb

dfNew <- map_dfr(comb, function(ii){
  res2 <- sig_resDP.peakAnno%>%dplyr::filter(comb==ii)
  DEG <- sig_resDE%>%dplyr::filter(comb==ii)%>%dplyr::pull(gene)
  res2 <- res2%>%mutate(is_DEG=ifelse(symbol%in%DEG,1,0))
  ###
  dx <- map_dfr(c(0,1),function(i){
     di <- res2%>%dplyr::filter(is_DEG==i)%>%arrange(p.value)
     ngene <- nrow(di)
     di <- di%>%mutate(observed=-log10(p.value), expected=-log10(ppoints(ngene)))
     di
  })
  dx
})

head(dfNew)
nrow(dfNew)

###
###
## lab1 <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX", "PHA"="PHA", "PHA-DEX"="PHA+DEX")

## p1 <- ggplot(dfNew, aes(x=expected,y=observed, colour=factor(is_DEG)))+
##    ggrastr::rasterise(geom_point(size=0.3),dpi=300)+
##    geom_abline(colour="red")+
##    scale_colour_manual(values=c("0"="grey40","1"="green"),
##       labels=c("0"="Not DEG", "1"="DEG"),
##       guide=guide_legend(override.aes=list(size=1)))+
##    facet_grid(MCls~contrast, scales="free", labeller=labeller(contrast=lab1))+
##    xlab(bquote("Expected"~ -log[10]~"("~italic(plain(P))~")"))+
##    ylab(bquote("Observed"~-log[10]~"("~italic(plain(P))~")"))+
##    theme_bw()+
##    theme(legend.title=element_blank(),strip.text=element_text(size=12))

## figfn <- "./2.2_compareRNAandATAC.outs/Figure1.1_qq.signac.png"
## png(figfn, width=750, height=750, res=120)
## print(p1)
## dev.off()


lab1 <- c("caffeine"="caffeine", "nicotine"="nicotine", "zinc"="zinc",
          "vitA"="vitA", "vitD"="vitD", "vitE"="vitE")
          ##"water"="water")

lab2 <- c("4-Bcell"="4-Bcell", "6-Monocyte"="6-Monocyte",
          "2-NKcell"="2-NKcell", "0-CD4Naive"="0-CD4Naive", "1-TCM"="1-TCM",
          "3-TEM"="3-TEM", "5-CD8Naive"="5-CD8Naive", "7-dnT"="7-dnT")

p <- ggplot(dfNew, aes(x=expected,y=observed, colour=factor(is_DEG)))+
    ggrastr::rasterise(geom_point(size=0.3),dpi=300)+
    geom_abline(colour="red")+
    scale_colour_manual(values=c("1"="green", "0"="grey40"),
                        labels=c("1"="sig DAR, sig DEG", "0"="sig DAR, non-sig DEG"),
                        guide=guide_legend(override.aes=list(size=1)))+
    facet_grid(MCls~contrast, scales="free",
               labeller=labeller(contrast=lab1, MCls=lab2))+
    xlab(bquote("Expected"~ -log[10]~"("~italic(plain(P))~")"))+
    ylab(bquote("Observed"~-log[10]~"("~italic(plain(P))~")"))+
    theme_bw()+
    theme(strip.text=element_text(size=12))

figfn <- "./2.2_compareRNAandATAC_3.outs/Figure1.1_qq.signac_control.png"
png(figfn, width=1600, height=1600, res=175)
print(p)
dev.off()



## #######################################
## ### peak annotation from ChIPseeker ###
## #######################################
## peakAnno <- read_rds("./2.2_compareRNAandATAC.outs/2_annot.ChIPseeker.rds")

## ##
## figfn <- "./2.2_compareRNAandATAC.outs/Figure2.1_annot.pie.png"
## png(figfn, width=480, height=320)
## plotAnnoPie(peakAnno)
## dev.off()

## ###
## figfn <- "./2.2_compareRNAandATAC.outs/Figure2.2_annot.ven.png"
## png(figfn, width=700, height=500)
## vennpie(peakAnno)
## dev.off()

## ###
## peakAnno <- as.data.frame(peakAnno)
## peakAnno <- peakAnno%>%
##    mutate(chr=gsub("chr","",seqnames),
##           peak_region=paste(chr,start,end,sep="-"))%>%
##    dplyr::select(peak_region, geneId, flank_geneIds) 


## ###
## res <- read_rds("./1.2_DiffPeak.outs/2.0_DESeq.results.rds")%>%as.data.frame()
## res <- res%>%dplyr::rename("peak_region"="gene")%>%drop_na(p.value)

## res <- res%>%
##    left_join(peakAnno, by="peak_region")%>%
##    mutate(comb=paste(MCls, contrast, sep="_"))


## ### previous identified DEGs
## fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
## resDE <- read_rds(fn)%>%dplyr::filter(qval<0.1, abs(beta)>0.5)
## resDE <- resDE%>%
##    mutate(comb=paste(MCls, contrast, sep="_"))


## ###
## ### qq for cloest genes
## comb <- sort(unique(resDE$comb))
## dfNew <- map_dfr(comb, function(ii){
##   res2 <- res%>%dplyr::filter(comb==ii)
##   DEG <- resDE%>%dplyr::filter(comb==ii)%>%dplyr::pull(gene)
##   res2 <- res2%>%mutate(is_DEG=ifelse(geneId%in%DEG,1,0))  
##   ###
##   dx <- map_dfr(c(0,1),function(i){
##      di <- res2%>%dplyr::filter(is_DEG==i)%>%arrange(p.value)
##      ngene <- nrow(di)
##      di <- di%>%mutate(observed=-log10(p.value), expected=-log10(ppoints(ngene)))
##      di
##   })
##   dx
## })

## ###
## ###
## lab1 <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX", "PHA"="PHA", "PHA-DEX"="PHA+DEX")
## p1 <- ggplot(dfNew, aes(x=expected,y=observed, colour=factor(is_DEG)))+
##    ggrastr::rasterise(geom_point(size=0.3),dpi=300)+
##    geom_abline(colour="red")+
##    scale_colour_manual(values=c("0"="grey40","1"="green"),
##       labels=c("0"="Not DEG", "1"="DEG"),
##       guide=guide_legend(override.aes=list(size=1)))+
##    facet_grid(MCls~contrast, scales="free", labeller=labeller(contrast=lab1))+
##    xlab(bquote("Expected"~ -log[10]~"("~italic(plain(P))~")"))+
##    ylab(bquote("Observed"~-log[10]~"("~italic(plain(P))~")"))+
##    theme_bw()+
##    theme(legend.title=element_blank(),strip.text=element_text(size=12))

## figfn <- "./2.2_compareRNAandATAC.outs/Figure2.3_qq.DEG.cloest.png"
## png(figfn, width=750, height=750, res=120)
## print(p1)
## dev.off()



## ###
## ### flanking region
## ### generate data for qq plot
## comb <- sort(unique(resDE$comb))
## dfNew <- map_dfr(comb, function(ii){
##   res2 <- res%>%dplyr::filter(comb==ii)
##   DEG <- resDE%>%dplyr::filter(comb==ii)%>%dplyr::pull(gene)
##   geneList <- str_split(res2$flank_geneIds, ";")
##   is_DEG <- map_dbl(geneList, function(i){
##      x <- ifelse(any(i%in%DEG),1,0)
##      x
##   })
##   res2$is_DEG <- is_DEG
  
##   ###
##   dx <- map_dfr(c(0,1),function(i){
##      di <- res2%>%dplyr::filter(is_DEG==i)%>%arrange(p.value)
##      ngene <- nrow(di)
##      di <- di%>%mutate(observed=-log10(p.value), expected=-log10(ppoints(ngene)))
##      di
##   })
##   dx
## })


## ###
## lab1 <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX", "PHA"="PHA", "PHA-DEX"="PHA+DEX")
## p2 <- ggplot(dfNew, aes(x=expected,y=observed, colour=factor(is_DEG)))+
##    ggrastr::rasterise(geom_point(size=0.3),dpi=300)+
##    geom_abline(colour="red")+
##    scale_colour_manual(values=c("0"="grey40","1"="green"),
##       labels=c("0"="Not DEG", "1"="DEG"),
##       guide=guide_legend(override.aes=list(size=1)))+
##    facet_grid(MCls~contrast, scales="free", labeller=labeller(contrast=lab1))+
##    xlab(bquote("Expected"~ -log[10]~"("~italic(plain(P))~")"))+
##    ylab(bquote("Observed"~-log[10]~"("~italic(plain(P))~")"))+
##    theme_bw()+
##    theme(legend.title=element_blank(),strip.text=element_text(size=12))

## figfn <- "./2.2_compareRNAandATAC.outs/Figure2.4_qq.DEG.png"
## png(figfn, width=750, height=750, res=120)
## print(p2)
## dev.off()




## ############################
## ### if enrich are in DVG ###
## ############################
## ### previous identified DEGs
## fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/10_RNA.Variance_output/tmp9/3_phiNew.meta"
## resDE <- read.table(fn,header=T)%>%dplyr::filter(qval<0.1, abs(beta)>0.5)
## resDE <- resDE%>%
##    mutate(comb=paste(MCls, contrast, sep="_"))


## ###
## ### qq for cloest genes
## comb <- sort(unique(resDE$comb))
## dfNew <- map_dfr(comb, function(ii){
##   res2 <- res%>%dplyr::filter(comb==ii)
##   DEG <- resDE%>%dplyr::filter(comb==ii)%>%dplyr::pull(gene)
##   res2 <- res2%>%mutate(is_DEG=ifelse(geneId%in%DEG,1,0))  
##   ###
##   dx <- map_dfr(c(0,1),function(i){
##      di <- res2%>%dplyr::filter(is_DEG==i)%>%arrange(p.value)
##      ngene <- nrow(di)
##      di <- di%>%mutate(observed=-log10(p.value), expected=-log10(ppoints(ngene)))
##      di
##   })
##   dx
## })

## ###
## ###
## lab1 <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX", "PHA"="PHA", "PHA-DEX"="PHA+DEX")
## p1 <- ggplot(dfNew, aes(x=expected,y=observed, colour=factor(is_DEG)))+
##    ggrastr::rasterise(geom_point(size=0.3),dpi=300)+
##    geom_abline(colour="red")+
##    scale_colour_manual(values=c("0"="grey40","1"="green"),
##       labels=c("0"="Not DVG", "1"="DVG"),
##       guide=guide_legend(override.aes=list(size=1)))+
##    facet_grid(MCls~contrast, scales="free", labeller=labeller(contrast=lab1))+
##    xlab(bquote("Expected"~ -log[10]~"("~italic(plain(P))~")"))+
##    ylab(bquote("Observed"~-log[10]~"("~italic(plain(P))~")"))+
##    theme_bw()+
##    theme(legend.title=element_blank(),strip.text=element_text(size=12))

## figfn <- "./2.2_compareRNAandATAC.outs/Figure2.5_qq.DVG.cloest.png"
## png(figfn, width=750, height=750, res=120)
## print(p1)
## dev.off()



###############################################
### forest plots DEG if is enriched in DARs ###
###############################################
### fisher test function
cal.fisher <- function(df){
   ###
 resfisher <- map_dfr(1:nrow(df),function(i){
      dmat <- matrix(as.numeric(df[i,]), 2, 2)
      colnames(dmat) <- c("interest", "not.interest")
      rownames(dmat) <- c("in.DAR", "not.DAR")
      ##rownames(dmat) <- c("in.DEG", "not.DEG")
      res <- fisher.test(dmat)
      res2 <- data.frame(odds=res$estimate, pval.fisher=res$p.value,
         lower=res$conf.int[1], upper=res$conf.int[2])
      res2
   })
   resfisher
}

## ##
## df.fisher <- function(resDE, resDP, peakAnno){
##    ### 
##    comb <- sort(unique(resDE$comb))
##    df.ls <- lapply(1:length(comb), function(i){
##       ii <- comb[i] 
##       cell <- gsub("_.*", "", ii)
##       contrast <- gsub(".*_", "", ii)
##       resDE2 <- resDE%>%dplyr::filter(comb==ii) ## dplyr::filter(qval<0.1, abs(beta)>0.5)
##       resDP2 <- resDP%>%
##          dplyr::filter(comb==ii, p.adjusted<0.1, abs(estimate)>0.5)%>%
##          dplyr::left_join(peakAnno, by="peak")
##       if( nrow(resDP2)>5){
##          sigGene <- resDE2%>%
##             dplyr::filter(qval<0.1, abs(beta)>0.5)%>%
##             dplyr::pull(gene)
##          ###
##          interest.in.DARs <- sum(sigGene%in%resDP2$geneID)
##          interest.not.DARs <- length(sigGene)-interest.in.DARs
##       ###
##          notSig <- setdiff(res2$gene, sigGene)
##          not.interest.in.DARs <- sum(notSig%in%resDP2$geneId)
##          not.interest.not.DARs <- length(notSig)-not.interest.in.DARs
##          df2 <- data.frame("interest.in.DARs"=interest.in.DARs,
##             "interest.not.DARs"=interest.not.DARs,
##             "not.interest.in.DARs"=not.interest.in.DARs,
##             "not.interest.not.DARs"=not.interest.not.DARs)
##          df2$cell <- cell
##          df2$contrast <- contrast
##          df2$comb <- ii
##       }else{
##         df2 <- NA    
##       }
##       df2
##    })    
##    df.ls <- df.ls[!is.na(df.ls)]
##    df <- do.call(rbind, df.ls)    
##    res <- cal.fisher(df[,1:4])
##    res2 <- cbind(df, res)
##    ### 
##    as.data.frame(res2)
## }    

head(resDP.peakAnno)
nrow(resDP.peakAnno)
colnames(resDP.peakAnno)

df.fisher <- function(resDE, resDP.peakAnno){
   ### 
   comb <- sort(unique(resDE$comb))
   df.ls <- lapply(1:length(comb), function(i){
      ii <- comb[i] 
      contrast <- gsub("_.*", "", ii)
      cell <- gsub(".*_", "", ii)
      ##
      resDE2 <- resDE%>%dplyr::filter(comb==ii) ## dplyr::filter(qval<0.1, abs(beta)>0.5)
      ##resDP2 <- resDP%>%
      ##   dplyr::filter(comb==ii, p.adjusted<0.1, abs(estimate)>0.5)%>%
      ##   dplyr::left_join(peakAnno, by="peak")
      resDP2 <- resDP.peakAnno%>%dplyr::filter(comb==ii)
      ##
      ##
      if(nrow(resDP2)>5){
          ##
          ##
          sigGene <- resDE2%>%dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)%>%
              dplyr::pull(gene)
          sigPeak <- unique(resDP2%>%dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)%>%dplyr::pull(symbol))
          ##nonsigPeak <- unique(resDP2%>%dplyr::filter(p.adjusted>0.1 | (p.adjusted<0.1 & abs(estimate)>0.5))%>%dplyr::pull(symbol))
          notSig <- setdiff(resDE2$gene, sigGene)
          ##
          ##
          interest.in.DARs <- sum(sigGene%in%sigPeak)  #%in%resDP2$symbol)
          interest.not.DARs <- length(sigGene)-interest.in.DARs
          ##interest.in.DEGs <- sum(sigPeak%in%sigGene)
          ##interest.not.DEGs <- sum(sigPeak%in%notSig)
          ##interest.in.DEGs <- sum(sigPeak%in%sigGene)
          ##interest.not.DEGs <- length(sigPeak)-interest.in.DEGs
          ##
          ##
          not.interest.in.DARs <- sum(notSig%in%sigPeak) #%in%resDP2$symbol)
          not.interest.not.DARs <- length(notSig)-not.interest.in.DARs
          ##not.interest.in.DEGs <- sum(nonsigPeak%in%sigGene)
          ##not.interest.not.DEGs <- sum(nonsigPeak%in%notSig)
          ##not.interest.in.DEGs <- sum(nonsigPeak%in%sigGene)
          ##not.interest.not.DEGs <- length(nonsigPeak)-not.interest.in.DEGs
          ##
          ##
          df2 <- data.frame("interest.in.DARs"=interest.in.DARs,
                            "interest.not.DARs"=interest.not.DARs,
                            "not.interest.in.DARs"=not.interest.in.DARs,
                            "not.interest.not.DARs"=not.interest.not.DARs)
          #df2 <- data.frame("interest.in.DEGs"=interest.in.DEGs,
          #                  "interest.not.DEGs"=interest.not.DEGs,
          #                  "not.interest.in.DEGs"=not.interest.in.DEGs,
          #                  "not.interest.not.DEGs"=not.interest.not.DEGs)
          ##
          ##
          df2$cell <- cell
          df2$contrast <- contrast
          df2$comb <- ii
      }else{
          df2 <- NA    
      }
      df2
   })    
    df.ls <- df.ls[!is.na(df.ls)]
    df <- do.call(rbind, df.ls)    
    res <- cal.fisher(df[,1:4])
    res2 <- cbind(df, res)
    ##
    as.data.frame(res2)
}
    

resFisher <- df.fisher(resDE, resDP.peakAnno)
head(resFisher)

#opfn <- "./2.2_compareRNAandATAC_3.outs/3.1_enrich.DEG_control.csv"
#write.csv(resFisher, opfn, row.names=F)

opfn <- "./2.2_compareRNAandATAC_3.outs/3.1_enrich.DEG_howmanygenes_control.csv"
write.csv(resFisher, opfn, row.names=F)



#################
### read data ###
#################
## fn <- "./1.2_DiffPeak.outs/2.0_DESeq.results_control.rds"
## resDP <- read_rds(fn)%>%drop_na(p.value)%>%
##    mutate(comb=paste(MCls, contrast, sep="_"))%>%
##    dplyr::rename("peak"="gene") 


## fn <- "./2.2_compareRNAandATAC.outs/2_annot.ChIPseeker.rds"
## peakAnno <- read_rds(fn)%>%as.data.frame()%>%
##    mutate(peak=paste(gsub("chr", "", seqnames), start, end, sep="-")) 
## peakAnno2 <- peakAnno%>%dplyr::select(peak, geneId, SYMBOL)


## ### test if DEG is DARs 
## fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
## res <- read_rds(fn)%>%drop_na(pval)%>%
##    mutate(comb=paste(MCls, contrast, sep="_"))%>%
##    as.data.frame()
## ##
## resFisher <- df.fisher(res, resDP, peakAnno2)
## opfn <- "./2.2_compareRNAandATAC.outs/3.1_enrich.DEG.csv"
## write.csv(resFisher, opfn, row.names=F)


## ### test if DVGs is DARs
## fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/10_RNA.Variance_output/tmp9/3_phiNew.meta"
## res <- read.table(fn, header=T) %>%drop_na(pval)%>%
##    mutate(comb=paste(MCls, contrast, sep="_"))

## resFisher2 <- df.fisher(res, resDP, peakAnno2)
## opfn <- "./2.2_compareRNAandATAC.outs/3.2_enrich.DVG.csv"
## write.csv(resFisher2, opfn, row.names=F)







#######################
### error bar plots ###
#######################
## col.treat <- c("LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
##              "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
## col.MCl <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
##              "NKcell"="#aa4b56", "Tcell"="#ffaa00")




col.treat <- c("caffeine"="red", "nicotine"="tan",
          "vitA"="tan4", "vitD"="seagreen4",
          "vitE"="salmon3",
          ##"water"="grey",
          "zinc"="maroon3")
col.MCl <- c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
          "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00", "1-TCM"="pink",
          "3-TEM"="blue", "5-CD8Naive"="green", "7-dnT"="black")



## lab1 <- c("caffeine"="caffeine", "nicotine"="nicotine", "zinc"="zinc",
##           "vitA"="vitA", "vitD"="vitD", "vitE"="vitE", "water"="water")

## lab2 <- c("4_Bcell"="4_Bcell", "6_Monocyte"="6_Monocyte",
##           "2_NKcell"="2_NKcell", "0_CD4Naive"="0_CD4Naive", "1_TCM"="1_TCM",
##           "3_TEM"="3_TEM", "5_CD8Naive"="5_CD8Naive", "7_dnT"="7_dnT")




fn <- "./2.2_compareRNAandATAC_3.outs/3.1_enrich.DEG_control.csv"
df <- read.csv(fn, header=T)
head(df)
nrow(df)

df <- resFisher
head(df)

## no need to calculate qval as only one row for each contrast-cell combination
#df <- df%>%group_by(cell, contrast)%>%
#       mutate(qvalue.fisher=p.adjust(pval.fisher, "BH"))%>%
#       ungroup()%>%
#       as.data.frame()
#head(df)
#nrow(df)

#df <- df %>% dplyr::filter(pval.fisher<0.1)
#nrow(df)


df%>%dplyr::filter(contrast=="caffeine")
df%>%dplyr::filter(contrast=="zinc")


df1 <- data.frame(odds=df$odds,
   CI.low=df$lower, CI.high=df$upper,
   comb=gsub("_", ".", df$comb),
   MCls=df$cell, contrast=df$contrast, gene="DEG")


## df1.notinf <- df1 %>% dplyr::filter(CI.high!=Inf)
## head(df1.notinf)
## error.max <- max(df1.notinf$CI.high)
## error.max
## high.quant <- quantile(abs(df1.notinf$CI.high),probs=0.75)
## high.quant

head(df1)

## fn <- "./2.2_compareRNAandATAC.outs/3.2_enrich.DVG.csv"
## df <- read.csv(fn, header=T)
## df2 <- data.frame(odds=df$odds,
##    CI.low=df$lower, CI.high=df$upper,
##    comb=gsub("_", ".", df$comb),
##    MCls=df$cell, contrast=df$contrast, gene="DVG")

## df2 <- rbind(df1, df2)

df2 <- df1

comb2 <- sort(unique(df2$comb))
comb2

##ylab2 <- gsub("-", "+", gsub(".*\\.", "", df2$comb))
ylab2 <- gsub(".*\\.", "", df2$comb)
names(ylab2) <- comb2
ylab2

p <- ggplot(df2, aes(x=odds, y=comb))+
   ##geom_errorbarh(aes(xmax=CI.high, xmin=CI.low, colour=MCls),
   ##    size=0.5, height=0.2)+
   geom_point(aes(colour=contrast), shape=19, size=1.5)+
   scale_colour_manual(values=col.treat)+
   geom_vline(aes(xintercept=1), size=0.25, linetype="dashed")+ 
   xlab("Odds ratio")+
   scale_y_discrete(labels=ylab2)+
   facet_wrap(~factor(gene), ncol=2)+ 
   theme_bw()+
   theme(plot.title=element_text(hjust=0.5, size=9),
         axis.title.y=element_blank(),
         axis.title.x=element_text(size=9),
         ## axis.text=element_text(size=8),
         legend.position="right")

###
figfn <- "./2.2_compareRNAandATAC_3.outs/Figure3.0_enrich.DAR_howmanygenes_control.png"
png(figfn, width=1500, height=1500, res=175)
p
dev.off()



##
### poster
p <- ggplot(df2%>%dplyr::filter(gene=="DEG"), aes(y=odds, x=comb))+
   geom_errorbar(aes(ymax=CI.high, ymin=CI.low, colour=MCls),
       size=0.5, width=0.2)+
   geom_point(aes(colour=MCls), shape=19, size=1.5)+
   scale_colour_manual(values=col.MCl)+
   geom_hline(aes(yintercept=1), size=0.25, linetype="dashed")+ 
   ylab("Odds ratio")+
   scale_x_discrete(labels=ylab2)+
   ## facet_wrap(~factor(gene), ncol=2)+ 
   theme_bw()+
   theme(plot.title=element_text(hjust=0.5, size=9),
         axis.title.x=element_blank(),
         axis.title.y=element_text(size=12),
         axis.text.x=element_text(size=10, angle=-90, hjust=0, vjust=0.5),
         axis.text.y=element_text(size=10),
         legend.position="none")

###
figfn <- "./2.2_compareRNAandATAC.outs/Figure3.1_enrich.DAR.png"
png(figfn, width=620, height=480, res=120)
p
dev.off()


###
p2 <- ggplot(df2%>%filter(gene=="DEG"), aes(x=odds, y=comb))+
   geom_errorbarh(aes(xmax=CI.high, xmin=CI.low, colour=MCls),
       size=0.5, height=0.2)+
   geom_point(aes(colour=MCls), shape=19, size=1.5)+
   scale_colour_manual(values=col.MCl)+
   geom_vline(aes(xintercept=1), size=0.25, linetype="dashed")+ 
   xlab("Odds ratio")+
   scale_y_discrete(labels=ylab2)+
   ## facet_wrap(~factor(gene), ncol=2)+ 
   theme_bw()+
   theme(plot.title=element_text(hjust=0.5, size=9),
         axis.title.y=element_blank(),
         axis.title.x=element_text(size=9),
         ## axis.text=element_text(size=8),
         legend.position="none")
###
figfn <- "./2.2_compareRNAandATAC.outs/Figure3.4_enrich.DAR.png"
png(figfn, width=380, height=480, res=120)
p2
dev.off()

##                  
## p1 <- ggplot(df2, aes(x=odds, y=comb))+
##    geom_errorbarh(aes(xmax=CI.high, xmin=CI.low, colour=MCls),
##        size=0.5, height=0.2)+
##    geom_point(aes(colour=MCls), shape=19, size=1.5)+
##    scale_colour_manual(values=col.MCl)+
##    geom_vline(aes(xintercept=1), size=0.25, linetype="dashed")+ 
##    xlab("Odds ratio")+
##    ## scale_y_discrete(labels=ylab2)+ 
##    ggtitle("DEGs")+    
##    theme_bw()+
##    theme(plot.title=element_text(hjust=0.5, size=9),
##          axis.title.y=element_blank(),
##          axis.title.x=element_text(size=9),
##          axis.text=element_text(size=8),
##          legend.position="none")

## ###
## fn <- "./2.2_compareRNAandATAC.outs/3.2_enrich.DVG.csv"
## df <- read.csv(fn, header=T)
## df2 <- data.frame(odds=df$odds,
##    CI.low=df$lower, CI.high=df$upper,
##    comb=gsub("_", ".", df$comb),
##    MCls=df$cell, contrast=df$contrast)

## p2 <- ggplot(df2, aes(x=odds, y=comb))+
##    geom_errorbarh(aes(xmax=CI.high, xmin=CI.low, colour=MCls),
##        size=0.5, height=0.2)+
##    geom_point(aes(colour=MCls), shape=19, size=1.5)+
##    scale_colour_manual(values=col.MCl)+
##    geom_vline(aes(xintercept=1), size=0.25, linetype="dashed")+ 
##    xlab("Odds ratio")+
##    ## scale_y_discrete(labels=ylab2)+   
##    ggtitle("DVGs")+    
##    theme_bw()+
##    theme(plot.title=element_text(hjust=0.5, size=9),
##          axis.title.y=element_blank(),
##          axis.title.x=element_text(size=9),
##          axis.text=element_text(size=8),
##          legend.position="none")

## ###
## figfn <- "./2.2_compareRNAandATAC.outs/Figure3_enrich.DAR.png"
## png(figfn, width=680, height=480, res=120)
## plot_grid(p1, p2, ncol=2)
## dev.off()



#######################
### Distance to TSS ###
#######################

###
### DEGs or DVGs are Distance to peaks

###
### differential peaks
fn <- "./1.2_DiffPeak.outs/2.0_DESeq.results.rds"
resDP <- read_rds(fn)%>%drop_na(p.value)%>%
   mutate(comb=paste(MCls, contrast, sep="_"))%>%
   dplyr::rename("peak"="gene")%>%
   dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)%>%
   as.data.frame()

peakAll <- unique(resDP$peak)


###
### annotation
fn <- "./2.2_compareRNAandATAC.outs/2_annot.ChIPseeker.rds"
peakAnno <- read_rds(fn)%>%as.data.frame()%>%
   mutate(peak=paste(gsub("chr", "", seqnames), start, end, sep="-")) 
peakAnno2 <- peakAnno%>%dplyr::filter(peak%in%peakAll)%>%
   dplyr::select(peak, geneId, SYMBOL, dtss=distanceToTSS)%>%mutate(dtss=abs(dtss))


 ## df2gene <- peakAnno2%>%
 ##     group_by(geneId)%>%summarise(npeaks=n(), .groups="drop")%>%
 ##     ungroup()%>%as.data.frame()


###
### DEG is DARs 
fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
resDEG <- read_rds(fn)%>%drop_na(pval)%>%
   mutate(comb=paste(MCls, contrast, sep="_"))%>%
   dplyr::filter(qval<0.1, abs(beta)>0.5)%>% 
   as.data.frame()

## DEG <- res%>%dplyr::filter(qval<0.1, abs(beta)>0.5)%>%
##     dplyr::pull(gene)%>%unique()
##

###
### DVG
fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/10_RNA.Variance_output/tmp9/3_phiNew.meta"
resDVG <- read.table(fn, header=T) %>%drop_na(pval)%>%
   mutate(comb=paste(MCls, contrast, sep="_"))%>%
   dplyr::filter(qval<0.1, abs(beta)>0.5)%>%
   as.data.frame()


###
###
comb <- sort(unique(resDEG$comb))
geneList <- lapply(comb, function(ii){
###
   resDP2 <- resDP%>%dplyr::filter(comb==ii)
   resDE2 <- resDEG%>%dplyr::filter(comb==ii)
   resDV2 <- resDVG%>%dplyr::filter(comb==ii)
###
   DP <- unique(resDP2$peak)
   DEG <- unique(resDE2$gene)
   DVG <- unique(resDV2$gene)
    
   DEGunq <- setdiff(DEG, DVG)
   shared <- intersect(DEG, DVG)
   DVGunq <- setdiff(DVG, DEG)
   DGall <- union(DEG, DVG)

   ## DEG 
   x <- peakAnno2%>%dplyr::filter(peak%in%DP, geneId%in%DEG)
   if(nrow(x)==0){
      res1 <- NA
    }else{
       res1 <- unique(x$geneId)
    }
    ##
    x <- peakAnno2%>%dplyr::filter(peak%in%DP, geneId%in%DVG)
    if(nrow(x)==0){
      res2 <- NA
    }else{
       res2 <- unique(x$geneId)
    }
    res <- list(DEG=res1, DVG=res2)
})

gene1 <- lapply(geneList, function(x) x[[1]])
gene1 <- unique(unlist(gene1[!is.na(gene1)]))

gene2 <- lapply(geneList, function(x) x[[2]])
gene2 <- unique(unlist(gene2[!is.na(gene2)]))

DEG <- unique(resDEG$gene)
DVG <- unique(resDVG$gene)
## DEG, 3115/6571, 47%
## DVG, 579/1409, 41%




###
### compute the number of differential peaks nearby DEG or DVG genes and distance to DPs of DEG or DVG 
comb <- sort(unique(resDEG$comb))
dd <- map_dfr(comb, function(ii){
###
   resDP2 <- resDP%>%dplyr::filter(comb==ii)
   resDE2 <- resDEG%>%dplyr::filter(comb==ii)
   resDV2 <- resDVG%>%dplyr::filter(comb==ii)
###
   DP <- unique(resDP2$peak)
   DEG <- unique(resDE2$gene)
   DVG <- unique(resDV2$gene)
    
   DEGunq <- setdiff(DEG, DVG)
   shared <- intersect(DEG, DVG)
   DVGunq <- setdiff(DVG, DEG)
   DGall <- union(DEG, DVG)
   
###
  
   ## d0 <- peakAnno2%>%dplyr::filter(peak%in%DP, !geneId%in%DGall)
   ## d0 <- d0%>%group_by(geneId)%>%
   ##    summarise(npeaks=n(), dtss=median(dtss), .groups="drop")%>%
   ##    ungroup()%>%as.data.frame()%>%
   ##    mutate(conditions=ii, grp="0")
    
  ###  
  d1 <- peakAnno2%>%dplyr::filter(peak%in%DP, geneId%in%DEG)  
  d1 <- d1%>%group_by(geneId)%>%
     summarise(npeaks=n(), dtss=min(dtss), .groups="drop")%>%
     ungroup()%>%as.data.frame()%>%
     mutate(conditions=ii, grp="1")
  ##
  ## d2 <- peakAnno2%>%dplyr::filter(peak%in%DP, geneId%in%shared)  
  ## d2 <- d2%>%group_by(geneId)%>%
  ##    summarise(npeaks=n(), dtss.min=min(dtss),.groups="drop")%>%
  ##    ungroup()%>%as.data.frame()%>%
  ##    mutate(conditions=ii, grp="2")
  ##
  d2 <- peakAnno2%>%dplyr::filter(peak%in%DP, geneId%in%DVG)  
  d2 <- d2%>%group_by(geneId)%>%
     summarise(npeaks=n(), dtss=min(dtss), .groups="drop")%>%
     ungroup()%>%as.data.frame()%>%
     mutate(conditions=ii, grp="2")
  ##
  dd <- rbind(d1, d2)  
#
  dd    
})    


          
###
p1 <- ggplot(dd, aes(x=factor(grp), y=log2(npeaks), fill=factor(grp)))+
   geom_violin()+
   scale_x_discrete(labels=c("0"="background", "1"="DEG", "3"="DVG-only"))+ 
   ## scale_fill_manual(values=c("DEG"="#d95f02", "DVG"="#1b9e77"))+
   ## scale_color_manual(values=c("DEG"="#d95f02", "DVG"="#1b9e77"))+
   ylab(bquote("Complexity ("~log[2]~" peaks per gene)"))+ 
   theme_bw()+
   theme(legend.position="none",
         axis.title.x=element_blank(),
         legend.background=element_blank()) 

figfn <- "./2.2_compareRNAandATAC.outs/Figure4.1_complexity_violin.png"
png(figfn, width=480, height=320, res=120)
print(p1)
dev.off()

          
###
p2 <- ggplot(dd, aes(x=factor(grp), y=log2(dtss+0.1), fill=factor(grp)))+
   geom_violin()+
   scale_x_discrete(labels=c("0"="background", "1"="DEG", "3"="DVG-only"))+ 
   ## scale_fill_manual(values=c("DEG"="#d95f02", "DVG"="#1b9e77"))+
   ## scale_color_manual(values=c("DEG"="#d95f02", "DVG"="#1b9e77"))+
   ylab(bquote("Distance to neary peaks"~log[2]))+ 
   theme_bw()+
   theme(legend.position="none",
         axis.title.x=element_blank(),
         legend.background=element_blank()) 

figfn <- "./2.2_compareRNAandATAC.outs/Figure4.2_dtss_violin.png"
png(figfn, width=480, height=320, res=120)
print(p2)
dev.off()


###
### summary number of peaks
dd1 <- dd%>%dplyr::filter(grp==1)
sum(dd1$npeaks<=10)/nrow(dd1)
sum(abs(dd1$dtss)<3000)/nrow(dd1)
## dd1 <- dd1%>%group_by(npeaks)%>%summarise(freq=n(),.groups="drop")%>%
##    mutate(prop=freq/sum(freq))%>%
##    ungroup()%>%as.data.frame()

## x1 <- data.frame(order=1, freq=sum(dd1[dd1$npeaks>=1&dd1$npeaks<8,2]), prop=sum(dd1[dd1$npeaks>=1&dd1$npeaks<8,3]))
## ## x2 <- data.frame(order=2, freq=sum(dd1[dd1$npeaks>=5&dd1$npeaks<11,2]), prop=sum(dd1[dd1$npeaks>=5&dd1$npeaks<11,3]))
## x2 <- data.frame(order=2, freq=sum(dd1[dd1$npeaks>=8,2]), prop=sum(dd1[dd1$npeaks>=8,3]))                            
## a <- rbind(x1,x2)
## a$grp <- 1

dd2 <- dd%>%dplyr::filter(grp==2)
sum(dd2$npeaks<=10)/nrow(dd2)
sum(abs(dd2$dtss)<3000)/nrow(dd2)


###
