###
library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(SeuratObject)
library(Signac)
library(SeuratWrappers)
library(SeuratData)
library(qqman)
library(JASPAR2020)
library(JASPAR2022)
library(TFBSTools)
## library(JASPAR2022, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/") 
## library(TFBSTools, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
library(BSgenome.Hsapiens.1000genomes.hs37d5, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(ChIPseeker, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
###
library(ggplot2)
library(cowplot) ##lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(grid)
library(gridExtra)
library(ggExtra)
library(RColorBrewer)
library(ggsci)
library(viridis)
library(ComplexHeatmap)
library(circlize)
theme_set(theme_grey())

outdir <- "./1_motif.outs/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


###########################
### obtain motif object ###
###########################
pfm <- getMatrixSet(x=JASPAR2022,
   opts=list(species=9606, all_version=F))

pfm <- getMatrixSet(x = JASPAR2022, opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))


##
atac <- read_rds("../1_processing/5.1_reCallPeak.outs/3_scATAC.annot.mitofilt.ndup.rds")
atac
atac[["ATAC"]]

atac.2 <- keepStandardChromosomes(atac, pruning.mode = "coarse")
atac.2
atac.2[["ATAC"]]

x <- atac@assays$ATAC@ranges
x
x@seqinfo@seqnames
autosome <- paste("chr", as.character(1:22), sep="")

## anno <- GetGRangesFromEnsDb(ensdb=EnsDb.Hsapiens.v75)
## seqlevelsStyle(anno) <- "UCSC"
## genome(anno) <- "hg19"
## Annotation(atac) <- anno

## atac <- AddMotifs(object=atac,
##    genome=BSgenome.Hsapiens.1000genomes.hs37d5,
##    pfm=pfm)




#----ANNOTATION----
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"
Annotation(atac) <- annotations


#----ADDMOTIFS----
atac <- AddMotifs(object=atac,
   genome=BSgenome.Hsapiens.UCSC.hg38,
   pfm=pfm)


#----SAVE FILE----
opfn <- "./1_motif.outs/1_scATAC.motif.rds"
write_rds(atac, opfn)








########################################
### MOTIF EXAMPLE                    ###
########################################
#----MOTIF FILE SAVED IN THIS SCRIPT----
fn <- "./1_motif.outs/1_scATAC.motif.rds" 
atac <- read_rds(fn)

motif <- Motifs(atac)
motif
#head(motif)

motif.col <- colnames(motif)
motif.col
length(motif.col)

#----count----
motif.matrix <- GetMotifData(object = atac)
motif.matrix
#motif_data <- SetMotifData(object = atac, assay = 'ATAC', slot = 'data', new.data = motif.matrix)
#motif_data

motif.matrix.df <- as.data.frame(motif.matrix)
head(motif.matrix.df)

motif.matrix.peak <- rownames(motif.matrix) 
motif.matrix.peak

cvt <- str_split(rownames(motif.matrix), "-", simplify=T)
table(cvt[,1])

## x <- atac@assays$ATAC@ranges
## #x

## xrange <- ranges(x)
## #xrange

## count <- motif.matrix
## #head(count)

## anno <- data.frame(rn=rownames(count), rnz=rowSums(count),
##                                       chr=as.character(seqnames(x)), start=start(xrange), end=end(xrange))
## head(anno)

## autosome <- paste("chr", as.character(1:22), sep="")
## autosome

## annoSel <- anno%>%dplyr::filter(rnz>0, chr%in%autosome)
## head(annoSel)

## Y <- count[annoSel$rn,]
## head(Y)

## cvt <- str_split(rownames(Y), "-", simplify=T)
## table(cvt[,1])









#----example----
## example <- c("MA0137.3", "MA0517.1", "MA0144.2",
##    "MA0105.4", "MA0778.1",
##    "MA0099.3", "MA1126.1", "MA1134.1", "MA1141.1")

example <- c("MA0137.3", "MA0517.1", "MA0144.2",
             "MA0105.4", "MA0778.1", "MA0099.3",
             "MA1563.2", "MA1134.1", "MA1141.1")


example %in% motif.col

#motif.matrix <- as.matrix(motif)

fig <- MotifPlot(object=atac,motifs=example, facet="wrap",ncol=3,nrow=3)

figfn <- "./1_motif.outs/Figure1.2_motif.example.png"
png(figfn, width=800, height=600,res=120)
print(fig)
dev.off()







############################################################################
############################################################################
### enrichment analysis for motif ###
############################################################################
############################################################################
rm(list=ls())

#----MOTIF FILE SAVED IN THIS SCRIPT----
fn <- "./1_motif.outs/1_scATAC.motif.rds" 
atac <- read_rds(fn)

atac
atac[["ATAC"]]

motif <- Motifs(atac)
motif

#Motifs(atac[["ATAC"]]) <- motif
#atac


head(atac@meta.data)

table(atac@meta.data$MCls)

## motif <- Motifs(atac)
## pfm <- GetMotifData(object=motif, slot="pwm")
## motif <- SetMotifData(object=motif, slot="pwm", new.data=pfm)
## differential peaks


#----RESDP FILE----
#fn <- "../2_Differential/1.2_DiffPeak.outs/2.0_DESeq.results.rds"
#fn <- "../2_Differential/1.2_DiffPeak.outs/2.0_DESeq.results_0.1_cn.rds"
#fn <- "../2_Differential/1_DiffPeak.outs/2.0_DESeq.results_clusters_treatandind_filtered_0.1_cn.rds"
fn <- "../2_Differential/1.2_DiffPeak.outs/2.0_DESeq.results_0.1_cn_sampleID+new_treat.rds"
resDP <- read_rds(fn)%>%as.data.frame()
head(resDP)
nrow(resDP)
table(resDP$MCls) 

#----fisher test----
cal.fisher <- function(df){
  ###  
  resfisher <- map_dfr(1:nrow(df),function(i){  
     dmat <- matrix(as.numeric(df[i,]),2,2)
     colnames(dmat) <- c("interest", "not.interest")
     rownames(dmat) <- c("in.motif", "not.motif")
     res <- fisher.test(dmat, alternative="greater")
     res2 <- data.frame(odds=res$estimate, pval.fisher=res$p.value)
     res2
  })
  resfisher
}


#----data----
## MCls <- rep(c("Bcell", "Monocyte", "NKcell", "Tcell"), each=4)
## contrast <- rep(c("LPS", "LPS-DEX", "PHA", "PHA-DEX"), times=4)
## dataset <- data.frame(MCls=MCls, contrast=contrast)

MCls <- rep(unique(resDP$MCls), each=length(unique(resDP$contrast)))
MCls

contrast <- rep(unique(resDP$contrast), times=length(unique(resDP$MCls)))
contrast

total.comb <- length(unique(resDP$contrast))*length(unique(resDP$MCls))
total.comb

dataset <- data.frame(MCls=MCls, contrast=contrast)
dataset
dataset[1,1]




#----enrichment motif analysis-individual step----
## #enriched.motif <- lapply(1:total.comb, function(i){
## ###
## cell0 <- dataset[1,1]
## cell0

## contrast0 <- dataset[1,2]

## cat(1, cell0, contrast0, "\n") 


## ###
## res2 <- resDP%>%dplyr::filter(MCls==cell0, contrast==contrast0)

## top.DP <- res2%>%dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)%>%dplyr::pull(gene)%>%as.character()
## top.DP

## bg.DP <- res2%>%dplyr::pull(gene)%>%as.character() 
## bg.DP

## n.interest <- length(top.DP)
## n.interest

## n.not <- length(bg.DP)-n.interest
## n.not

## ###    
## #if(n.interest>0){

## enrich2 <- FindMotifs(object=atac, features=top.DP, background=bg.DP)

## head(enrich2)
## nrow(enrich2)

## df <- data.frame("interest.in.motif"=enrich2$observed,
##                      "interest.not.motif"=n.interest-enrich2$observed,
##                      "not.interest.in.motif"=enrich2$background-enrich2$observed)

## head(df)

## df <- df%>%mutate("not.interest.not.motif"=n.not-not.interest.in.motif)

## head(df)

## ## cal.fisher <- function(df){
## ## ###  
## ## resfisher <- map_dfr(1:nrow(df),function(i){  

## ## dmat <- matrix(as.numeric(df[1,]),2,2)
## ## #dmat
## ## colnames(dmat) <- c("interest", "not.interest")
## ## rownames(dmat) <- c("in.motif", "not.motif")
## ## #dmat

## ## res <- fisher.test(dmat, alternative="greater")
## ## #res

## ## res2 <- data.frame(odds=res$estimate, pval.fisher=res$p.value)
## ## res2
## ## })
## ## resfisher
## ## }

## fisher <- cal.fisher(df)
## head(fisher)

## enrich2 <- cbind(enrich2,fisher)
## head(enrich2)

## enrich2$qvalue.hyper <- p.adjust(enrich2$pvalue,"BH")
## enrich2$qvalue.fisher <- p.adjust(enrich2$pval.fisher,"BH")
## enrich2$MCls <- cell0
## enrich2$contrast <- contrast0
## head(enrich2)

## #   }else{
## #     enrich2 <- NA   
## #   }
## #   enrich2 
## #})

## enriched.motif <- enriched.motif[!is.na(enriched.motif)]
## head(enriched.motif)

## enriched.motif <- do.call(rbind,enriched.motif)
## head(enriched.motif)


#----enrichment motif analysis----
enriched.motif <- lapply(1:total.comb, function(i){
###
   cell0 <- dataset[i,1]
   contrast0 <- dataset[i,2]
   cat(i, cell0, contrast0, "\n") 
###
   res2 <- resDP%>%dplyr::filter(MCls==cell0, contrast==contrast0)
   top.DP <- res2%>%
       dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)%>%
       dplyr::pull(gene)%>%as.character()
   bg.DP <- res2%>%dplyr::pull(gene)%>%as.character() 
   n.interest <- length(top.DP)
   n.not <- length(bg.DP)-n.interest
###    
   if(n.interest>0){
     enrich2 <- FindMotifs(
        object=atac,
        features=top.DP,
        background=bg.DP)
     df <- data.frame("interest.in.motif"=enrich2$observed,
                     "interest.not.motif"=n.interest-enrich2$observed,
                     "not.interest.in.motif"=enrich2$background-enrich2$observed)
     df <- df%>%mutate("not.interest.not.motif"=n.not-not.interest.in.motif)
     fisher <- cal.fisher(df)
     enrich2 <- cbind(enrich2,fisher)
     enrich2$qvalue.hyper <- p.adjust(enrich2$pvalue,"BH")
     enrich2$qvalue.fisher <- p.adjust(enrich2$pval.fisher,"BH")
     enrich2$MCls <- cell0
     enrich2$contrast <- contrast0
   }else{
     enrich2 <- NA   
   }
   enrich2 
})

enriched.motif <- enriched.motif[!is.na(enriched.motif)]
enriched.motif <- do.call(rbind,enriched.motif)

opfn <- "./1_motif.outs/2_motif.enrich_macs2_0.1_cn_control.rds"
write_rds(enriched.motif, opfn)

enriched.motif <- read_rds("./1_motif.outs/2_motif.enrich_macs2_0.1_cn_control.rds")

#enriched.motif[1,]

head(enriched.motif)

dim(enriched.motif)
colnames(enriched.motif)

dim(enriched.motif%>%dplyr::filter(MCls=="0_CD4Naive", contrast=="caffeine"))

max(enriched.motif$qvalue.fisher)
min(enriched.motif$qvalue.fisher)

max(enriched.motif$fold.enrichment)
min(enriched.motif$fold.enrichment)

nrow(enriched.motif)
nrow(enriched.motif%>%dplyr::filter(qvalue.fisher<=0.1))








#-----------------MOHAMMED HUSAIN----------------------------------
######################################################
# significant motifs
######################################################
#----enriched.motif----
enriched.motif <- read_rds("./1_motif.outs/2_motif.enrich_macs2_0.1_cn_control.rds")

res <- enriched.motif %>% as.data.frame()
head(res)
nrow(res)

max(res$fold.enrichment)
min(res$fold.enrichment)

res <- res %>% drop_na(pvalue) %>% mutate(MCls=gsub("_", "-", MCls)) %>% mutate(comb=paste(MCls, contrast, sep="_"))

res <- res %>% mutate(logFC=log2(fold.enrichment)) 
max(res$logFC)
min(res$logFC)

## nrow(res %>% dplyr::filter(fold.enrichment==0))
## res0 <- res %>% dplyr::filter(fold.enrichment==0)
## head(res0)

res.not0 <- res %>% dplyr::filter(fold.enrichment>0)
head(res.not0)
max(res.not0$logFC)
min(res.not0$logFC)
## res %>% dplyr::filter(logFC < -3.93, logFC != -Inf)

res$logFC <- ifelse(res$logFC==-Inf, minlogFC, res$logFC) 
head(res)
nrow(res)
max(res$logFC)
min(res$logFC)

table(res$MCls, res$contrast)

res_MCls <- unique(res$MCls)
#res_MCls
res_contrast <- unique(res$contrast)
#res_contrast
res_comb <- unique(sort(res$comb))
#res_comb


res_MCls.comb <- rep(unique(res$MCls), each=length(unique(res$contrast)))
#res_MCls.comb
res_contrast.comb <- rep(unique(res$contrast), times=length(unique(res$MCls)))
#res_contrast.comb
res_total.comb <- length(unique(res$contrast))*length(unique(res$MCls))
#res_total.comb
res_dataset_res <- data.frame(MCls=res_MCls.comb, contrast=res_contrast.comb)
#res_dataset_res
res_dataset_res[1,1]

#res_vitE_0_CD4Naive <- res %>% dplyr::filter(res$contrast == "vitE",res$MCls=="0_CD4Naive", qvalue.fisher<0.1)
#head(res_vitE_0_CD4Naive)
#res_water_0_CD4Naive <- res %>% dplyr::filter(res$contrast == "water",res$MCls=="0_CD4Naive", qvalue.fisher<0.1)
#head(res_water_0_CD4Naive)


head(res)

#----sig_res_annotation-----
padj_cutoff <- 0.1
logFC_cutoff <- 0.5

##FC <- 1.41
##res_sig_anno <- res %>% mutate(sig=ifelse(qvalue.fisher < padj_cutoff | abs(fold.enrichment)>= FC, 1, 0))
####%>% dplyr::arrange(qvalue.fisher)
##nrow(res_sig_anno)


res_sig_anno <- res %>% mutate(sig=ifelse(qvalue.fisher <= padj_cutoff | abs(logFC)>= logFC_cutoff, 1, 0))
head(res_sig_anno)
nrow(res_sig_anno)

res_sig_anno_MCls <- unique(res_sig_anno$MCls)
#res_sig_anno_MCls
res_sig_anno_contrast <- unique(res_sig_anno$contrast)
#res_sig_anno_contrast
res_sig_anno_comb <- unique(sort(res_sig_anno$comb))
#res_sig_anno_comb


res_sig_anno_MCls.comb <- rep(unique(res_sig_anno$MCls), each=length(unique(res_sig_anno$contrast)))
#res_sig_anno_MCls.comb
res_sig_anno_contrast.comb <- rep(unique(res_sig_anno$contrast), times=length(unique(res_sig_anno$MCls)))
#res_sig_anno_contrast.comb
res_sig_anno_total.comb <- length(unique(res_sig_anno$contrast))*length(unique(res_sig_anno$MCls))
#res_sig_anno_total.comb
res_sig_anno_dataset <- data.frame(MCls=res_sig_anno_MCls.comb, contrast=res_sig_anno_contrast.comb)
#res_sig_anno_dataset
#res_sig_anno_dataset[1,1]



#----sig_res----
## sig_res <- res %>% dplyr::filter(qvalue.fisher < padj_cutoff, abs(fold.enrichment)> FC)
## #%>% dplyr::arrange(qvalue.fisher)
## head(sig_res)
## nrow(sig_res)
## #table(sig_res$MCls, sig_res$contrast)
## #table(sig_res$sig)

head(res)

#max(res$fold.enrichment)
#min(res$fold.enrichment)

sig_res <- res %>% dplyr::filter(qvalue.hyper < padj_cutoff, abs(logFC) > logFC_cutoff)
#%>% dplyr::arrange(qvalue.fisher)
#just checking
#sig_res <- res %>% dplyr::filter(logFC  < 0, fold.enrichment != 0) %>% dplyr::arrange(fold.enrichment)
head(sig_res)
nrow(sig_res)

sig_res_maxlogFC <- max(sig_res$logFC)
sig_res_maxlogFC
sig_res_minlogFC <- min(sig_res$logFC)
sig_res_minlogFC
table(sig_res$MCls, sig_res$contrast)

sig_res_MCls <- unique(sig_res$MCls)
#sig_res_MCls
sig_res_contrast <- unique(sig_res$contrast)
#sig_res_contrast
sig_res_comb <- unique(sort(sig_res$comb))
#sig_res_comb


sig_res_MCls.comb <- rep(unique(sig_res$MCls), each=length(unique(sig_res$contrast)))
#sig_res_MCls.comb
sig_res_contrast.comb <- rep(unique(sig_res$contrast), times=length(unique(sig_res$MCls)))
#sig_res_contrast.comb
sig_res_total.comb <- length(unique(sig_res$contrast))*length(unique(sig_res$MCls))
#sig_res_total.comb
sig_res_dataset <- data.frame(MCls=sig_res_MCls.comb, contrast=sig_res_contrast.comb)
#sig_res_dataset
#sig_res_dataset[1,1]


#sig_res_vitE_0_CD4Naive <- sig_res %>% dplyr::filter(sig_res$contrast == "vitE",sig_res$MCls=="0_CD4Naive")
#sig_res_vitE_0_CD4Naive

#write.csv(sig_res, paste0("./1_DiffRNA_2.outs/sig_genes.csv"), quote = FALSE, row.names = FALSE)
#sig_res <- read.csv("./1_DiffRNA_2.outs/sig_genes.csv")




## #----sig_res_split----
## #sig_res_split <- split(sig_res, f = list(sig_res$contrast,sig_res$MCls))
## sig_res_split <- split(sig_res, f = sig_res$comb)
## head(sig_res_split)
## length(sig_res_split)
## names(sig_res_split)
## #sig_res_split[["vitD.3_TEM"]]




## #----sig_res_split2 -> list----
## sig_res_split2 = list()
## for (i in names(sig_res_split)){
##     print(i)
##     sig_res_split2[[i]] <- sig_res_split[[i]]  %>% dplyr::arrange(qvalue.fisher) %>% dplyr::pull(motif.name)
## }

## ## for (i in names(sig_res_split)){
## ##     print(i)
## ##     sig_res_split2[[i]] <- sig_res_split[[i]]  %>% dplyr::arrange(qvalue.fisher) %>% dplyr::pull(motif.name) %>% head(n=20)
## ## }

## sig_res_split2

## #write.table(x, file, append = FALSE, sep = " ", dec = ".", row.names = TRUE, col.names = TRUE)
## #write.table(sig_res_split2, file="./1_motif.outs/sig_enrichment__macs2_0.1_cn.txt")


#----combine_data----
#head(sig_res)


## #######################################################
## #----dataframe & heatmap for significant annotation----
## #######################################################
## head(res_sig_anno)
## nrow(res_sig_anno)

## table(res_sig_anno$motif.name)
## res_sig_anno_comb

## sig_motifs <- res_sig_anno$motif
## names(sig_motifs) <- res_sig_anno$motif.name
## sig_motifs <- sig_motifs[!duplicated(sig_motifs)]
## head(sig_motifs)
## length(sig_motifs)

## mat <- lapply(res_sig_anno_comb, function(ii){
##    enrich2 <- res_sig_anno%>%dplyr::filter(comb==ii)
##    z <- enrich2$sig
## #   names(z) <- enrich2$motif
## #   z[motif]
##    names(z) <- enrich2$motif.name
##    z[names(sig_motifs)]
##    })

## head(mat)

## length(mat[[1]])
## #min(mat)
## #max(mat)

## #****
## mat <- do.call(cbind, mat)
## colnames(mat) <- res_sig_anno_comb

## head(mat)
## length(mat)
## nrow(mat)
## min(mat)
## max(mat)
## #colnames(mat)

## #rownames(mat) <- motif
## #head(mat)

## mat[!is.na(mat)] = 1
## head(mat)

## mat[is.na(mat)] = 0
## head(mat)


## column_ha <- HeatmapAnnotation(
## #    celltype=res_sig_anno$MCls,
## #    treatment=res_sig_anno$contrast,
##     celltype=rep(MCls, each=length(unique(res_sig_anno$contrast))),
##     treatment=rep(contrast, times=length(unique(res_sig_anno$MCls))),
##     col=list(celltype=c("4_Bcell"="#4daf4a", "6_Monocyte"="#984ea3",
##                     "2_NKcell"="#aa4b56", "0_CD4Naive"="#ffaa00", "1_TCM"="pink",
##                     "3_TEM"="blue", "5_CD8Naive"="green", "7_dnT"="black"),
##              treatment=c("caffeine"="red", "nicotine"="tan",
##                          "vitA"="tan4", "vitD"="seagreen4",
##                          "vitE"="salmon3", "water"="grey", "zinc"="maroon3")))


## fig <- Heatmap(mat,
## #   Rowv = NA,
## #   Colv = NA, 
##    cluster_rows=T, cluster_columns=F,
##    top_annotation=column_ha,
##    heatmap_legend_param=list(title="significance",
##       title_gp=gpar(fontsize=10),
##       labels_gp=gpar(fontsize=10)),
##    show_row_names=T, show_column_names=T,
##    row_names_gp=gpar(fontsize=2),
##    column_names_gp=gpar(fontsize=5),
##    raster_device="png")


## figfn <- "./1_motif.outs/Figure4.1_significant_heatmap_macs2_0.1_cn.png"
## png(figfn, height=6000, width=3000, res=350)
## #set.seed(0)
## fig <- draw(fig)
## dev.off()

   

## enrich%>%mutate(LFC=log2(fold.enrichment))%>%
##    dplyr::filter(qvalue.hyper<0.1, LFC>0.5)%>%
##    group_by(condition)%>%summarise(ny=n(),.groups="drop")

## head(enrich)


#################################################
#----dataframe & heatmap for only significant----
#################################################
## head(sig_res)
## ncol(sig_res)
## nrow(sig_res)
## #table(sig_res$motif.name)
## max(sig_res$logFC)
## min(sig_res$logFC)

## sig_res_comb

## sig_motifs <- sig_res$motif
## names(sig_motifs) <- sig_res$motif.name
## sig_motifs <- sig_motifs[!duplicated(sig_motifs)]
## head(sig_motifs)
## length(sig_motifs)

## mat <- lapply(sig_res_comb, function(ii){
##    enrich2 <- sig_res%>%dplyr::filter(comb==ii)
##    z <- enrich2$logFC
## #   names(z) <- enrich2$motif
## #   z[motif]
##    names(z) <- enrich2$motif.name
##    z[names(sig_motifs)]
##    })
## head(mat)
## length(mat[[1]])
## min(mat)
## max(mat)







## #****
## mat <- do.call(cbind, mat)
## colnames(mat) <- sig_res_comb
## rownames(mat) <- names(sig_motifs)
## head(mat)
## length(mat)
## nrow(mat)

## sig_res_maxlogFC
## sig_res_minlogFC
## ## mat_df <- as.data.frame(mat)
## ## #mat_df <- na.omit(mat_df)
## ## mat_df[is.na(mat)] = 0
## ## head(mat_df)
## ## max(mat_df)
## ## min(mat_df)


## #doing this will indicate non significants as downregulated
## #mat[is.na(mat)] = -max(mat_df)
## #head(mat)
## #max(mat)
## #min(mat)







##----mat_estimate----
sig_motifs <- sig_res$motif
names(sig_motifs) <- sig_res$motif.name
sig_motifs <- sig_motifs[!duplicated(sig_motifs)]
head(sig_motifs)
length(sig_motifs)


mat <- map_dfc(comb, function(ii){
    enrich2 <- res%>%dplyr::filter(comb==ii)
    z <- enrich2$logFC
    #res2 <- res %>% dplyr::filter(comb==ii)
    #b <- res2$estimate
    #names(b) <- res2$peak
    #b[sigdps.qval]
    names(z) <- enrich2$motif.name
    z[names(sig_motifs)]
})
print("Number of significant in mat")
nrow(mat)


mat <- as.matrix(mat)
rownames(mat) <- sigdps.qval #DP
colnames(mat) <- comb
head(mat)
max(mat)
min(mat)









#this is the reight way!!
#mat[!is.na(mat)] = 1
#head(mat)
#mat[is.na(mat)] = 1 #since this is a ratio and 1 indicates no change
mat[is.na(mat)] = 0  #if using for logFC
head(mat)
max(mat)
min(mat)


mat.ct <- mat
head(mat.ct)

#----annotation as cell type----
## colnames(mat)
## cvt0 <- str_split(colnames(mat), "_", simplify=T)
## cvt0
## cvt0[,1]

comb.celltype <- colnames(mat.ct)
comb.celltype

cvt.celltype <- as.data.frame(str_split(comb.celltype, "_", simplify=T)) %>% mutate(comb=paste0(V1, "_", V2))
cvt.celltype
#cvt.celltype[,1]

column_ha_ct <- HeatmapAnnotation(
    celltype=cvt.celltype[,1],
    treatment=cvt.celltype[,2],
#    celltype=rep(MCls, each=length(unique(sig_res$contrast))),
#    treatment=rep(contrast, times=length(unique(sig_res$MCls))),
    col=list(celltype=c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
                    "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00", "1-TCM"="pink",
                    "3-TEM"="blue", "5-CD8Naive"="green", "7-dnT"="black"),
             treatment=c("caffeine"="red", "nicotine"="tan",
                         "vitA"="tan4", "vitD"="seagreen4",
                         "vitE"="salmon3", "water"="grey", "zinc"="maroon3")))


#----annotation and data as treatment----
cvt.treat <- cvt.celltype[order(cvt.celltype$V2),]
cvt.treat

mat.treat <- mat.ct[, cvt.treat$comb]
head(mat.treat)

#annotation as treatment
column_ha_treat <- HeatmapAnnotation(
    celltype=cvt.treat[,1],
    treatment=cvt.treat[,2],
    col=list(celltype=c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
                    "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00", "1-TCM"="pink",
                    "3-TEM"="blue", "5-CD8Naive"="green", "7-dnT"="black"),
             treatment=c("caffeine"="red", "nicotine"="tan",
                         "vitA"="tan4", "vitD"="seagreen4",
                         "vitE"="salmon3", "water"="grey", "zinc"="maroon3")))


#----breaks and color----
#breaks <- c(seq(0,1,length.out=50), quantile(b2, probs=seq(0, 1, length.out=49)), 38.93) 
#breaks
#col_fun <-  colorRamp2(breaks,
#   colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100))


#breaks <- (seq(-6, 6, length.out=13)) 
#breaks


#col_fun <-  colorRamp2(breaks,
#   colorRampPalette(rev(brewer.pal(n=11, name="RdBu")))(13))


#col_fun <-  colorRampPalette(brewer.pal(n=9, name="Reds"))(100)

y <- as.numeric(mat.treat)
max(y)
abs(min(y))
#length(y)

scale <- quantile(abs(y),probs=0.99)
scale


mybreaks <- seq(-scale, scale, length.out=100)
mybreaks
names(mybreaks) <- NULL

col_fun <- colorRampPalette(rev(brewer.pal(n=11, name="RdBu")))(100)


#----figure as celltype----
fig <- Heatmap(mat.ct,
    col_fun,
#   Rowv = NA,
#   Colv = NA, 
   cluster_rows=T, cluster_columns=F,
   top_annotation=column_ha_ct,
   heatmap_legend_param=list(title="fold.enrichment",
      title_gp=gpar(fontsize=10),
      labels_gp=gpar(fontsize=10)),
   show_row_names=T, show_column_names=T,
   row_names_gp=gpar(fontsize=2),
   column_names_gp=gpar(fontsize=10),
   raster_device="png")


figfn <- "./1_motif.outs/Figure4.2_significant_heatmap_macs2_0.1_cn_celltype_control.png"
png(figfn, height=6000, width=3000, res=350)
#set.seed(0)
fig <- draw(fig)
dev.off()


#----figure as treat----
fig <- Heatmap(mat.treat,
    col_fun,
#   Rowv = NA,
#   Colv = NA, 
   cluster_rows=T, cluster_columns=F,
   top_annotation=column_ha_treat,
   heatmap_legend_param=list(title="fold.enrichment",
      title_gp=gpar(fontsize=10),
      labels_gp=gpar(fontsize=10)),
   show_row_names=T, show_column_names=T,
   row_names_gp=gpar(fontsize=2),
   column_names_gp=gpar(fontsize=10),
   raster_device="png")


figfn <- "./1_motif.outs/Figure4.2_significant_heatmap_macs2_0.1_cn_treat_control.png"
png(figfn, height=6000, width=3000, res=350)
#set.seed(0)
fig <- draw(fig)
dev.off()

   

enrich%>%mutate(LFC=log2(fold.enrichment))%>%
   dplyr::filter(qvalue.hyper<0.1, LFC>0.5)%>%
   group_by(condition)%>%summarise(ny=n(),.groups="drop")
head(enrich)




#######################################################
#----dataframe & heatmap for only significant top10----
#######################################################
##head(sig_res)
##table(sig_res$motif.name)

##sig_res_comb

## sig_motifs <- sig_res$motif
## names(sig_motifs) <- sig_res$motif.name
## sig_motifs <- sig_motifs[!duplicated(sig_motifs)]
## head(sig_motifs)
## length(sig_motifs)

## mat <- lapply(sig_res_comb, function(ii){
##    enrich2 <- sig_res%>%dplyr::filter(comb==ii) %>% dplyr::arrange(qvalue.fisher) %>% head(n=10)
##    z <- enrich2$logFC
## #   names(z) <- enrich2$motif
## #   z[motif]
##    names(z) <- enrich2$motif.name
##    z[names(sig_motifs)]
##    })
## head(mat)
## length(mat[[1]])
## #min(mat)
## #max(mat)

## #****
## mat <- do.call(cbind, mat)
## colnames(mat) <- sig_res_comb
## rownames(mat) <- names(sig_motifs)
## head(mat)
## length(mat)
## nrow(mat)

## #----sum filt----
## ii <- rowSums(is.na(mat))
## head(ii)
## length(ii)

## #EITHER
## mat2 <- mat[ii==0,]

## #OR
## mat2 <- mat[ii<length(sig_res_comb),]
## head(mat2)
## nrow(mat2)
## max(mat2)
## min(mat2)
## mean(mat2)






head(res)
nrow(res)

head(sig_res)
nrow(sig_res)

sig_res_comb

sig_res_filt <- map_dfr(sig_res_comb, function(ii){
    enrich2 <- sig_res %>% dplyr::filter(comb==ii)
    enrich2 <- enrich2 %>% dplyr::arrange(qvalue.fisher)
    enrich2 <- enrich2[1:10,]
})
head(sig_res_filt)
nrow(sig_res_filt)

sig_motifs <- sig_res_filt$motif
names(sig_motifs) <- sig_res_filt$motif.name
sig_motifs <- sig_motifs[!duplicated(sig_motifs)]
head(sig_motifs)
length(sig_motifs)

mat <- map_dfc(sig_res_comb, function(ii){
    enrich2 <- res %>% dplyr::filter(comb==ii)
    z <- enrich2$logFC
    names(z) <- enrich2$motif.name
    z[names(sig_motifs)]
})
    
head(mat)
length(mat[[1]])
#min(mat)
#max(mat)

#****
mat <- do.call(cbind, mat)
colnames(mat) <- sig_res_comb
rownames(mat) <- names(sig_motifs)
head(mat)
length(mat)
nrow(mat)

max(mat)
min(mat)

#----sum filt----
## ii <- rowSums(is.na(mat))
## head(ii)
## length(ii)

## #EITHER
## mat2 <- mat[ii==0,]

## #OR
## mat2 <- mat[ii<length(sig_res_comb),]
## head(mat2)
## nrow(mat2)
## max(mat2)
## min(mat2)
## mean(mat2)

#OR
## mat[is.na(mat)] = 0
## mat2 <- mat

## sig_res_comb.2 <- sort(sig_res_comb)
## mat3 <- mat2[, sig_res_comb.2]

## max(mat3)
## min(mat3)

mat2 <- mat

#--------
sig_res_maxlogFC
sig_res_minlogFC
## mat_df <- as.data.frame(mat)
## #mat_df <- na.omit(mat_df)
## mat_df[is.na(mat)] = 0
## head(mat_df)
## max(mat_df)
## min(mat_df)


#doing this will indicate non significants as downregulated
#mat[is.na(mat)] = -max(mat_df)
#head(mat)
#max(mat)
#min(mat)


#this is the reight way!!
#mat[!is.na(mat)] = 1
#head(mat)
#mat[is.na(mat)] = 1 #since this is a ratio and 1 indicates no change
mat2[is.na(mat2)] = 0  #if using for logFC
head(mat2)
nrow(mat2)
max(mat2)
min(mat2)


mat.ct <- mat2
head(mat.ct)

#----annotation as cell type----
## colnames(mat)
## cvt0 <- str_split(colnames(mat), "_", simplify=T)
## cvt0
## cvt0[,1]

comb.celltype <- colnames(mat.ct)
comb.celltype

cvt.celltype <- as.data.frame(str_split(comb.celltype, "_", simplify=T)) %>% mutate(comb=paste0(V1, "_", V2))
cvt.celltype
#cvt.celltype[,1]

column_ha_ct <- HeatmapAnnotation(
    celltype=cvt.celltype[,1],
    treatment=cvt.celltype[,2],
#    celltype=rep(MCls, each=length(unique(sig_res$contrast))),
#    treatment=rep(contrast, times=length(unique(sig_res$MCls))),
    col=list(celltype=c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
                    "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00", "1-TCM"="pink",
                    "3-TEM"="blue", "5-CD8Naive"="green", "7-dnT"="black"),
             treatment=c("caffeine"="red", "nicotine"="tan",
                         "vitA"="tan4", "vitD"="seagreen4",
                         "vitE"="salmon3", "water"="grey", "zinc"="maroon3")))


#----annotation and dataas treatment----
cvt.treat <- cvt.celltype[order(cvt.celltype$V2),]
cvt.treat

mat.treat <- mat.ct[, cvt.treat$comb]
head(mat.treat)

#annotation as treatment
column_ha_treat <- HeatmapAnnotation(
    celltype=cvt.treat[,1],
    treatment=cvt.treat[,2],
    col=list(celltype=c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
                    "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00", "1-TCM"="pink",
                    "3-TEM"="blue", "5-CD8Naive"="green", "7-dnT"="black"),
             treatment=c("caffeine"="red", "nicotine"="tan",
                         "vitA"="tan4", "vitD"="seagreen4",
                         "vitE"="salmon3", "water"="grey", "zinc"="maroon3")))


#----breaks and color----
#breaks <- c(seq(0,1,length.out=50), quantile(b2, probs=seq(0, 1, length.out=49)), 38.93) 
#breaks
#col_fun <-  colorRamp2(breaks,
#   colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100))


#breaks <- (seq(-6, 6, length.out=13)) 
#breaks


#col_fun <-  colorRamp2(breaks,
#   colorRampPalette(rev(brewer.pal(n=11, name="RdBu")))(13))

##col_fun <-  colorRampPalette(brewer.pal(n=9, name="Reds"))(100)


y <- as.numeric(mat.ct)
length(y)
max(y)
min(y)


scale <- quantile(abs(y),probs=0.99)
scale
length(y[abs(y)>scale])

mybreaks <- seq(-scale, scale, length.out=100)
mybreaks
names(mybreaks) <- NULL

col_fun <- colorRampPalette(rev(brewer.pal(n=11, name="RdBu")))(100)



#----figure as celltype----
fig <- Heatmap(mat.ct,
    col_fun,
#   Rowv = NA,
#   Colv = NA, 
   cluster_rows=T, cluster_columns=F,
   top_annotation=column_ha_ct,
   heatmap_legend_param=list(title="fold.enrichment",
      title_gp=gpar(fontsize=10),
      labels_gp=gpar(fontsize=10)),
   show_row_names=T, show_column_names=T,
   row_names_gp=gpar(fontsize=5),
   column_names_gp=gpar(fontsize=10),
   raster_device="png")


figfn <- "./1_motif.outs/Figure5.1_top10_heatmap_macs2_0.1_cn_celltype_control.png"
png(figfn, height=6000, width=3000, res=350)
#set.seed(0)
fig <- draw(fig)
dev.off()


#----figure as treat----
fig <- Heatmap(mat.treat,
    col_fun,
#   Rowv = NA,
#   Colv = NA, 
   cluster_rows=T, cluster_columns=F,
   top_annotation=column_ha_treat,
   heatmap_legend_param=list(title="fold.enrichment",
      title_gp=gpar(fontsize=10),
      labels_gp=gpar(fontsize=10)),
   show_row_names=T, show_column_names=T,
   row_names_gp=gpar(fontsize=5),
   column_names_gp=gpar(fontsize=10),
   raster_device="png")


figfn <- "./1_motif.outs/Figure5.1_top10_heatmap_macs2_0.1_cn_treat_control.png"
png(figfn, height=6000, width=3000, res=350)
#set.seed(0)
fig <- draw(fig)
dev.off()

   

enrich%>%mutate(LFC=log2(fold.enrichment))%>%
   dplyr::filter(qvalue.hyper<0.1, LFC>0.5)%>%
   group_by(condition)%>%summarise(ny=n(),.groups="drop")
head(enrich)








########################################
#----sig_res_unique----
########################################
sig_res_comb

sig_res_unique <- map_dfr(sig_res_comb, function(oneX){
          time0 <- Sys.time()
                  ###
                  res.2 <- sig_res%>%dplyr::filter(comb==oneX)
                  ##
                  res.3 <- sig_res%>%dplyr::filter(comb!=oneX)
                      #mutate(count=ifelse(nrow(cvt.2)>=8,1,0))
                  ##
                  other.motif.names <- unique(res.3$motif.name)
                  res.2.2 <- res.2%>%dplyr::filter(!motif.name %in% other.motif.names)
                  time1 <- Sys.time()
                  elapsed <- difftime(time1, time0, units="mins")
                  cat(oneX, elapsed, "Done\n")
                  res.2.2
                })
head(sig_res_unique)

table(sig_res_unique$motif.name, sig_res_unique$contrast, sig_res_unique$MCls)

sig_res_unique_split <- split(sig_res_unique, f = sig_res_unique$comb)
head(sig_res_unique_split)
length(sig_res_unique_split)
names(sig_res_unique_split)

sig_res_unique_split2 = list()
for (i in names(sig_res_unique_split)){
    print(i)
    sig_res_unique_split2[[i]] <- sig_res_unique_split[[i]]  %>% dplyr::arrange(qvalue.fisher) %>% dplyr::pull(motif.name)
}
sig_res_unique_split2



#MCls
sig_res_unique_MCls <- map_dfr(sig_res_MCls, function(oneX){
          time0 <- Sys.time()
                  ###
                  res.2 <- sig_res%>%dplyr::filter(MCls==oneX)
                  ##
                  res.3 <- sig_res%>%dplyr::filter(MCls!=oneX)
                      #mutate(count=ifelse(nrow(cvt.2)>=8,1,0))
                  ##
                  other.motif.names <- unique(res.3$motif.name)
                  res.2.2 <- res.2%>%dplyr::filter(!motif.name %in% other.motif.names)
                  time1 <- Sys.time()
                  elapsed <- difftime(time1, time0, units="mins")
                  cat(oneX, elapsed, "Done\n")
                  res.2.2
                })
head(sig_res_unique_MCls)

table(sig_res_unique_MCls$motif.name, sig_res_unique_MCls$contrast, sig_res_unique_MCls$MCls)




#contrast
sig_res_unique_contrast <- map_dfr(sig_res_contrast, function(oneX){
          time0 <- Sys.time()
                  ###
                  res.2 <- sig_res%>%dplyr::filter(contrast==oneX)
                  ##
                  res.3 <- sig_res%>%dplyr::filter(contrast!=oneX)
                      #mutate(count=ifelse(nrow(cvt.2)>=8,1,0))
                  ##
                  other.motif.names <- unique(res.3$motif.name)
                  res.2.2 <- res.2%>%dplyr::filter(!motif.name %in% other.motif.names)
                  time1 <- Sys.time()
                  elapsed <- difftime(time1, time0, units="mins")
                  cat(oneX, elapsed, "Done\n")
                  res.2.2
                })
head(sig_res_unique_contrast)

table(sig_res_unique_contrast$motif.name, sig_res_unique_contrast$MCls, sig_res_unique_contrast$contrast)




end

#--------------------------JULONG SCRIPT---------------------------------
##
## fn <- "../2_Differential/2.2_compareRNAandATAC.outs/2_annot.ChIPseeker.rds"
## peakAnno <- read_rds(fn)%>%
##    as.data.frame()%>%
##    mutate(chr=gsub("chr", "", seqnames),
##                    peak_region=paste(chr,start,end,sep="-"))%>%
##    dplyr::select(peak_region, geneId, SYMBOL, distanceToTSS)

## ###
## resDP2 <- resDP2%>%left_join(peakAnno, by=c("gene"="peak_region"))

## top.DP <- as.character(resDP2$gene)

### test enrichment
## build a set of background peaks
## find peaks open in T cell
## atac2 <- subset(atac, subset=MCls=="Tcell")
## cell0 <- Cells(atac2)
## open.peaks <- AccessiblePeaks(atac, cells=cell0)
## ## match the overall GC content in the peak set
## meta.feature <- GetAssayData(atac, assay="ATAC", slot="meta.features")
## peak.matched <- MatchRegionStats(
##    meta.feature=meta.feature[open.peaks,],
##    query.feature=meta.feature[top.DP,],
##    n=15000)

## ### test enrichment
## enriched <- FindMotifs(
##    object=atac,
##    features=top.DP,
##    background=peak.matched)

##




####################
### heatmap
####################
rm(list=ls())

#enrich <- read_rds("./1_motif.outs/2_motif.enrich_macs2_0.1_cn.rds")

enrich <- enriched.motif
head(enrich)
nrow(enrich)

enrich <- enrich%>%
   drop_na(fold.enrichment)%>%
   mutate(condition=paste(MCls, contrast, sep="_"))
head(enrich)
#nrow(enrich%>%dplyr::filter(fold.enrichment<3))

condition <- unique(enrich$condition)
condition

#----build matrix for heatmap----
motif <- enrich$motif
head(motif)
length(motif)

names(motif) <- enrich$motif.name
head(motif)
#length(motif)

motif <- motif[!duplicated(motif)]
length(motif)


#----individual step----
## enrich2 <- enrich%>%dplyr::filter(condition=="7_dnT_caffeine")
## head(enrich2)
## nrow(enrich2)

## z <- enrich2$fold.enrichment
## z

## names(z) <- enrich2$motif
## z
## length(z)

## z[motif] 
## length(z[motif])


#----function----
mat <- lapply(condition, function(ii){
   enrich2 <- enrich%>%dplyr::filter(condition==ii)
   z <- enrich2$fold.enrichment
#   names(z) <- enrich2$motif
#   z[motif]
   names(z) <- enrich2$motif.name
   z[names(motif)]
   })

head(mat)
length(mat[[1]])
#min(mat)
#max(mat)

#--------
mat <- do.call(cbind, mat)
colnames(mat) <- condition
head(mat)
length(mat)
min(mat)
max(mat)
#colnames(mat)

#rownames(mat) <- motif
#head(mat)

## sort(mat[,"7_dnT_caffeine"])
## head(sort(mat[,"7_dnT_caffeine"]))
## length(mat[,"7_dnT_caffeine"])

###
b <- as.vector(mat)
#head(b)
#length(b)

b2 <- b[b>1&b<2.66]
#head(b2)
#length(b2)
#min(b2)
#max(b2)

#seq(0,1,length.out=50)
#quantile(b2, probs=seq(0, 1, length.out=49))

breaks <- c(seq(0,1,length.out=50), quantile(b2, probs=seq(0, 1, length.out=49)), 38.93) 
breaks


col_fun <-  colorRamp2(breaks,
   colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100))
#head(col_fun)

#head(enrich)
#head(motif)

MCls <- unique(enrich$MCls)
MCls

treats <- unique(enrich$contrast)
treats

#done above
#condition <- unique(enrich$condition)
#condition


## column_ha <- HeatmapAnnotation(
##     celltype=c(rep("Bcell",each=3),rep(c("Monocyte", "NKcell", "Tcell"),each=4)),
##     treatment=c("LPS-DEX","PHA","PHA-DEX",
##                rep(c("LPS", "LPS-DEX", "PHA", "PHA-DEX"),times=3)),
##     col=list(celltype=c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
##                     "NKcell"="#aa4b56", "Tcell"="#ffaa00"),
##              treatment=c("LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
##                          "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")))

column_ha <- HeatmapAnnotation(
    celltype=rep(MCls, each=length(unique(enrich$contrast))),
    treatment=rep(treats, times=length(unique(enrich$MCls))),
    col=list(celltype=c("4_Bcell"="#4daf4a", "6_Monocyte"="#984ea3",
                    "2_NKcell"="#aa4b56", "0_CD4Naive"="#ffaa00", "1_TCM"="pink",
                    "3_TEM"="blue", "5_CD8Naive"="green", "7_dnT"="black"),
             treatment=c("caffeine"="red", "nicotine"="tan",
                         "vitA"="tan4", "vitD"="seagreen4",
                         "vitE"="salmon3", "water"="grey", "zinc"="maroon3")))


## fig <- Heatmap(mat, col=col_fun,
##    cluster_rows=T, cluster_columns=T,
##    top_annotation=column_ha,
##    heatmap_legend_param=list(title="fold.enrichment",
##       title_gp=gpar(fontsize=10),
##       labels_gp=gpar(fontsize=10)),
##    show_row_names=F, show_column_names=T,
##    column_names_gp=gpar(fontsize=5),
##    raster_device="png")


## figfn <- "./1_motif.outs/Figure2.1_heatmap_macs2_0.1_cn_3.png"
## png(figfn, height=800, width=600, res=125)
## #set.seed(0)
## fig <- draw(fig)
## dev.off()


fig <- Heatmap(mat,
   cluster_rows=T, cluster_columns=T,
   top_annotation=column_ha,
   heatmap_legend_param=list(title="fold.enrichment",
      title_gp=gpar(fontsize=10),
      labels_gp=gpar(fontsize=10)),
   show_row_names=T, show_column_names=T,
   row_names_gp=gpar(fontsize=2),
   column_names_gp=gpar(fontsize=5),
   raster_device="png")


figfn <- "./1_motif.outs/Figure2.1_heatmap_macs2_0.1_cn_3.png"
png(figfn, height=6000, width=3000, res=350)
#set.seed(0)
fig <- draw(fig)
dev.off()

   

enrich%>%mutate(LFC=log2(fold.enrichment))%>%
   dplyr::filter(qvalue.hyper<0.1, LFC>0.5)%>%
   group_by(condition)%>%summarise(ny=n(),.groups="drop")

head(enrich)




################
### qq plots ###
################
res <- read_rds("./1_motif.outs/2_motif.enrich_macs2_0.1_cn_control.rds")%>%as.data.frame() %>% mutate(condition=paste(MCls, contrast, sep="_"))

head(res)
nrow(res)

#condition <- sort(unique(res$condition))
#condition

res <- res %>% drop_na(pvalue)
head(res)
nrow(res)

condition <- sort(unique(res$condition))
condition

dfNew <- map_dfr(condition, function(ii){
  res2 <- res%>%dplyr::filter(condition==ii)
  ngene <- nrow(res2)
  res2 <- res2%>%
      arrange(pvalue)%>%
      mutate(observed=-log10(pvalue), expected=-log10(ppoints(ngene)))
  res2
})

head(dfNew)

table(dfNew$contrast)

lab1 <- c("caffeine"="caffeine", "nicotine"="nicotine", "zinc"="zinc",
          "vitA"="vitA", "vitD"="vitD", "vitE"="vitE", "water"="water")

## lab2 <- c("Bcell"="B cell", "Monocyte"= "Monocyte", "NKcell"="NK cell",
##    "Tcell"="T cell")

lab2 <- c("4_Bcell"="4_Bcell", "6_Monocyte"="6_Monocyte",
          "2_NKcell"="2_NKcell", "0_CD4Naive"="0_CD4Naive", "1_TCM"="1_TCM",
          "3_TEM"="3_TEM", "5_CD8Naive"="5_CD8Naive", "7_dnT"="7_dnT")

p <- ggplot(dfNew, aes(x=expected,y=observed))+
   ggrastr::rasterise(geom_point(size=0.3, color="grey40"),dpi=300)+
   geom_abline(colour="red")+
   facet_grid(MCls~contrast, scales="free",
      labeller=labeller(contrast=lab1, MCls=lab2))+
   xlab(bquote("Expected"~ -log[10]~"("~italic(plain(P))~")"))+
   ylab(bquote("Observed"~-log[10]~"("~italic(plain(P))~")"))+
   theme_bw()+
   theme(strip.text=element_text(size=12))

figfn <- "./1_motif.outs/Figure2.2_qq_macs2_0.1_cn_control.png"
png(figfn, width=1600, height=1600, res=175)
print(p)
dev.off()







############################################################################
############################################################################
### motif enrichment analysis directionally, up and down ###
############################################################################
############################################################################
#----MOTIF FILE SAVED IN THIS SCRIPT----
#fn <- "./1_motif.outs/1_scATAC.motif.rds" 
#atac <- read_rds(fn)
## motif <- Motifs(atac)
## pfm <- GetMotifData(object=motif, slot="pwm")
## motif <- SetMotifData(object=motif, slot="pwm", new.data=pfm)

#----RESDP FILE----
## differential peaks

fn <- "../2_Differential/1.2_DiffPeak.outs/2.0_DESeq.results_0.1_cn_sampleID+new_treat.rds"
resDP <- read_rds(fn) %>% as.data.frame() %>% mutate(direction=ifelse(estimate>0, 1, 0))

resDP <- resDP %>% as.data.frame() %>% mutate(direction=ifelse(estimate>0, 1, 0))
head(resDP)

### fisher test
cal.fisher <- function(df){
  ###  
  resfisher <- map_dfr(1:nrow(df),function(i){  
     dmat <- matrix(as.numeric(df[i,]),2,2)
     colnames(dmat) <- c("interest", "not.interest")
     rownames(dmat) <- c("in.motif", "not.motif")
     res <- fisher.test(dmat, alternative="greater")
     res2 <- data.frame(odds=res$estimate, pval.fisher=res$p.value)
     res2
  })
  resfisher
}

    
#MCls <- rep(c("Bcell", "Monocyte", "NKcell", "Tcell"), each=4)
#contrast <- rep(c("LPS", "LPS-DEX", "PHA", "PHA-DEX"), times=4)

MCls <- rep(unique(resDP$MCls), each=length(unique(resDP$contrast)))
#MCls
contrast <- rep(unique(resDP$contrast), times=length(unique(resDP$MCls)))
#contrast
total.comb <- length(unique(resDP$contrast))*length(unique(resDP$MCls))
#total.comb
dataset <- data.frame(MCls=MCls, contrast=contrast)
head(dataset)
dataset[1,1]




###enrichment motif analysis
enriched.motif <- lapply(1:total.comb, function(i){
###
    cell0 <- dataset[i,1]
    contrast0 <- dataset[i,2]
    cat(i, cell0, contrast0, "\n") 
###
    enrich <- lapply(c(0,1),function(ii){ 
        res2 <- resDP%>%
            dplyr::filter(MCls==cell0, contrast==contrast0)
        top.DP <- res2%>%
            dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5, direction==ii)%>%
            dplyr::pull(gene)%>%as.character()
        bg.DP <- res2%>%dplyr::pull(gene)%>%as.character() 
        n.interest <- length(top.DP)
        n.not <- length(bg.DP)-n.interest
### 
        if(n.interest>5){
            enrich2 <- FindMotifs(
                object=atac,
                features=top.DP,
                background=bg.DP)
            df <- data.frame("interest.in.motif"=enrich2$observed,
                             "interest.not.motif"=n.interest-enrich2$observed,
                             "not.interest.in.motif"=enrich2$background-enrich2$observed)
            df <- df%>%mutate("not.interest.not.motif"=n.not-not.interest.in.motif)
            fisher <- cal.fisher(df)
            enrich2 <- cbind(enrich2, fisher)
            enrich2$qvalue.hyper <- p.adjust(enrich2$pvalue,"BH")
            enrich2$qvalue.fisher <- p.adjust(enrich2$pval.fisher,"BH")
            enrich2$MCls <- cell0
            enrich2$contrast <- contrast0
            enrich2$direction <- ii
        }else{
            enrich2 <- NA   
        }
        enrich2
    })    
   ##return results 
    if ( sum(is.na(enrich))==2){
        enrich <- NA
    }else{   
        enrich <-enrich[!is.na(enrich)]
        enrich <- do.call(rbind,enrich)
    }   
    enrich 
})

head(enriched.motif)

enriched.motif <- enriched.motif[!is.na(enriched.motif)]
enriched.motif <- do.call(rbind,enriched.motif)
head(enriched.motif)

opfn <- "./1_motif.outs/3_motif.enrich.direction_macs2_0.1_cn_control.rds"
write_rds(enriched.motif, opfn)

enriched.motif <- read_rds("./1_motif.outs/3_motif.enrich.direction_macs2_0.1_cn.rds")
head(enriched.motif)



######################
###    heatmap     ###
######################
rm(list=ls())

direction2 <- c("0"="Down","1"="Up")

enrich <- read_rds("./1_motif.outs/3_motif.enrich.direction_macs2_0.1_cn.rds")
#enrich <- enriched.motif

enrich <- enrich%>%
   drop_na(fold.enrichment)%>%
   mutate(dir2=direction2[as.character(direction)],
          condition=paste(MCls, contrast, dir2, sep="_"))

head(enrich)

condition <- unique(enrich$condition)
condition

#----build matrix for heatmap----
motif <- enrich$motif
names(motif) <- enrich$motif.name
motif <- motif[!duplicated(motif)]


###
mat <- lapply(condition, function(ii){
   enrich2 <- enrich%>%dplyr::filter(condition==ii)
   z <- enrich2$fold.enrichment
#   names(z) <- enrich2$motif
#   z[motif]
   names(z) <- enrich2$motif.name
   z[names(motif)]
})

mat <- do.call(cbind, mat)
colnames(mat) <- condition
#head(mat)

#----not using this step----
#rownames(mat) <- motif
#----not using this step----

###
b <- as.vector(mat)

b2 <- b[b>1&b<2.62]

breaks <- c(seq(0, 1, length.out=50),
            quantile(b2, probs=seq(0, 1, length.out=49)),7.71) 

col_fun <-  colorRamp2(breaks,
   colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100))

## column_ha <- HeatmapAnnotation(
##     celltype=gsub("_.*", "", condition),
##     treatment=gsub(".*cell_|.*cyte_|_(D|U).*", "", condition),
##     col=list(celltype=c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
##                     "NKcell"="#aa4b56", "Tcell"="#ffaa00"),
##              treatment=c("LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
##                         "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")))

MCls <- unique(enrich$MCls)
MCls

treats <- unique(enrich$contrast)
treats

#done above
## condition <- unique(enrich$condition)
## condition

celltype <- str_split(condition, "_", simplify=T)
celltype

celltype.2 <- paste(celltype[,1], "_", celltype[,2], sep="")
celltype.2

#treatment=gsub(".*cell_|.*cyte_|_(D|U).*", "", condition)
#treatment
treatment.2 <- celltype[,3]
treatment.2

column_ha <- HeatmapAnnotation(
    celltype=celltype.2,
    treatment=treatment.2,
    col=list(celltype=c("4_Bcell"="#4daf4a", "6_Monocyte"="#984ea3",
                    "2_NKcell"="#aa4b56", "0_CD4Naive"="#ffaa00", "1_TCM"="pink",
                    "3_TEM"="blue", "5_CD8Naive"="green", "7_dnT"="black"),
             treatment=c("caffeine"="red", "nicotine"="tan",
                         "vitA"="tan4", "vitD"="seagreen4",
                         "vitE"="salmon3", "water"="grey", "zinc"="maroon3"))
)


## fig <- Heatmap(mat, col=col_fun,
##    cluster_rows=T, cluster_columns=F,
##    top_annotation=column_ha,
##    heatmap_legend_param=list(title="fold.enrichment",
##       title_gp=gpar(fontsize=10),
##       labels_gp=gpar(fontsize=10)),
##    show_row_names=F, show_column_names=T,
##    column_names_gp=gpar(fontsize=5),
##    raster_device="png")

## figfn <- "./1_motif.outs/Figure3.1_heatmap_macs2_0.1_cn.png"
## png(figfn, height=800, width=700, res=120)
## set.seed(0)
## fig <- draw(fig)
## dev.off()

fig <- Heatmap(mat,
   cluster_rows=T, cluster_columns=F,
   top_annotation=column_ha,
   heatmap_legend_param=list(title="fold.enrichment",
      title_gp=gpar(fontsize=10),
      labels_gp=gpar(fontsize=10)),
   show_row_names=F, show_column_names=T,
   column_names_gp=gpar(fontsize=5),
   raster_device="png")

figfn <- "./1_motif.outs/Figure3.1_heatmap_macs2_0.1_cn_2.png"
png(figfn, height=3200, width=2800, res=225)
#set.seed(0)
fig <- draw(fig)
dev.off()



enrich%>%mutate(LFC=log2(fold.enrichment))%>%
   dplyr::filter(qvalue.hyper<0.1, LFC>0.5)%>%
   group_by(condition)%>%summarise(ny=n(),.groups="drop")




##############################
### barplot enriched motif ###
##############################
direction2 <- c("0"="Down","1"="Up")

enrich <- read_rds("./1_motif.outs/3_motif.enrich.direction_macs2_0.1_cn_control.rds")
head(enrich)
nrow(enrich)
##table(enrich$MCls)
##table(enrich$contrast)

enrich <- enrich%>%
   drop_na(fold.enrichment)%>%
   mutate(LFC=log2(fold.enrichment))

###
res <- enrich%>%dplyr::filter(qvalue.hyper<0.1, fold.enrichment>1.41)

sigs <- res%>%group_by(MCls, contrast, direction) %>% summarise(ny=n(), .groups="drop")
head(sigs)

###
sig4 <- sigs%>%mutate(ny2=ifelse(direction==0, -ny, ny))
sig4

breaks_value <- pretty(c(-600, 600), 5)
breaks_value

## p <- ggplot(sig4, aes(x=MCls, y=ny2))+
##    geom_bar(aes(fill=factor(MCls), alpha=factor(direction)), stat="identity")+
##    scale_fill_manual(values=c("Bcell"="#4daf4a",
##       "Monocyte"="#984ea3", "NKcell"="#aa4b56", "Tcell"="#ffaa00"))+
##    scale_alpha_manual(values=c("0"=0.5, "1"=1))+
##    geom_hline(yintercept=0, color="grey60")+
##    geom_text(aes(x=MCls, y=ny2, label=abs(ny2),
##       vjust=ifelse(direction==1, -0.2, 1.2)), size=3)+
##    scale_y_continuous("", breaks=breaks_value, limits=c(-600,650),
##                       labels=abs(breaks_value))+
##    facet_grid(~contrast,
##       labeller=labeller(contrast=c("LPS"="LPS","LPS-DEX"="LPS+DEX",
##                                    "PHA"="PHA", "PHA-DEX"="PHA+DEX")))+
##    theme_bw()+
##    theme(legend.position="none",
##          axis.title.x=element_blank(),
##          axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5))


p <- ggplot(sig4, aes(x=MCls, y=ny2))+
    geom_bar(aes(fill=factor(MCls), alpha=factor(direction)), stat="identity")+
    scale_fill_manual(values=c("4_Bcell"="#4daf4a", "6_Monocyte"="#984ea3",
                               "2_NKcell"="#aa4b56", "0_CD4Naive"="#ffaa00",
                               "1_TCM"="#ffaa00", "3_TEM"="#ffaa00",
                               "5_CD8Naive"="#ffaa00", "7_dnT"="#ffaa00"))+
    scale_alpha_manual(values=c("0"=0.5, "1"=1))+
    geom_hline(yintercept=0, color="grey60")+
    geom_text(aes(x=MCls, y=ny2, label=abs(ny2),
                  vjust=ifelse(direction==1, -0.2, 1.2)), size=3)+
    scale_y_continuous("", breaks=breaks_value, limits=c(-600,650),
                       labels=abs(breaks_value))+
    facet_grid(~contrast,
               ## labeller=labeller(contrast=c("4_Bcell"="4_Bcell",
               ##                              "6_Monocyte"="6_Monocyte",
               ##                              "2_NKcell"="2_NKcell",
               ##                              "0_CD4Naive"="0_CD4Naive",
               ##                              "1_TCM"="1_TCM",
               ##                              "3_TEM"="3_TEM",
               ##                              "5_CD8Naive"="5_CD8Naive",
               ##                              "7_dnT"="7_dnT")))+
               labeller=labeller(contrast=c("caffeine"="caffeine",
                                            "nicotine"="nicotine",
                                            "zinc"="zinc",
                                            "vitA"="vitA",
                                            "vitD"="vitD",
                                            "vitE"="vitE")))+
    theme_bw()+
    theme(legend.position="none",
          axis.title.x=element_blank(),
          axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5))


###
figfn <- "./1_motif.outs/Figure3.2_barplot_macs2_0.1_cn_control.png"
png(filename=figfn, width=3000, height=500, pointsize=24, res=225)
print(p)
dev.off()
                      







#######################################
### dot plots show motif enrichment ###
#######################################
#----PREPARE DATA FOR DOT PLOT----
#----ENRICHED.MOTIF----
enrich <- read_rds("./1_motif.outs/3_motif.enrich.direction_macs2_0.1_cn.rds")
head(enrich)

#### label
## clst <- paste(rep(c("LPS", "PHA", "LPS-DEX", "PHA-DEX"), each=4),
##    rep(c("Bcell", "Monocyte", "NKcell", "Tcell"), times=4), sep=".")


MCls <- unique(enrich$MCls)
MCls
length(MCls)

treats <- unique(enrich$contrast)
treats
length(treats)

## enrich <- enrich%>%
##    drop_na(fold.enrichment)%>%
##    mutate(dir2=direction2[as.character(direction)],
##           condition=paste(MCls, contrast, dir2, sep="_"))

#condition <- unique(enrich$condition)
#condition


clst <- paste(rep(unique(enrich$contrast), each=length(unique(enrich$MCls))),
   rep(unique(enrich$MCls), times=length(unique(enrich$contrast))), sep=".")
clst
length(clst)

clst <- paste(rep(clst,times=2), rep(c("Up","Down"), each=length(clst)), sep=".")
clst
length(clst)


###
ii <- paste(rep(c("A","B","C","D", "E", "F", "G"),each=length(unique(enrich$MCls))),rep(1:length(unique(enrich$MCls)),times=length(unique(enrich$contrast))), sep="")
ii

clst2 <- paste(rep(c("X","Y"),each=length(clst)/2), rep(ii,times=2), sep=".")
clst2

names(clst2) <- clst
clst2

lab2 <- gsub("-", "+", clst)
lab2

names(lab2) <- clst2 
lab2


###
drt2 <- c("0"="Down", "1"="Up")

enrich <- enrich%>%mutate(direction2=drt2[as.character(direction)],
   cluster=paste(contrast, MCls, direction2, sep="."),
   newCluster=clst2[cluster])
head(enrich)


###
## topmotif <- enrich%>%
##    dplyr::filter(fold.enrichment>1.41, qvalue.fisher<0.1)%>%dplyr::pull(motif)%>%unique()
##    ## group_by(cluster)%>%
##    ## top_n(n=5, wt=fold.enrichment)%>%ungroup()%>%as.data.frame()
## head(topmotif)
## length(topmotif)

## head(enrich2)
## topmotif <- enrich2%>%dplyr::pull(motif)
## head(topmotif)

## enrich3 <- enrich%>%
##    dplyr::filter(fold.enrichment>1.41, qvalue.fisher<0.1)
## head(enrich3)
## nrow(enrich3)


topmotif <- enrich%>%
       ## dplyr::filter(fold.enrichment>1.41, qvalue.fisher<0.1)%>%
       ## dplyr::pull(motif)%>%unique()
       group_by(cluster)%>%
       top_n(n=6, wt=fold.enrichment)%>%ungroup()%>%dplyr::pull(motif)%>%unique()
head(topmotif)

enrich3 <- enrich%>%
   dplyr::filter(motif%in%topmotif, fold.enrichment>1.41, qvalue.fisher<0.1)
head(enrich3)
nrow(enrich3)




#----MOTIF ACTIVITIES/CHROMVAR----
res <- read_rds("./2_motif.activities.outs/3_motif.diff.results_macs2_0.1_cn.rds")
head(res)


#motif <- read_rds("./2_motif.activities.outs/2_motif.ave_macs2_0.1_cn.rds")
#head(motif)

topmotif2 <- res%>%dplyr::filter(qval<0.1, abs(beta)>1.41)%>%dplyr::pull(motif)%>%unique()
head(topmotif2)
length(topmotif2)





#----PLOT----
p <- ggplot(enrich3, aes(x=newCluster, y=motif.name))+
   geom_point(aes(size=fold.enrichment, colour=qvalue.fisher))+
   scale_colour_gradient(name="FDR",
      low="blue", high="red", na.value=NA, trans="reverse", n.breaks=5,
      guide=guide_colourbar(order=1))+
   scale_size_binned("fold enrichment",
      guide=guide_bins(show.limits=T, axis=T,
          axis.show=arrow(length=unit(1.5,"mm"), ends="both"), order=2),
      n.breaks=4)+
   scale_x_discrete(labels=lab2)+ 
   theme_bw()+
   theme(axis.title=element_blank(),
         axis.text.x=element_text(angle=60, hjust=1, size=10),
         axis.text.y=element_text(size=8))

###
figfn <- "./1_motif.outs/Figure3.3_dotplot_macs2_0.1_cn.png"
png(figfn, width=4000, height=4000, res=250)
print(p)
dev.off()




################
### qq plots ###
################
drt2 <- c("0"="Down", "1"="Up")

res <- read_rds("./1_motif.outs/3_motif.enrich.direction_macs2_0.1_cn.rds")%>%
   as.data.frame()%>%
   drop_na(pvalue)%>%
   mutate(direction2=drt2[as.character(direction)],
          comb=paste(MCls, contrast, direction2, sep="_"))
head(res)

comb <- sort(unique(res$comb))
comb

dfNew <- map_dfr(comb, function(ii){
  res2 <- res%>%dplyr::filter(comb==ii)
  ngene <- nrow(res2)
  res2 <- res2%>%
     arrange(pvalue)%>%
      mutate(observed=-log10(pvalue), expected=-log10(ppoints(ngene)))
  res2
})

head(dfNew)

unique(dfNew$MCls)

unique(dfNew$contrast)

unique(dfNew$comb)

vars(comb)

p <- ggplot(dfNew, aes(x=expected,y=observed))+
   ggrastr::rasterise(geom_point(size=0.3, color="grey40"),dpi=300)+
   geom_abline(colour="red")+
   facet_wrap(vars(comb), scales="free", nrow=112, ncol=7)+
   xlab(bquote("Expected"~ -log[10]~"("~italic(plain(P))~")"))+
   ylab(bquote("Observed"~-log[10]~"("~italic(plain(P))~")"))+
   theme_bw()+
   theme(axis.text=element_text(size=7),
         strip.text=element_text(size=8))

figfn <- "./1_motif.outs/Figure3.4_qq_macs2_0.1_cn.png"
png(figfn, width=3200, height=3600, res=225)
print(p)
dev.off()




