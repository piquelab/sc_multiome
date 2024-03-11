###
##
library(Matrix)
library(tidyverse)
library(data.table)
library(irlba)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(Signac)
library(SeuratWrappers)
library(SeuratObject)
library(GenomicRanges)


##
library(cowplot)
library(RColorBrewer)
library(ComplexHeatmap)
library(viridis)
library(circlize)
library(ggrepel)
library(ggrastr)
library(openxlsx)

rm(list=ls())





###########################
### Plots 3D,  examples ###
###########################

outdir2 <- "./Plots_pub/3_compare_DEG_DARs_plots/"
if ( !file.exists(outdir2) ) dir.create(outdir2, showWarnings=F, recursive=T)

 


###
###
## fn <- "./sc_multiome_data_rs/1_seurat_combined_all_data.rds"
fn <- "./sc_multiome_data/1_processing/5.1_reCallPeak.outs/3_scATAC.annot.mitofilt.ndup.rds"
sc <- read_rds(fn)

MCls_name <- c("0"="0_CD4Naive", "1"="1_TCM", "2"="2_NKcell", "3"="3_TEM",
    "4"="4_Bcell", "5"="5_CD8Naive", "6"="6_Monocyte", "7"="7_dnT",
    "8"="8_MAIT", "9"="9_Platelets", "10"="10_DC")

x <- sc@meta.data
x <- x%>%mutate(MCls=MCls_name[as.character(wsnn_res.0.1)],
                treat2=ifelse(treats%in%c("etOH", "water"), "control", treats))
sc <- AddMetaData(sc, x)


### used for AnnotationPlot plot
FindRegion <- function(
   object,
   region,
   sep = c("-", "-"),
   assay = NULL,
   extend.upstream = 0,
   extend.downstream = 0
) {
  if (!is(object = region, class2 = "GRanges")) {
      # first try to convert to coordinates, if not lookup gene
      region <- tryCatch(
         expr = suppressWarnings(
           expr = StringToGRanges(regions = region, sep = sep)
           ),
           error = function(x) {
              region <- LookupGeneCoords(
                 object = object,
                 assay = assay,
                 gene = region
                 )
                 return(region)
            }
        )
        if (is.null(x = region)) {
           stop("Gene not found")
         }
      }
      region <- suppressWarnings(expr = Extend(
          x = region,
          upstream = extend.upstream,
          downstream = extend.downstream
       )
       )
       return(region)
 }


###
###
MCls_sel <- c("0_CD4Naive", "1_TCM", "2_NKcell", "3_TEM", "4_Bcell", "5_CD8Naive", "6_Monocyte")

col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
  "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
   "6_Monocyte"="#984ea3", "7_dnT"="black")
col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
       "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3",
       "control"="grey70")


sc2 <- subset(sc, MCls=="6_Monocyte"&treat2%in%c("zinc", "control"))
DefaultAssay(sc2) <- "ATAC"


genes <- c("MT1G", "MT1H", "MT1A", "MT1M", "MT1X", "MT1F", "MT1E")


### differential ATAC results 
#######################################

fn <- "./sc_multiome_data/2_Differential/1.2_DiffPeak.outs/2.0_DESeq.results_0.1_cn_sampleID+new_treat.rds"
res_DAR <- read_rds(fn)%>%
    dplyr::filter(MCls=="6_Monocyte", contrast=="zinc")%>%
    mutate(comb=paste(MCls, contrast, sep="_"),
           is_sig=ifelse(p.adjusted<0.1&abs(estimate)>0.5, 1, 0))
    ## filter(is_sig==1)
###
res_DAR <- res_DAR%>%dplyr::select(comb, peak=gene, baseMean, estimate, p.adjusted, is_sig)
 

### peak annotation
fn <- "./sc_multiome_data/2_Differential/2.2_compareRNAandATAC.outs/1_annot.signac_macs2_0.1_cn.rds"
peakAnno <- read_rds(fn)
peakAnno <- peakAnno%>%dplyr::select(peak=query_region, gene_id=gene_name, dtss=distance)%>%
    mutate(peak=gsub(":", "-", peak), dtss=abs(dtss))

res_DAR <- res_DAR%>%left_join(peakAnno, by="peak")

res2_DAR <- res_DAR%>%dplyr::filter(gene_id%in%genes, dtss<1e+05, p.adjusted<0.1) ##, is_sig==1) 

##
DAR <- unique(res2_DAR$peak) 




####
#### Find regions 

genes <- c("MT1G", "MT1H", "MT1A", "MT1M", "MT1X", "MT1F", "MT1E")
region_df <- map_dfr(genes, function(geneId){
###    
   region <- FindRegion(sc2, region=geneId, assay="ATAC", extend.upstream=0, extend.downstream=0)
   chr <- paste0("chr", as.character(seqnames(region)))
   s0 <- start(region)
   s1 <- end(region)
   ###
   df0 <- data.frame("chr"=chr, "s0"=s0, "s1"=s1)
   df0
})

chr <- region_df$chr[1]
s0 <- min(region_df$s0)
s1 <- max(region_df$s1)

region2 <- paste(chr, s0-0.5e+03, s1+3e+03, sep="-") ##"chr16-56666731-56668065"
 

###
### get peaks
peaks <- granges(sc2)
region0 <- FindRegion(sc2, region=region2, assay="ATAC", extend.upstream=0, extend.downstream=0)

peak_df <- subsetByOverlaps(x=peaks, ranges=region0)%>% 
    as.data.frame()%>%mutate(peak=paste(seqnames, start, end, sep="-")) 
peakSel <- peak_df$peak

peakSel2 <- intersect(DAR, peakSel)




##########################
### plots 
##########################


###
### coverage plots (1)
p0 <- CoveragePlot(sc2, region=region2, group.by="treat2", region.highlight=StringToGRanges(peakSel2),
                   extend.upstream=0e+03, extend.downstream=0e+03, annotation=F, peaks=F, tile=F, links=F,
                   max.downsample=3000,
                   )&
    scale_fill_manual(values=col2)&
   ## ylab("Normalized accessibility")&       
   ggtitle("MT family in Monocyte")&
   theme(plot.title=element_text(hjust=0.5, size=10),
          axis.title.y=element_text(size=10),
          axis.text.x=element_blank(),
          axis.ticks=element_blank()) 

## p0$layers[[5]]    

### Tile plots
## p1 <- TilePlot(sc2, region=region2, group.by="treat2", extend.upstream=5e+03, extend.downstream=5e+03)

### 
### gene annotation, (2)  
region2a <- gsub("chr", "", region2)  
p2 <- AnnotationPlot(sc2, region=region2a)& ##, extend.upstream=2e+03, extend.downstream=2e+03)&
    theme(## text=element_text(size=4),
          axis.title.y=element_text(size=10))
nn <- length(p2$layers)

p2$layers[[nn]]$aes_params$size <- 1.8
p2$layers[[nn]]$aes_params$fontface <- "italic"


###
### peak plots, (3)
p3 <- PeakPlot(sc2, region=region2)& ##, extend.upstream=2e+03, extend.downstream=2e+03)&
    theme(axis.text.x=element_blank(),
          axis.ticks=element_blank(),
          axis.title=element_text(size=10))

###
comb <- CombineTracks(plotlist=list(p0, p2, p3), heights=c(5, 2.5, 2))


### output
figfn <- paste(outdir2, "Figure3D_coverage.png", sep="")
ggsave(figfn, comb, width=700, height=550, units="px", dpi=120)






## fn <- "./sc_multiome_data/2_Differential/1.2_DiffPeak.outs/2.0_DESeq.results_0.1_cn_sampleID+new_treat.rds"
## res <- read_rds(fn)%>%as.data.frame()%>%
##     mutate(comb=paste(MCls, contrast, sep="_"),
##            is_sig=ifelse(p.adjusted<0.1&abs(estimate)>0.5, 1, 0))

## res2 <- res%>%dplyr::filter(comb=="6_Monocyte_zinc")%>%
##     dplyr::select(peak=gene, comb,  baseMean, estimate, zscore=statistic, p.adjusted, is_sig)

## peak_mat <- str_split(res2$peak, "-", simplify=T)
## colnames(peak_mat) <- c("chr", "start", "end")

## res2 <- cbind(res2, peak_mat)

## res3 <- res2%>%dplyr::filter(chr=="chr16", as.numeric(start)>(s0-2e+03), as.numeric(end)<(s1+2e+03))


##############################################
### box and violin plots of peaks and gene ###
##############################################

library(limma)
library(edgeR)
library(DESeq2)


###
### function 
getData <- function(gene, mat, cvt){
   ##
   ##gene <- "MA0098.3" 
   cvt$y <- as.numeric(mat[gene,])
   MCls <- sort(unique(cvt$MCls))

   ## control has replicates infividuals from etOH and water 
   cvt <- cvt%>%group_by(MCls, treats, ind)%>%summarise(y=mean(y, na.rm=T), .groups="drop")%>%ungroup()   

   ##
   cvt2 <- map_dfr(MCls, function(oneMCl){
      ##
      x <- cvt%>%dplyr::filter(MCls==oneMCl)
      x0 <- x%>%dplyr::filter(treats=="control")
      y0 <- x0$y
      names(y0) <- x0$ind
      ### 
      x <- x%>%mutate("y0"=y0[as.character(ind)])
      x
   })
   cvt2 <- cvt2%>%mutate(y2=y-y0)
   cvt2
}    


###
### setting colors
col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
  "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
   "6_Monocyte"="#984ea3", "7_dnT"="black")
col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
       "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3", "control"="grey70")


### genes and peakSel from the above script


###
### differential gene expression data
fn <- "./sc_multiome_data/2_Differential/1_DiffRNA_2.outs/1_YtX.sel_clusters_0.1_cn.rds"
data <- read_rds(fn)

### normalized  
dge <- DGEList(data)
dge <- calcNormFactors(dge)
v <- voom(dge)
mat2 <- v$E

### manual normalization 
## depth <- colSums(mat)
## mat2 <- sweep(mat, 2, depth, "/")
## mat2 <- log2(mat2*1e+06+1)


###
x <- str_split(colnames(mat2), "_", simplify=T)
cvt <- data.frame(MCls=paste(x[,1], x[,2], sep="_"), treats=x[,3], ind=x[,4])
cvt <- cvt%>%mutate(treats=ifelse(treats%in%c("etOH", "water"), "control", treats))


## genes <- c("MT1G", "MT1H", "MT1A", "MT1M", "MT1X", "MT1F", "MT1E")

MCls_sel <- "6_Monocyte"

## sort gene
gene2 <- c("MT1E", "MT1M", "MT1A", "MT1F", "MT1G", "MT1H", "MT1X")
plotDF <- map_dfr(1:length(gene2), function(i){
###
gene0 <- gene2[i]    
tmp <- getData(gene0, mat2, cvt)
tmp2 <- tmp%>%dplyr::filter(treats%in%c("control", "zinc"), MCls%in%MCls_sel)
tmp2$genes <- gene0
tmp2$gene_val <- i
tmp2    
})    

plotDF <- plotDF%>%mutate(gene2=fct_reorder(genes, gene_val))

###
###
## comb_sigs <- res%>%filter(is_sig==1, gene==gene0)%>%pull(comb) 

## sigs_df <- plotDF2%>%group_by(MCls, treats)%>%
##     summarize(ymax=max(y2), .groups="drop")%>%
##     mutate(MCls_value=as.numeric(MCls_val[MCls]),
##            MCl2=fct_reorder(MCls, MCls_value), comb=paste(MCls, treats, sep="_"),
##            sigs=ifelse(comb%in%comb_sigs, "Sig", "Not"))

   
p2 <- ggplot(plotDF, aes(x=factor(gene2), y=y, color=treats))+
   ## geom_violin(position="dodge", width=0.8)+
   geom_boxplot(outlier.shape=NA)+
   geom_point(position=position_jitterdodge(jitter.width=0.25, dodge.width=0.75),
              size=0.6)+ 
   ## geom_text(data=sigs_df, aes(x=MCl2, y=ymax+0.05, label=sigs), size=3)+
   scale_y_continuous("Normalized gene expression",
      expand=expansion(mult=c(0.2, 0.2)))+
   scale_color_manual("", values=col2)+
   ggtitle("Gene expression")+ 
   ## facet_wrap(~treats)+
   theme_bw()+
   theme(## axis.text.x=element_text(angle=-90, hjust=0, size=10),
         axis.text.x=element_text(size=9, face="italic"),
         axis.title=element_blank(),
         ## axis.title.y=element_text(size=10),
         axis.text.y=element_text(size=9),
         plot.title=element_text(size=10, hjust=0.5),
          legend.position="none")
         
         ###strip.text=element_text(size=10), 
         ##plot.title=element_text(hjust=0.5, size=12),
         ###legend.position="none")
  
figfn <- paste(outdir2, "Figure3D_gene_box.png", sep="")
ggsave(figfn, p2, width=480, height=300, units="px", dpi=120)






###
### peak data
   
fn <- "./sc_multiome_data/2_Differential/1.2_DiffPeak.outs/1_YtX.sel_0.1_cn.rds"
data <- read_rds(fn)

### normalized  
dge <- DGEList(data)
dge <- calcNormFactors(dge)
v <- voom(dge)
mat2 <- v$E

## depth <- colSums(mat)
## mat2 <- sweep(mat, 2, depth, "/")
## mat2 <- log2(mat2*1e+06+1)

###
x <- str_split(colnames(mat2), "_", simplify=T)
cvt <- data.frame(MCls=paste(x[,1], x[,2], sep="_"), treats=x[,3], ind=x[,4])
cvt <- cvt%>%mutate(treats=ifelse(treats%in%c("etOH", "water"), "control", treats))


###
### plot data
MCls_sel <- "6_Monocyte"
peakSel2 <- sort(peakSel2)
plotDF2 <- map_dfr(1:length(peakSel2), function(i){
###
peak0_id <- peakSel2[i]
peak0 <- paste("peak", i, sep="")    
tmp <- getData(peak0_id, mat2, cvt)
tmp2 <- tmp%>%dplyr::filter(treats%in%c("control", "zinc"), MCls%in%MCls_sel)
tmp2$peaks_id <- peak0_id
tmp2$peaks <- peak0    
tmp2
})



###
p3 <- ggplot(plotDF2, aes(x=factor(peaks), y=y, color=treats))+
   ## geom_violin(position="dodge", width=0.8)+
   geom_boxplot(outlier.shape=NA)+
   geom_point(position=position_jitterdodge(jitter.width=0.25, dodge.width=0.75),
              size=0.6)+ 
   ## geom_text(data=sigs_df, aes(x=MCl2, y=ymax+0.05, label=sigs), size=3)+
   scale_y_continuous("Normalized accessibility",
      expand=expansion(mult=c(0.2, 0.2)))+
   ## scale_x_discrete(
   ##     labels=c("Peak1"=bquote("Peak ("~italic(MT1G)~")"), "Peak2"=bquote("Peak ("~italic(MT1H)~")")))+ 
   scale_color_manual("", values=col2)+
   ggtitle("Chromatin accessibility")+ 
   ## facet_wrap(~treats)+
   theme_bw()+
   theme(## axis.text.x=element_text(angle=-90, hjust=0, size=10),
         axis.text=element_text(size=9),
         axis.title=element_blank(),
         plot.title=element_text(hjust=0.5, size=10),
         ## axis.title.y=element_text(size=10),
         ## axis.text.y=element_text(size=10),
          legend.position="none")
         
         ###strip.text=element_text(size=10), 
         ##plot.title=element_text(hjust=0.5, size=12),
         ###legend.position="none")
  
figfn <- paste(outdir2, "Figure3D_peak_box.png", sep="")
ggsave(figfn, p3,  width=480, height=300, units="px", dpi=120)




###
### DEGs
fn <- "./sc_multiome_data/2_Differential/1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn_sampleID+new_treat.rds"
res <- read_rds(fn)%>%as.data.frame()%>%
    mutate(comb=paste(MCls, contrast, sep="_"),
           is_sig=ifelse(p.adjusted<0.1&abs(estimate)>0.5, 1, 0))

###
res3 <- res%>%filter(MCls=="6_Monocyte", contrast=="zinc", grepl("^MT", gene))%>%arrange(desc(abs(estimate)))%>%
   dplyr::select(comb, gene, baseMean, estimate, stderror, zscore=statistic, p.adjusted) 


fn <- "./3_compare_plots.outs/1.2_LFC_gene_ATAC_comb.rds"
x <- read_rds(fn)
x2 <- x%>%filter(comb=="6_Monocyte_zinc")%>%filter(grepl("^MT", gene))
    
###
### DARs
 
fn <- "./sc_multiome_data/2_Differential/1.2_DiffPeak.outs/2.0_DESeq.results_0.1_cn_sampleID+new_treat.rds"
res_DAR <- read_rds(fn)%>%as.data.frame()%>%
    mutate(comb=paste(MCls, contrast, sep="_"),
           is_sig=ifelse(p.adjusted<0.1&abs(estimate)>0.5, 1, 0))

res2_DAR <- res_DAR%>%dplyr::filter(comb=="6_Monocyte_zinc")%>%
    dplyr::select(peak=gene, comb,  baseMean, estimate, zscore=statistic, p.adjusted, is_sig)


### peak annotation
fn <- "./sc_multiome_data/2_Differential/2.2_compareRNAandATAC.outs/1_annot.signac_macs2_0.1_cn.rds"
peakAnno <- read_rds(fn)
peakAnno <- peakAnno%>%dplyr::select(peak=query_region, gene_id=gene_name, dtss=distance)%>%
    mutate(peak=gsub(":", "-", peak), dtss=abs(dtss))

res2_DAR <- res2_DAR%>%left_join(peakAnno, by="peak")

###
res3_DAR <- res2_DAR%>%dplyr::filter(gene_id%in%genes, is_sig>0)

peakSel <- sort(unique(res3_DAR$peak))
