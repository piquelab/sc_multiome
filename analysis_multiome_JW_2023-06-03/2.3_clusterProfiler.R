##
library(Matrix)
library(tidyverse)
library(data.table)
##
library(annotables)
##
library(aplot) ##, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(igraph)

library(clusterProfiler) ##, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
## library(enrichplot, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(GOSemSim)
library(org.Hs.eg.db)
##
library(cowplot)
library(RColorBrewer)
library(ComplexHeatmap)
library(viridis)
library(circlize)
library(ggrepel)
library(ggrastr)
library(openxlsx)

library(ggtext)
library(glue)
library(dendextend)


rm(list=ls())

outdir <- "./Plots_pub/2_diff_plots//"
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)
 
 

###
### Differential results
fn <- "./sc_multiome_data/2_Differential/1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn_sampleID+new_treat.rds"
res <- read_rds(fn)%>%as.data.frame()%>%mutate(comb=paste(MCls, contrast, sep="_"))

gene0 <- unique(res$gene)
gene1 <- grch38%>%dplyr::filter(chr%in%as.character(1:22))%>%pull(symbol)%>%unique()
gene2 <- intersect(gene0, gene1)


res2 <- res%>%filter(p.adjusted<0.1, abs(estimate)>0.5, gene%in%gene2)%>%
   mutate(direction=ifelse(estimate>0, "Up", "Down")) # sign(estimate)))


### Background genes
BG <- gene2
geneBG <- bitr(BG, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)%>%
    distinct(SYMBOL, .keep_all=T)

 
geneCluster <- res2%>%dplyr::select(contrast, MCls, direction, SYMBOL=gene)%>%
    inner_join(geneBG, by="SYMBOL")
 

###
### GO enrichment
cg <- compareCluster(ENTREZID~contrast+MCls+direction,
    data=geneCluster,
    universe=geneBG$ENTREZID,
    fun="enrichGO",
    OrgDb="org.Hs.eg.db",
    pvalueCutoff=1,
    qvalueCutoff=1,
    ont="ALL",
    minGSSize=0,
    maxGSSize=1000)
 
###
###
opfn <- paste(outdir, "1.2_DEG_enrichGO.rds", sep="")
write_rds(cg, opfn)




##########################################
### visulization of tree plots version ###
##########################################

rm(list=ls())

outdir <- "./Plots_pub/2_diff_plots/"
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)



###
### calculate odds ratio
fn <- paste(outdir, "1.2_DEG_enrichGO.rds", sep="")
cg <- read_rds(fn)

x <- cg@compareClusterResult
x2 <- x%>%dplyr::select(Cluster, ID, GeneRatio, BgRatio)%>%
   mutate(n1=as.numeric(gsub("/.*", "", GeneRatio)),
          ndeg=as.numeric(gsub(".*/", "", GeneRatio)),
          n2=as.numeric(gsub("/.*", "", BgRatio)),
          nbg=as.numeric(gsub(".*/", "", BgRatio)))

###
### odds ratio
df2 <- map_dfr(1:nrow(x2), function(i){
   ###
   n1 <- x2$n1[i]
   ndeg <- x2$ndeg[i]
   n2 <- x2$n2[i]
   nbg <- x2$nbg[i]
    
   ng_1 <- n1  ## interest.in
   ng_2 <- ndeg-n1 ## interest.not 

   ng_3 <- n2-n1   ## not.interest.in
   ng_4 <- nbg-ndeg-ng_3  ## not.interest.not

   ##
   dmat <- matrix(c(ng_1, ng_2, ng_3, ng_4), 2, 2)
   colnames(dmat) <- c("interest.in", "not.interest.in")
   rownames(dmat) <- c("interest.not", "not.interest.not")
   ##
   enrich <- fisher.test(dmat)
   df0 <- data.frame("Cluster"=x2$Cluster[i], ID=x2$ID[i], log_odds=enrich$estimate) 
   df0
})    

df2$ID <- x2$ID
df2$log_odds <- log(df2$log_odds)

df3 <- df2%>%mutate(rn=paste(Cluster, ID, sep="_"))%>%dplyr::select(rn, log_odds)
xnew <- x%>%mutate(rn=paste(Cluster, ID, sep="_"))%>%left_join(df3, by="rn")

### output
cg@compareClusterResult <- xnew
fn <- paste(outdir, "1.2_DEG_enrichGO_odds.rds", sep="")
write_rds(cg, file=fn)


#################################################
### compare current vs old enrichment results ###
#################################################


## fn <- "./Plots_pub/2_diff_plots/1.2_DEG_enrichGO_odds.rds"
## cg <- read_rds(fn)
## xnew <- cg@compareClusterResult

## xnew2 <- xnew%>%dplyr::select(rn, MCls, contrast, direction, ID, pval_new=pvalue)

## ###
## fn <- "./2_diff_plots.outs/old/1.2_enrichGO_odds.rds"
## old <- read_rds(fn)@compareClusterResult
## old2 <- old%>%dplyr::select(rn, pval_old=pvalue) 


## res <- xnew2%>%inner_join(old2, by="rn")
## res <- res%>%mutate(log10p_new=-log10(pval_new), log10p_old=-log10(pval_old))

## col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
##        "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3")

 
## ###
## ### Up 
## ## plotDF <- res%>%filter(direction=="Up")
## p0 <- ggplot(res, aes(x=log10p_old, y=log10p_new, color=contrast))+
##    geom_point(size=0.4)+
##    geom_abline(slope=1, intercept=0, color="grey")+ 
##    facet_grid(contrast~MCls, scales="free")+
##    scale_color_manual(values=col2, guide=guide_legend(override.aes=list(size=2)))+
##    ggtitle("Enrichment comparison")+
##    xlab(bquote(~-log[10]~italic(p)~"Not filter MT and scaffold genes"))+
##    ylab(bquote(~-log[10]~italic(p)~"After filter"))+    
##    theme_bw()+
##    theme(plot.title=element_text(hjust=0.5, size=12),
##          legend.title=element_blank(),
##          legend.key=element_blank(),
##          legend.key.size=grid::unit(0.4, "cm"),
##          legend.background=element_blank(),
##          strip.text.x=element_text(size=10))

## figfn <- "./2_diff_plots.outs/old/Figure3_compare_enrich.png"
## ggsave(figfn, p0, width=1200, height=700, units="px", dpi=120)

   




##########################
### enrichment results ###
##########################

rm(list=ls())

outdir <- "./Plots_pub/2_diff_plots/"
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)


fn <- paste(outdir, "1.2_DEG_enrichGO_odds.rds", sep="")
cg <- read_rds(fn)

cg <- cg%>%mutate(Cluster2=paste(MCls, contrast, sep="_"),
   maxGSSize=as.numeric(gsub("/.*", "", BgRatio)),
   ngene=as.numeric(gsub(".*/", "", GeneRatio)) )

## x <- cg@compareClusterResult%>%dplyr::select(Cluster2, MCls, contrast)%>%distinct(Cluster2, .keep_all=T)%>%
##     arrange(contrast)
## Cl_value <- 1:nrow(x)
## names(Cl_value) <- x$Cluster2

## cg <- cg%>%mutate(Cluster_val=Cl_value[Cluster2],
##                   ClusterNew=fct_reorder(Cluster2, Cluster_val))
 

###
### tree plots for Up-regulated DEGs

cg2 <- cg%>%dplyr::filter(direction=="Up", maxGSSize>10, maxGSSize<500, ngene>5, p.adjust<0.05) 

x <- cg2@compareClusterResult

GOsel <- x%>%filter(ONTOLOGY=="BP")%>%
    group_by(Cluster2)%>%slice_min(order_by=pvalue, n=5)%>%pull(Description)%>%unique()


###
### similarity and cluster analysis based on distance
d <- GOSemSim::godata("org.Hs.eg.db", ont="BP")
cg2 <- enrichplot::pairwise_termsim(cg2, semData=d, method="Wang")
 
sim <- cg2@termsim
sim2 <- sim[GOsel, GOsel]
dsim2 <- as.dist(1-sim2)
row_dend <- as.dendrogram(hclust(dsim2)) ##, method="median"))

 
###
### Heatmap

###
### generate matrix for Heatmap plots
dfmat <- x%>%pivot_wider(id_cols=Description, names_from=Cluster2, values_from=log_odds, values_fill=0)
mat <- dfmat%>%column_to_rownames(var="Description")%>%as.matrix()

cvt <- str_split(colnames(mat), "_", simplify=T)
cvt2 <- data.frame(cluster=colnames(mat), MCls=paste(cvt[,1], cvt[,2], sep="_"), treats=cvt[,3])
cvt2 <- cvt2%>%arrange(treats)

mat2 <- mat[GOsel,cvt2$cluster]

### significance
not_sig <- x%>%
    pivot_wider(id_cols=Description, names_from=Cluster2, values_from=p.adjust, values_fill=1)%>%
    column_to_rownames(var="Description")%>%as.matrix()
not_sig <- not_sig[GOsel, cvt2$cluster]
not_sig <- not_sig>0.05


 
### colors
y <- as.numeric(mat2)
quantile(y, probs=c(0.01, 0.1, 0.9, 0.95, 0.98, 0.99, 0.999, 1))

mybreak <- seq(0, 6, length.out=50)
 
mycol <- colorRamp2(mybreak, colorRampPalette(brewer.pal(n=9,name="Reds"))(50))

 
###
### annotation columns
col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
  "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
   "6_Monocyte"="#984ea3", "7_dnT"="black")
col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
       "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3")



df_col <- data.frame(celltype=cvt2$MCls, treatment=cvt2$treats)
col_ha <- HeatmapAnnotation(df=df_col, col=list(celltype=col1, treatment=col2),
    annotation_name_gp=gpar(fontsize=10),
    annotation_legend_param=list(
  celltype=list(title_gp=gpar(fontsize=10), labels_gp=gpar(fontsize=10), title="celltype",
                grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm")),
  treatment=list(title_gp=gpar(fontsize=10), labels_gp=gpar(fontsize=10), title="treatment",
                grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm"))))
  ## show_legend=c(F,F))
  ##simple_anno_size=unit(0.3, "cm"))



 
mat3 <- mat2
mat3[not_sig] <- NA

rownames(mat3)[2] <- "transmembrane receptor protein tyrosine kinase"
rownames(mat3)[4] <- "transmembrane receptor protein serine/threonine kinase"
rownames(mat3)[58] <- "immune response-activating cell surface receptor"
rownames(mat3)[59] <- "immune response-regulating cell surface receptor"


p1 <- Heatmap(mat3, col=mycol, na_col="white",
   rect_gp=gpar(col="black", lwd=1),           
   cluster_rows=row_dend, cluster_columns=F,
   row_dend_width=grid::unit(2, "cm"), 
   show_row_names=T, row_names_gp=gpar(fontsize=10),
   row_names_max_width=unit(8, "cm"),
   show_column_names=T, column_names_gp=gpar(fontsize=10),
   column_names_rot=-45,   
   top_annotation=col_ha,
   heatmap_legend_param=list(title="log odds",
      title_gp=gpar(fontsize=10),
      at=seq(0, 6, 2), 
      labels_gp=gpar(fontsize=10),
      grid_width=grid::unit(0.5, "cm"),
      legend_height=grid::unit(7.8, "cm")))
 
###
figfn <- paste(outdir, "FigS2_4_enriched_up.heatmap.png", sep="")
## ggsave(figfn, plot=p2, device="png", width=680, height=700, units="px", dpi=300) ## ggsave 
png(figfn, width=1500, height=1100,res=120)
set.seed(0)
p1 <- draw(p1, heatmap_legend_side="left", padding=unit(c(0.2, 0.2, 0.2, 3.2), "cm"),
           annotation_legend_side="left")
dev.off()
 





############
### Down ###
############

###
### selected GO terms
cg2 <- cg%>%dplyr::filter(direction=="Down", maxGSSize>10, maxGSSize<500, ngene>5, p.adjust<0.05) 

x <- cg2@compareClusterResult

GOsel <- x%>%filter(ONTOLOGY=="BP")%>%
    group_by(Cluster2)%>%slice_min(order_by=pvalue, n=5)%>%pull(Description)%>%unique()


### similarity 
d <- GOSemSim::godata("org.Hs.eg.db", ont="BP")
cg2 <- enrichplot::pairwise_termsim(cg2, semData=d, method="Wang")

sim <- cg2@termsim
sim2 <- sim[GOsel, GOsel]
dsim2 <- as.dist(1-sim2)
row_dend <- as.dendrogram(hclust(dsim2))



###
### generate matrix for Heatmap plots
dfmat <- x%>%pivot_wider(id_cols=Description, names_from=Cluster2, values_from=log_odds, values_fill=0)
mat <- dfmat%>%column_to_rownames(var="Description")%>%as.matrix()

cvt <- str_split(colnames(mat), "_", simplify=T)
cvt2 <- data.frame(cluster=colnames(mat), MCls=paste(cvt[,1], cvt[,2], sep="_"), treats=cvt[,3])
cvt2 <- cvt2%>%arrange(treats)

mat2 <- mat[GOsel,cvt2$cluster]

### significance
not_sig <- x%>%
    pivot_wider(id_cols=Description, names_from=Cluster2, values_from=p.adjust, values_fill=1)%>%
    column_to_rownames(var="Description")%>%as.matrix()
not_sig <- not_sig[GOsel, cvt2$cluster]
not_sig <- not_sig>0.05


 

###
### setting colors
y <- as.numeric(mat2)
quantile(y, probs=c(0.01, 0.1, 0.9, 0.95, 0.98, 0.99, 0.999, 1))

mybreak <- seq(0, 6, length.out=100)
 
mycol <- colorRamp2(mybreak, colorRampPalette(brewer.pal(n=9,name="Reds"))(100))


 
###
### annotation columns
col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
  "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
   "6_Monocyte"="#984ea3", "7_dnT"="black")
col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
       "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3")



df_col <- data.frame(celltype=cvt2$MCls, treatment=cvt2$treats)
col_ha <- HeatmapAnnotation(df=df_col, col=list(celltype=col1, treatment=col2),
    annotation_name_gp=gpar(fontsize=10),
    annotation_legend_param=list(
  celltype=list(title_gp=gpar(fontsize=10), labels_gp=gpar(fontsize=10), title="celltype",
                grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm")),
  treatment=list(title_gp=gpar(fontsize=10), labels_gp=gpar(fontsize=10), title="treatment",
                grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm"))))
  ##show_legend=c(F,F),
  ##simple_anno_size=unit(0.3, "cm"))




### plot
mat3 <- mat2
mat3[not_sig] <- NA

nn <- sapply(rownames(mat3), nchar) 

rownames(mat3)[15] <- "adaptive immune response based on somatic recombination of immune receptors"
rownames(mat3)[16] <- "antigen processing and presentation of exogenous petide antigen"
rownames(mat3)[26] <- "CD4+ alpha-beta T differentiation in immune response"


p2 <- Heatmap(mat3, col=mycol, na_col="white",
   rect_gp=gpar(col="black", lwd=1),           
   cluster_rows=row_dend, cluster_columns=F,
   row_dend_width=grid::unit(2, "cm"), 
   show_row_names=T, row_names_gp=gpar(fontsize=10),
   row_names_max_width=unit(9, "cm"),
   show_column_names=T, column_names_gp=gpar(fontsize=9),
   column_names_rot=-45,   
   top_annotation=col_ha,
   heatmap_legend_param=list(title="log odds",
      title_gp=gpar(fontsize=10),
      at=seq(0, 6, 2), 
      labels_gp=gpar(fontsize=10),
      grid_width=grid::unit(0.5, "cm"),
      legend_height=grid::unit(7.8, "cm")))



###
figfn <- paste(outdir, "FigS2_5_enriched_down.heatmap.png", sep="")
## ggsave(figfn, plot=p2, device="png", width=680, height=700, units="px", dpi=300) ## ggsave 
png(figfn, width=1850, height=1400, res=120)
set.seed(0)
p2 <- draw(p2, heatmap_legend_side="left", padding=unit(c(0.2, 0.2, 0.2, 6), "cm"),
           annotation_legend_side="left")
dev.off()


