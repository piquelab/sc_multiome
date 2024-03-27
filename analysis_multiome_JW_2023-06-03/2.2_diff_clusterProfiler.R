##
library(Matrix)
library(tidyverse)
library(data.table)
##
library(annotables)
##
library(aplot, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
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

rm(list=ls())

outdir <- "./2_diff_plots.outs/enriched_outs/"
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)



###
###
### DE results

###
### Differential results
fn <- "./sc_multiome_data/2_Differential/1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn_sampleID+new_treat.rds"
res <- read_rds(fn)%>%as.data.frame()%>%mutate(comb=paste(MCls, contrast, sep="_"))
res2 <- res%>%filter(p.adjusted<0.1, abs(estimate)>0.5)%>%
   mutate(direction=ifelse(estimate>0, "Up", "Down")) # sign(estimate)))


### Background genes
BG <- unique(res$gene)
geneBG <- bitr(BG, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)


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
opfn <- paste(outdir, "1_enrichGO.rds", sep="")
write_rds(cg, opfn)


###
### KEGG enrichment analysis
## ck <- compareCluster(ENTREZID~contrast+MCls+direction,
##     data=geneCluster,
##     universe=geneBG$ENTREZID,
##     organism="hsa",
##     fun="enrichKEGG",
##     pvalueCutoff=1,
##     qvalueCutoff=1,
##     minGSSize=0,
##     maxGSSize=1000)
  
## ###
## ###
## opfn <- paste(outdir, "2_enrichKEGG.rds", sep="")
## write_rds(ck, opfn)


###############################
### visulization GO results ###
###############################

fn <- paste(outdir, "1_enrichGO.rds", sep="")
cg <- read_rds(fn)


### 
cg <- cg%>%mutate(Cluster2=paste(MCls, contrast, sep="_"),
   maxGSSize=as.numeric(gsub("/.*", "", BgRatio)),
   ngene=as.numeric(gsub(".*/", "", GeneRatio)) )

x <- cg@compareClusterResult%>%dplyr::select(Cluster2, MCls, contrast)%>%distinct(Cluster2, .keep_all=T)%>%
    arrange(contrast)
Cl_value <- 1:nrow(x)
names(Cl_value) <- x$Cluster2

cg <- cg%>%mutate(Cluster_val=Cl_value[Cluster2],
                  ClusterNew=fct_reorder(Cluster2, Cluster_val))


###
###
cg2 <- cg%>%dplyr::filter(direction=="Up", maxGSSize>10, maxGSSize<500, ngene>5, p.adjust<0.05) 
p1 <- enrichplot::dotplot(cg2, x="ClusterNew", showCategory=5)+
   scale_y_discrete(labels=function(y) str_wrap(y, width=100))+ 
   theme(axis.title=element_blank(),
         axis.text.x=element_markdown(angle=45, hjust=1,size=14),
         axis.text.y=element_text(size=12),
         legend.text=element_text(size=12),
         legend.title=element_text(size=12))

 
###
figfn <- paste(outdir, "Figure1.1_up_GO_enriched.pdf", sep="")
ggsave(figfn, p1, device="pdf", width=16, height=14)



###
### Down-regulated

cg2 <- cg%>%dplyr::filter(direction=="Down", maxGSSize>10, maxGSSize<500, ngene>5, p.adjust<0.05) 
p2 <- enrichplot::dotplot(cg2, x="ClusterNew", showCategory=5)+
   scale_y_discrete(labels=function(y) str_wrap(y, width=80))+ 
   theme(axis.title=element_blank(),
         axis.text.x=element_markdown(angle=45, hjust=1,size=14),
         axis.text.y=element_text(size=12),
         legend.text=element_text(size=12),
         legend.title=element_text(size=12))
###
figfn <- paste(outdir, "Figure1.2_down_GO_enriched.pdf", sep="")
ggsave(figfn, p2, device="pdf", width=18, height=20)




#############################
### visulization of plots ###
#############################

rm(list=ls())

outdir2 <- "./2_diff_plots.outs/enriched_trees/"
if ( !file.exists(outdir2) ) dir.create(outdir2, showWarnings=F, recursive=T)


####
### differential results 
fn <- "./sc_multiome_data/2_Differential/1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn_sampleID+new_treat.rds"
res <- read_rds(fn)%>%as.data.frame()%>%mutate(comb=paste(MCls, contrast, sep="_"))
res2 <- res%>%filter(p.adjusted<0.1, abs(estimate)>0.5)%>%
    mutate(direction=ifelse(estimate>0, "Up", "Down"))

###
### Background genes
BG <- unique(res$gene)
geneBG <- bitr(BG, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)



geneCluster <- res2%>%
    dplyr::select(contrast, MCls, comb, direction,  SYMBOL=gene)%>%
    inner_join(geneBG, by="SYMBOL")


###
### GO enrichment and tree plots for each combination 



geneCluster2 <- geneCluster%>%filter(direction=="Up")
comb <- sort(unique(geneCluster2$comb))
for (ii in comb){
##
    gene_test <- geneCluster2%>%filter(comb==ii)%>%pull(ENTREZID)
    cat(ii, "\n")
    ###
    cg <- enrichGO(gene_test,
        universe=geneBG$ENTREZID, 
        OrgDb="org.Hs.eg.db",
        pvalueCutoff=1, qvalueCutoff=1,
        ont="ALL",
        minGSSize=5,
        maxGSSize=500)
    
    cg2 <- cg%>%dplyr::filter(p.adjust<0.1)

    
    cg2 <- try(enrichplot::pairwise_termsim(cg2), silent=T)
    if (class(cg2)=="try-error") next
    
    p0 <- try(enrichplot::treeplot(cg2, color="p.adjust", nWords=0, cex_category=0.8), silent=T)
    if (class(p0)=="try-error") next
    
     p0 <- p0+ 
        ggtitle(ii)+
        theme(plot.title=element_text(hjust=0.5, size=12),
              ##axis.text=element_text(size=4),
              legend.title=element_text(size=8),
              legend.text=element_text(size=8),
              legend.key.size=unit(0.2, "cm"))
 
    
figfn <- paste(outdir2, "Figure1_", ii, "_up_tree.png", sep="")
ggsave(figfn, p0, width=900, height=580, units="px", dpi=100)

}


###
### down
geneCluster2 <- geneCluster%>%filter(direction=="Down")
comb <- sort(unique(geneCluster2$comb))
for (ii in comb){
##
    gene_test <- geneCluster2%>%filter(comb==ii)%>%pull(ENTREZID)
    cat(ii, "\n")
    ###
    cg <- enrichGO(gene_test,
        universe=geneBG$ENTREZID, 
        OrgDb="org.Hs.eg.db",
        pvalueCutoff=1, qvalueCutoff=1,
        ont="ALL",
        minGSSize=5,
        maxGSSize=500)
    
    cg2 <- cg%>%dplyr::filter(p.adjust<0.1)
    
    cg2 <- try(enrichplot::pairwise_termsim(cg2), silent=T)
    if ( class(cg2)=="try-error") next
    
    p0 <- try(enrichplot::treeplot(cg2, color="p.adjust", nWords=0, cex_category=0.8), silent=T)
    
    if (class(p0)=="try-error") next
    
     p0 <- p0+ 
        ggtitle(ii)+
        theme(plot.title=element_text(hjust=0.5, size=12),
              ##axis.text=element_text(size=4),
              legend.title=element_text(size=8),
              legend.text=element_text(size=8),
              legend.key.size=unit(0.2, "cm"))
 
    
figfn <- paste(outdir2, "Figure2_", ii, "_down_tree.png", sep="")
ggsave(figfn, p0, width=900, height=580, units="px", dpi=100)

}
 

##########################################
### visulization of tree plots version ###
##########################################

rm(list=ls())

outdir <- "./2_diff_plots.outs/enriched_outs/"
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)

library(dendextend)


###
### calculate odds ratio
fn <- paste(outdir, "1.2_enrichGO_odds.rds", sep="")
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
fn <- paste(outdir, "1.2_enrichGO_odds.rds", sep="")
write_rds(cg, file=fn)


##########################
### enrichment results ###
##########################

fn <- paste(outdir, "1.2_enrichGO_odds.rds", sep="")
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
row_dend <- as.dendrogram(hclust(dsim2))


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



df_col <- data.frame(celltype=cvt2$MCls, treats=cvt2$treats)
col_ha <- HeatmapAnnotation(df=df_col, col=list(celltype=col1, treats=col2),
    annotation_name_gp=gpar(fontsize=10),
    annotation_legend_param=list(
  celltype=list(title_gp=gpar(fontsize=10), labels_gp=gpar(fontsize=10), title="celltype",
                grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm")),
  treats=list(title_gp=gpar(fontsize=10), labels_gp=gpar(fontsize=10), title="treats",
                grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm"))))
  ## show_legend=c(F,F))
  ##simple_anno_size=unit(0.3, "cm"))



 
mat3 <- mat2
mat3[not_sig] <- NA
 
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
figfn <- paste(outdir, "Figure2.1_enriched_up.heatmap.png", sep="")
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



df_col <- data.frame(celltype=cvt2$MCls, treats=cvt2$treats)
col_ha <- HeatmapAnnotation(df=df_col, col=list(celltype=col1, treats=col2),
    annotation_name_gp=gpar(fontsize=10),
    annotation_legend_param=list(
  celltype=list(title_gp=gpar(fontsize=10), labels_gp=gpar(fontsize=10), title="celltype",
                grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm")),
  treats=list(title_gp=gpar(fontsize=10), labels_gp=gpar(fontsize=10), title="treats",
                grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm"))))
  ##show_legend=c(F,F),
  ##simple_anno_size=unit(0.3, "cm"))




### plot
mat3 <- mat2
mat3[not_sig] <- NA

rownames(mat3)[26] <- "adaptive immune response based on somatic recombination of immune receptors"

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
figfn <- paste(outdir, "Figure2.2_enriched_down.heatmap.png", sep="")
## ggsave(figfn, plot=p2, device="png", width=680, height=700, units="px", dpi=300) ## ggsave 
png(figfn, width=1850, height=1300, res=120)
set.seed(0)
p2 <- draw(p2, heatmap_legend_side="left", padding=unit(c(0.2, 0.2, 0.2, 6), "cm"),
           annotation_legend_side="left")
dev.off()




###################################################
### treeplots for comparison of clusterprofiler ###
###################################################


## oneMCl <- MCls[7]
## geneCluster3 <- geneCluster2%>%filter(MCls==oneMCl)
## cg <- compareCluster(ENTREZID~contrast,
##     data=geneCluster3,
##     universe=geneBG$ENTREZID,
##     fun="enrichGO",
##     OrgDb="org.Hs.eg.db",
##     pvalueCutoff=1,
##     qvalueCutoff=1,
##     ont="ALL",
##     minGSSize=5,
##     maxGSSize=500)

## cg2 <- enrichplot::pairwise_termsim(cg)
## p0 <- enrichplot::treeplot(cg2, color="p.adjust") ##, group_colors=col2)+
##    theme(plot.title=element_text(hjust=0.5, size=12),
##          legend.title=element_text(size=8),
##          legend.text=element_text(size=8),
##          legend.key.size=unit(0.2, "cm"))

## figfn <- paste(outdir, "Figure3.1_", oneMCl, "_up_tree.png", sep="")
## ggsave(figfn, pcomb, width=800, height=400, units="px", dpi=120)
         




## cluster.params=list(cluster="Cluster", method=cluster::pam)
## p2 <- emapplot(cg2, showCategory=100, pie="count", cex_category=1, group_category=T) #, cluster.params=cluster.params)

## figfn <- paste(outdir, "Figure2.1_emap.png", sep="")
## ggsave(figfn, p2, device="png", width=620, height=620, units="px", dpi=120)


###
### tree plots

### install packages
## library(devtools)
## library(aplot, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
## install_github("YuLab-SMU/enrichplot", lib="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")

