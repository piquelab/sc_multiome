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


##
library(cowplot)
library(RColorBrewer)
library(ComplexHeatmap)
library(viridis)
library(circlize)
library(ggrepel)
library(ggrastr)
library(openxlsx)


### Library CCA package
library("rainbow") ##, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages")
library("fds") ###, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages")
library(fda) ### lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages")
library(CCA) ### lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages")
library(expm)

library("GGally") 


rm(list=ls())

outdir2 <- "./4_motif_plots.outs/2_all_response/"
if ( !file.exists(outdir2) ) dir.create(outdir2, showWarnings=F, recursive=T)




###
### compare response changes in TF motif, TF encoding genes and TF-regulated genes 
### By JW, Nov-21-2023
### Last motified, Jan 29-2024


### File 1.0_cl6, and 1.2_*, 111 top 10 response motifs with cluster, use for old analysis
### File 3.0_TF, response motifs with gene name
### File 3.2_corr, correlation of 187 response motif between response changes in TF activity and TF gene expression
### File 3.2_corr_reorder.xlsx, 
### File 3_TF_activity_Gene.rds, response change in TF activity and TF gene expression



#################################################################################
### Response change in TF activity and TF gene expression
###################################################################################


####
#### union of response motifs 


rm(list=ls())
outdir2 <- "./4_motif_plots.outs/2_all_response/"


###
### motifs-genes
fn <- "./4_motif_plots.outs/response_motif_correct/2.2_response_motif_th0.1.txt"
DF <- read.table(fn, header=T)
DF <- DF%>%distinct(motif_name, .keep_all=T)%>%dplyr::select(motif_ID, motif_name)

###
nmotif <- nrow(DF)
DF2 <- map_dfr(1:nmotif, function(i){
   ###
   gene0 <- DF$motif_name[i]
   if ( grepl("::", gene0)){
      ##
      gene2 <- unique(unlist(str_split(gene0, "::")))
   }else{
      ### 
      gene2 <- unique(unlist(str_split(gene0, "-")))
   }
   ###
   df2 <- data.frame(motif_name=gene0, gene=gene2)
   df2
})

opfn <- paste(outdir2, "3.0_TF_gene.xlsx", sep="")
write.xlsx(DF2, file=opfn)

###
fn <- paste(outdir2, "3.0_TF_gene.xlsx", sep="")
DF2 <- read.xlsx(fn)
geneSel <- unique(DF2$gene)
motifSel <- unique(DF2$motif_name)




###
### response motif activity
fn <- "./4_motif_plots.outs/3.2_motif_diff_fitInd_JW.rds"
res_motif <- read_rds(fn)%>%
    mutate(comb=paste(MCls, contrast, sep="_"), zscore=beta/stderr,
                      is_sig=ifelse(qval<0.1, 1, 0))
res_motif <- res_motif%>%dplyr::select(comb, MCls, contrast, motif_name, zscore_motif=zscore, beta_motif=beta)
res_motif2 <- res_motif%>%filter(motif_name%in%motifSel)



#### gene expression
### Differential results
fn <- "./sc_multiome_data/2_Differential/1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn_sampleID+new_treat.rds"
res_gene <- read_rds(fn)%>%
    mutate(comb=paste(MCls, contrast, sep="_"),
           is_sig=ifelse(p.adjusted<0.1, 1, 0), zscore=statistic)
###
res_gene2 <- res_gene%>%filter(gene%in%geneSel)%>%
    dplyr::select(comb, gene, zscore_rna=zscore, beta_rna=estimate)%>%as.data.frame()


###
### combine diff gene and diff motif
comb2 <- sort(unique(res_motif$comb))
res <- map_dfr(comb2, function(ii){
    ### diff gene results
    res0 <- res_gene2%>%filter(comb==ii)%>%dplyr::select(-comb)
    res0 <- DF2%>%left_join(res0, by="gene")%>%drop_na(zscore_rna, beta_rna)
    res0 <- res0%>%group_by(motif_name)%>%slice_max(order_by=abs(zscore_rna), n=1)%>%
        as.data.frame()
    
    ### diff motifs
    res2 <- res_motif2%>%filter(comb==ii)%>%group_by(motif_name)%>%
        slice_max(order_by=abs(zscore_motif), n=1)%>%
        as.data.frame()
    ### combine
    res2 <- res2%>%inner_join(res0, by="motif_name")
    res2
})

res <- res%>%drop_na(beta_rna, beta_motif, zscore_rna, zscore_motif)
opfn <- paste(outdir2, "3_TF_activity_Gene.rds", sep="")
write_rds(res, file=opfn)



##############################################################
### correlation between TF activity and TF gene expression ###
##############################################################

rm(list=ls())
outdir2 <- "./4_motif_plots.outs/2_all_response/"

###
fn <- paste(outdir2, "3_TF_activity_Gene.rds", sep="")
x <- read_rds(fn)

 
## motif changes
### df1 <- x%>%dplyr::select(comb, motif_name, beta_motif)
mat1 <- x%>%pivot_wider(id_cols=motif_name, names_from=comb, values_from=zscore_motif, values_fill=NA)%>%
    column_to_rownames(var="motif_name")
id1 <- rownames(mat1)
comb1 <- names(mat1)

id1_sel <- rowSums(is.na(mat1))==0


### TF gene expression changes
## df2 <- x%>%dplyr::select(comb, motif_name, beta_rna)
mat2 <- x%>%pivot_wider(id_cols=motif_name, names_from=comb, values_from=zscore_rna, values_fill=NA)%>%
    column_to_rownames(var="motif_name")
id2 <- rownames(mat2)
comb2 <- names(mat2)

id2_sel <- rowSums(is.na(mat2))==0

id_sel <- id1_sel&id2_sel

identical(id1, id2)
identical(comb1, comb2)



###
###
mat1 <- mat1[id_sel,]
mat2 <- mat2[id_sel,]
nm <- nrow(mat1)
id <- rownames(mat1)

df2 <- map_dfr(1:nm, function(i){
   ##
   id0 <- id[i] 
   b1 <- as.numeric(mat1[i,])
   b2 <- as.numeric(mat2[i,])
   rr <- cor(b1, b2, method="spearman")
   df0 <- data.frame(motif_name=id0, "rr"=rr)
   df0 
})

fn <- paste(outdir2, "3.2_corr.xlsx", sep="")
tmp <- read.xlsx(fn)
names(tmp) <- c("motif_name", "rr_beta")

dfcomb <- df2%>%left_join(tmp, by="motif_name")%>%rename("rr_zscore"="rr")

dfcomb <- dfcomb%>%arrange(desc(rr_zscore))
opfn <- paste(outdir2, "3.2_corr.xlsx", sep="")
write.xlsx(dfcomb, file=opfn)



####
#### box plots show correlation

## fn <- paste(outdir2, "3.2_corr.xlsx", sep="")
## df2 <- read.xlsx(fn)

## ### motif cluster
## fn <- paste(outdir2, "1.0_cl6_gene_reorder.xlsx", sep="")
## df_cl <- read.xlsx(fn)

## dfcomb <- df2%>%inner_join(df_cl, by=c("motif_name"="gene"))

## dfcomb <- dfcomb%>%mutate(cluster2=paste("CL", cluster, sep=""))
## ###

## p <- ggplot(dfcomb, aes(x=cluster2, y=rr, color=cluster2))+
##    geom_boxplot(outlier.shape=NA, lwd=0.25)+
##    geom_jitter(width=0.2, size=0.5)+ 
##    stat_summary(fun=median, geom="point", shape=23, size=1.8, stroke=0.9)+ylim(-1,1)+
##    ggtitle("SCC of LFC between TFs activity and TFs gene")+  
##    theme_bw()+
##    theme(legend.position="none",
##          axis.title=element_blank(),
##          axis.text=element_text(size=10),
##          plot.title=element_text(size=10,hjust=0.5))

## figfn <- paste(outdir2, "Figure3.0_corr_box.png", sep="")
## ggsave(figfn, p, width=450, height=350, units="px", dpi=120)






#####################################################################################
#### Depict response changes of TF activity and TF gene expression seperately  
######################################################################################



###########################################################
### We use all the response motif to make heatmap plots ###
###########################################################


## rm(list=ls())

## outdir <- "./4_motif_plots.outs/2_all_response/"

## ###
## ### get motif_name and motif_ID
## fn <- "./sc_multiome_data/3_motif/2_motif.activities.outs/1_scATAC.motifActivities.rds"
## atac <- read_rds(fn)
## x <- Motifs(atac)
## motifs <- unlist(x@motif.names)
 

###
### differential motifs 
fn <- "./4_motif_plots.outs/3.2_motif_diff_fitInd_JW.rds"
res <- read_rds(fn)%>%
    mutate(comb=paste(MCls, contrast, sep="_"), zscore=beta/stderr,
                      is_sig=ifelse(qval<0.1, 1, 0))

###
### this file contains correlation between response change in TF activity and TF gene
### These motif are response motifs and filtered by removing missing value in activity and gene expression 
fn <- "./4_motif_plots.outs/response_motif_correct/2.2_response_motif_th0.1.txt"
resp_motif <- read.table(fn, header=TRUE)$motif_name%>%unique()


### top 10
res2 <- res%>%filter(motif_name%in%resp_motif)
top10 <- res2%>%group_by(comb)%>%slice_max(order_by=abs(zscore), n=10)%>%pull(motif_name)%>%unique()

res3 <- res2%>%filter(motif_name%in%top10)%>%dplyr::select(comb, motif_name, zscore, is_sig)
res3 <- res3%>%group_by(comb, motif_name)%>%slice_max(order_by=abs(zscore), n=1)%>%as.data.frame()

## res2 <- res%>%filter(gene%in%motif_zinc, contrast=="zinc", is_sig==1)%>%
##     dplyr::select(comb, motif_name, beta, zscore) 

## res3 <- res%>%filter(motif_name=="MTF1", contrast=="zinc")%>%dplyr::select(comb, motif_name, beta, zscore)


### motif matrix
dfmat <- res3%>%pivot_wider(id_cols=motif_name, names_from=comb, values_from=zscore, values_fill=NA)
mat <- dfmat%>%column_to_rownames(var="motif_name")%>%as.matrix()
## rnSel <- rowSums(is.na(mat))==0
## sum(rnSel)

###
### is significance matrix 
dfmat_sig <- res3%>%pivot_wider(id_cols=motif_name, names_from=comb, values_from=is_sig, values_fill=0)
mat_sig <- dfmat_sig%>%column_to_rownames(var="motif_name")%>%as.matrix()
##
## rnSel <- rowSums(is.na(mat_sig))==0
## sum(rnSel)
## length(rnSel)

identical(rownames(mat), rownames(mat_sig))
identical(colnames(mat), colnames(mat_sig))



###
### get colnames and re-order by treats
rn <- colnames(mat)
x <- str_split(rn, "_", simplify=T)
cvt <- data.frame(comb=rn, MCls=paste(x[,1], x[,2], sep="_"), contrast=x[,3])
cvt <- cvt%>%arrange(MCls)

mat2 <- mat[,cvt$comb]


###
### color for heatmap value
y <- as.vector(mat2)
##y0 <- y[abs(y)<2]
## quantile(abs(y), probs=0.99)

mybreak <- c(min(y,na.rm=T), seq(-6, 6, length.out=98), max(y,na.rm=T))

## quantile(abs(y), probs=c(0.9,0.95,0.99))
## range(y)

mycol <- colorRamp2(mybreak, colorRampPalette(rev(brewer.pal(n=7,name="RdBu")))(100))

 
###
### annotation columns
col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
  "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
   "6_Monocyte"="#984ea3", "7_dnT"="black")
col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
       "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3")


x <- str_split(colnames(mat2), "_", simplify=T)
df_col <- data.frame(celltype=paste(x[,1], x[,2], sep="_"), contrast=x[,3])
col_ha <- HeatmapAnnotation(df=df_col, col=list(celltype=col1, contrast=col2),
    annotation_name_gp=gpar(fontsize=10),
    annotation_legend_param=list(
  celltype=list(title_gp=gpar(fontsize=9), labels_gp=gpar(fontsize=9), title="celltype",
                grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm")),
  contrast=list(title_gp=gpar(fontsize=9), labels_gp=gpar(fontsize=9), title="contrast",
                grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm"))))

###
### 111 motifs

p0 <- Heatmap(mat2, col=mycol,
   cluster_rows=T, km=6,
   cluster_columns=F, column_km=2, ##row_split=c(rep("1", 11), rep("2", 20), rep("3", 46)), 
   show_row_names=T, row_names_gp=gpar(fontsize=5.5),
   show_column_names=T, column_names_gp=gpar(fontsize=7),
   show_row_dend=T, row_dend_width=grid::unit(2, "cm"), show_column_dend=T,
   ##column_names_rot=-45,   
   top_annotation=col_ha,
   heatmap_legend_param=list(title=bquote(~italic(Z)~"-score"),
      title_gp=gpar(fontsize=9),
      at=seq(-9, 9, by=3), 
      labels_gp=gpar(fontsize=9),
      grid_width=grid::unit(0.45, "cm"),
      legend_height=grid::unit(7.8, "cm")))
 
### 1.0-using raw z-score, mat3 reroder
### 1.1-using motified z-score, not sig asign to 0, mat3_sig reorder but not sig assign 0
### mat2, set cluster=T
figfn <- paste(outdir2, "Figure1.0_MCls_top10_motif.heatmap.png", sep="")
png(figfn, width=1100, height=1050,res=120)
set.seed(0)
p0 <- draw(p0)
dev.off()
 


## ### keep the order
 
hmap <- Heatmap(mat2, cluster_rows=T, row_km=6, cluster_columns=F)
set.seed(0)
hmap <- draw(hmap)
cl <- row_order(hmap)
lengths(cl)

geneSel <- rownames(mat2)
DF_cl <- NULL
for (i in names(cl)){
   cl_tmp <- data.frame(cluster=i, gene=geneSel[cl[[i]]])
   DF_cl <- rbind(DF_cl, cl_tmp)
}

## newCL <- c("3"="1", "2"="2", "1"="3")
## DF_cl <- DF_cl%>%mutate(cluster2=newCL[cluster])
## DF2 <- DF_cl%>%arrange(cluster2)

opfn <- paste(outdir2, "1.0_cl6_gene_reorder.xlsx", sep="")
write.xlsx(DF_cl, file=opfn, overwrite=T)

## ###
## opfn <- paste(outdir, "1_zscore_motif.rds", sep="")
## write_rds(mat2, file=opfn)


####
#### correlation in TF activity between conditions
## fn <- "./4_motif_plots.outs/3.2_motif_diff_fitInd_JW.rds"
## res <- read_rds(fn)%>%
##     mutate(comb=paste(MCls, contrast, sep="_"), zscore=beta/stderr,
##                       is_sig=ifelse(qval<0.1, 1, 0))

## res2 <- res%>%group_by(comb, motif_name)%>%summarize(zscore=max(zscore, na.rm=TRUE),.groups="drop")%>%ungroup()
## mat <- res2%>%
##     pivot_wider(id_cols=motif_name, names_from=comb, values_from=zscore)%>%
##     column_to_rownames(var="motif_name")%>%
##     as.matrix()

## fn <- "4_motif_plots.outs/response_motif_correct/2.2_response_motif_th0.1.txt"
## resp_motif <- read.table(fn, header=TRUE)$motif_name%>%unique()

## mat2 <- mat[resp_motif,]
## ### correlation data for plot
## corr_mat <- cor(mat2, method="spearman")

## mycol <- colorRamp2(seq(-1, 1, length.out=100), colorRampPalette(rev(brewer.pal(n=7,name="RdBu")))(100))

 
## ###
## ### annotation columns
## col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
##   "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
##    "6_Monocyte"="#984ea3", "7_dnT"="black")
## col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
##        "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3")

## x <- str_split(colnames(corr_mat), "_", simplify=T)
## df_col <- data.frame(celltype=paste(x[,1], x[,2], sep="_"), contrast=x[,3])
## col_ha <- HeatmapAnnotation(df=df_col, col=list(celltype=col1, contrast=col2),
##     annotation_name_gp=gpar(fontsize=10),
##     annotation_legend_param=list(
##   celltype=list(title_gp=gpar(fontsize=10), labels_gp=gpar(fontsize=8), title="celltype",
##                 grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm")),
##   contrast=list(title_gp=gpar(fontsize=10), labels_gp=gpar(fontsize=8), title="contrast",
##                 grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm")))) 


## ###
## ### main heatmap  
## p0 <- Heatmap(corr_mat, col=mycol,
##    cluster_rows=F, cluster_columns=F,
##    show_row_names=F,
##    show_column_names=T, column_names_gp=gpar(fontsize=6),
##    ##column_names_rot=-45,   
##    top_annotation=col_ha,
##    heatmap_legend_param=list(title="rho",
##       title_gp=gpar(fontsize=10),
##       at=seq(-1, 1, by=0.25), 
##       labels_gp=gpar(fontsize=8),
##       grid_width=grid::unit(0.45, "cm"),
##       legend_height=grid::unit(7.8, "cm")), 
##    use_raster=T, raster_device="png")
 
## ###
## figfn <- paste(outdir2, "Figure1.1_corr.heatmap.png", sep="")
## ## ggsave(figfn, plot=p2, device="png", width=680, height=700, units="px", dpi=300) ## ggsave 
## png(figfn, width=700, height=780,res=120)
## p0 <- draw(p0)
## dev.off()




## ##########################
## ### TF gene expression ###
## ##########################

rm(list=ls())
outdir2 <- "./4_motif_plots.outs/2_all_response/"

###
### motifs, cluster
fn <- "./4_motif_plots.outs/2_all_response/1.0_cl6_gene_reorder.xlsx"
DF <- read.xlsx(fn)
names(DF)[2] <- "motif_name"
motifSel <- DF$motif_name

###
### resp motif inlcude gene name
fn2 <- "./4_motif_plots.outs/2_all_response/3.0_TF_gene.xlsx"
DF2 <- read.xlsx(fn2) 

geneSel <- unique(DF2$gene)


#### gene expression
### Differential results
fn <- "./sc_multiome_data/2_Differential/1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn_sampleID+new_treat.rds"
resDiff <- read_rds(fn)%>%
    mutate(comb=paste(MCls, contrast, sep="_"),
           is_sig=ifelse(p.adjusted<0.1, 1, 0), zscore=statistic)
###
resDiff2 <- resDiff%>%filter(gene%in%geneSel)%>%
    dplyr::select(comb, gene, beta_rna=estimate, zscore_rna=zscore, is_sig)

###
### plot data frame
comb <- sort(unique(resDiff2$comb))
res2 <- map_dfr(comb, function(ii){
   ### 
   res0 <- resDiff2%>%filter(comb==ii)

   res0 <- DF2%>%inner_join(res0, by="gene") 
   res0 <- res0%>%group_by(motif_name)%>%slice_max(order_by=abs(zscore_rna), n=1)%>%
       as.data.frame()
   res0
})    


###
### Differential matrix
dfmat <- res2%>%pivot_wider(id_cols=motif_name, names_from=comb, values_from=zscore_rna, values_fill=NA)
mat <- dfmat%>%column_to_rownames(var="motif_name")%>%as.matrix()
mat <- mat[motifSel,]

rnSel <- rowSums(is.na(mat))==0
sum(rnSel)
length(rnSel)


## ###
## ### is significance matrix
## dfmat_sig <- resDiff2%>%pivot_wider(id_cols=gene, names_from=comb, values_from=is_sig)
## is_sig <- dfmat_sig%>%column_to_rownames(var="gene")%>%as.matrix()
## ## is_sig <- is_sig[rnSel, comb]
## ## mat_sig <- mat*is_sig


###
### get colnames and re-order by treats
rn <- colnames(mat)
x <- str_split(rn, "_", simplify=T)
cvt <- data.frame(comb=rn, MCls=paste(x[,1], x[,2], sep="_"), contrast=x[,3])
cvt <- cvt%>%arrange(contrast)

mat2 <- mat[,cvt$comb] 



###
### setting color for heatmap value
y <- as.vector(mat2)
##y0 <- y[abs(y)<2] ### quantile(abs(y), probs=0.99)
mybreak <- c(min(y,na.rm=T), seq(-6, 6, length.out=98), max(y,na.rm=T))

quantile(abs(y), probs=c(0.9,0.95,0.99), na.rm=T)
range(y)

mycol <- colorRamp2(mybreak, colorRampPalette(rev(brewer.pal(n=7,name="RdBu")))(100))

 
###
### annotation columns
col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
  "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
   "6_Monocyte"="#984ea3", "7_dnT"="black")
col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
       "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3")

###
### column annotation
x <- str_split(colnames(mat2), "_", simplify=T)
df_col <- data.frame(celltype=paste(x[,1], x[,2], sep="_"), contrast=x[,3])
col_ha <- HeatmapAnnotation(df=df_col, col=list(celltype=col1, contrast=col2),
    annotation_name_gp=gpar(fontsize=10),
    annotation_legend_param=list(
  celltype=list(title_gp=gpar(fontsize=9), labels_gp=gpar(fontsize=9), title="celltype",
                grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm")),
  contrast=list(title_gp=gpar(fontsize=9), labels_gp=gpar(fontsize=9), title="contrast",
                grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm"))))

###
### row annotation
df_row <- DF%>%select(cluster)
row_ha <- rowAnnotation(cl=anno_block(gp=gpar(fill=2:7), labels=c("CL2", "CL1", "CL3", "CL4", "CL5", "CL6")))
 
table(DF$cluster)
##mat[is.na(mat)] <- 0
p0 <- Heatmap(mat2, col=mycol, 
   cluster_rows=F, 
   cluster_columns=F,
   row_split=c(rep("1",33), rep("2", 8), rep("3", 22), rep("4", 15), rep("5",12), rep("6", 21)),
   ##row_split=c(rep("1", 9), rep("2", 17), rep("3", 38)), 
   show_row_names=T, row_names_gp=gpar(fontsize=5),
   show_column_names=T, column_names_gp=gpar(fontsize=7),
   show_row_dend=F, show_column_dend=F,
   row_title=NULL,
   ##column_names_rot=-45,   
   top_annotation=col_ha,
   left_annotation=row_ha,
   heatmap_legend_param=list(title=bquote(~italic(Z)~"-score"),
      title_gp=gpar(fontsize=9),
      at=seq(-9, 9, by=3), 
      labels_gp=gpar(fontsize=9),
      grid_width=grid::unit(0.45, "cm"),
      legend_height=grid::unit(7.8, "cm")))

###
figfn <- paste(outdir2, "Figure2.0_top10_geneExp.heatmap.png", sep="")
png(figfn, width=1000, height=1050,res=120)
set.seed(0)
p0 <- draw(p0)
dev.off()


## opfn <- paste(outdir, "2_zscore_TF_genes.rds", sep="")
## write_rds(mat, file=opfn)

###########################################################################
### TF activity and TF gene expression together 
##########################################################################


###################################################################
### Heatmap together, combine TF activity and TF gene expression ###
###################################################################


rm(list=ls())

outdir2 <- "./4_motif_plots.outs/2_all_response/"
if ( !file.exists(outdir2) ) dir.create(outdir2, showWarnings=F, recursive=T)


###
### top 10 motifs
fn <- "./4_motif_plots.outs/2_all_response/1.0_cl6_gene_reorder.xlsx"
df_motif <- read.xlsx(fn)
motifSel <- unique(df_motif$gene)

###
### data for TF activity and TF genes
fn <- "./4_motif_plots.outs/2_all_response/3_TF_activity_Gene.rds"
df2 <- read_rds(fn)

x1 <- df2%>%dplyr::select(comb, motif_name, zscore_motif)%>%
    pivot_wider(id_cols=motif_name, names_from=comb, values_from=zscore_motif, values_fill=NA)%>%
    column_to_rownames(var="motif_name")%>%
    as.matrix()

###
x2 <- df2%>%dplyr::select(comb, motif_name, zscore_rna)%>%
    pivot_wider(id_cols=motif_name, names_from=comb, values_from=zscore_rna, values_fill=NA)%>%
    column_to_rownames(var="motif_name")%>%
    as.matrix()

identical(colnames(x1), colnames(x2))

###
### rename x1
comb2 <- paste0("1_", colnames(x1), sep="")
colnames(x1) <- comb2
x1 <- x1[motifSel,]

###
### rename x2
comb2 <- paste0("2_", colnames(x2), sep="")
colnames(x2) <- comb2
x2 <- x2[motifSel,]

identical(rownames(x1), rownames(x2))

### combine matrix
mat <- cbind(x1, x2)


###
### data frame of colnames 
cvt <- str_split(colnames(mat), "_", simplify=T)

cvt2 <- data.frame(comb=colnames(mat), Feature=as.numeric(cvt[,1]),
   MCls=paste(cvt[,2], cvt[,3], sep="_"), contrast=cvt[,4])
cvt2 <- cvt2%>%arrange(Feature, contrast)

mat <- mat[, cvt2$comb]

rnz <- rowSums(is.na(mat))
mat2 <- mat[rnz==0,]

opfn <- paste(outdir2, "3_TF_activity_Gene.mat.rds", sep="")
write_rds(mat2, file=opfn)





#############
### plots ###
#############

rm(list=ls())

outdir2 <- "./4_motif_plots.outs/2_all_response/"
if ( !file.exists(outdir2) ) dir.create(outdir2, showWarnings=F, recursive=T)

fn <- "./4_motif_plots.outs/2_all_response/3_TF_activity_Gene.mat.rds"
mat <- read_rds(fn)

###
### setting color for heatmap value
y <- as.vector(mat)
##y0 <- y[abs(y)<2] ### quantile(abs(y), probs=0.99)
mybreak <- c(min(y), seq(-6, 6, length.out=98), max(y)) 
mycol <- colorRamp2(mybreak, colorRampPalette(rev(brewer.pal(n=7,name="RdBu")))(100))

quantile(abs(y), probs=c(0.9, 0.95, 0.99))


###
### annotation columns
col0 <- c("1"="#8c510a", "2"="#d8b365")
col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
  "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
   "6_Monocyte"="#984ea3", "7_dnT"="black")
col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
       "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3")


x <- str_split(colnames(mat), "_", simplify=T)
x2 <- data.frame(comb=colnames(mat), Feature=as.numeric(x[,1]),
   MCls=paste(x[,2], x[,3], sep="_"), contrast=x[,4])
x2 <- x2%>%arrange(Feature, contrast)


df_col <- x2%>%dplyr::select(celltype=MCls, contrast) ##%>%mutate(Feature=as.character(Feature))

col_ha <- HeatmapAnnotation(df=df_col, col=list(celltype=col1, contrast=col2),
    annotation_label=c("cell type", "contrast"),                        
    annotation_name_gp=gpar(fontsize=11),
    annotation_legend_param=list(       
  celltype=list(title_gp=gpar(fontsize=10), labels_gp=gpar(fontsize=10), title="cell type",
                grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm")),
  contrast=list(title_gp=gpar(fontsize=9), labels_gp=gpar(fontsize=9), title="contrast",
                grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm"))))

comb <- gsub("^[12]_", "", colnames(mat))
colnames(mat) <- comb


## show_legend=c(F,F))
### main plots
p2 <- Heatmap(mat, col=mycol,
   cluster_rows=T, cluster_columns=F, row_km=6, column_split=rep(c("TF activity", "TF genes"), each=48),
   column_title_gp=gpar(fontsize=15),
   show_row_names=T, row_names_gp=gpar(fontsize=5.5),
   show_column_names=T, column_names_gp=gpar(fontsize=8),
   show_row_dend=T, row_dend_width=grid::unit(2, "cm"),
   show_column_dend=F,
   ##column_names_rot=-45,   
   top_annotation=col_ha,
   heatmap_legend_param=list(title=bquote(~italic(Z)~"-score"),
      title_gp=gpar(fontsize=10),
      at=seq(-9, 9, by=3), 
      labels_gp=gpar(fontsize=10),
      grid_width=grid::unit(0.45, "cm"),
      legend_height=grid::unit(7.8, "cm")))
 
###
figfn <- paste(outdir2, "Figure3.1_top10_comb2_cl6.heatmap.png", sep="")
## ggsave(figfn, plot=p2, device="png", width=680, height=700, units="px", dpi=300) ## ggsave 
png(figfn, width=1600, height=1200,res=120)
set.seed(0)
p2 <- draw(p2)
dev.off()


###
### TF clusters
hmap <- Heatmap(mat, cluster_rows=T, row_km=6, cluster_columns=F)
set.seed(0)
hmap <- draw(hmap)
cl <- row_order(hmap)
lengths(cl)

geneSel <- rownames(mat)
DF_cl <- NULL
for (i in names(cl)){
   cl_tmp <- data.frame(cluster=i, gene=geneSel[cl[[i]]])
   DF_cl <- rbind(DF_cl, cl_tmp)
}

## newCL <- c("3"="1", "2"="2", "1"="3")
## DF_cl <- DF_cl%>%mutate(cluster2=newCL[cluster])
## DF2 <- DF_cl%>%arrange(cluster2)
fn <- paste(outdir2, "3.2_corr.xlsx", sep="")
DF_rr <- read.xlsx(fn)

DFcomb <- DF_cl%>%left_join(DF_rr, by=c("gene"="motif_name"))

opfn <- paste(outdir2, "3.2_corr_reorder.xlsx", sep="")
### order from combination TFs and TF gene expression
write.xlsx(DFcomb, file=opfn, overwrite=T)




######################################
#### boxplots correlation 
######################################


fn <- paste(outdir2, "3.2_corr_reorder.xlsx", sep="")
df0 <- read.xlsx(fn)
df2 <- df0%>%mutate(cluster2=paste("CL", cluster, sep=""))
###
  
p <- ggplot(df2, aes(x=cluster2, y=rr_zscore, color=cluster2))+
   geom_boxplot(outlier.shape=NA, lwd=0.25)+
   geom_jitter(width=0.2, size=0.5)+ 
   stat_summary(fun=median, geom="point", shape=23, size=1.8, stroke=0.9)+ylim(-1,1)+
   ggtitle("SCC of zscore between TFs activity and TFs gene")+  
   theme_bw()+
   theme(legend.position="none",
         axis.title=element_blank(),
         axis.text=element_text(size=10),
         plot.title=element_text(size=10,hjust=0.5))

figfn <- paste(outdir2, "Figure3.0_corr_box.png", sep="")
ggsave(figfn, p, width=450, height=350, units="px", dpi=120)





### TF motif vs TF genes

######################
#### CCA analysis ####
######################

rm(list=ls())

outdir2 <- "./4_motif_plots.outs/2_all_response/"

fn <- "./4_motif_plots.outs/2_all_response/3_TF_activity_Gene.mat.rds"
x <- read_rds(fn)
comb2 <- gsub("^[12]_", "", colnames(x))
colnames(x) <- comb2


###
x1 <- x[,1:48]
x1 <- scale(x1)

###
x2 <- x[,49:96]
x2 <- scale(x2)


###
### cca analysis
cc0 <- cc(x1, x2)
cc2 <- comput(x1, x2, cc0)



####
####

###
motif <- as.data.frame(cc2$xscores)
names(motif) <- paste("S", 1:48, sep="")
###
rna <- as.data.frame(cc2$yscores)
names(rna) <- paste("S", 1:48, sep="")


identical(rownames(motif), rownames(rna))



###
### custom function 
my_box <- function(data, mapping, shape,...){
   ##
   p <- ggplot(data=data, mapping=mapping)+
       geom_boxplot(outlier.shape=NA, lwd=0.25)+coord_flip()
   ## geom_jitter(data=data, mapping=mapping, width=0.2, size=0.5, shape=23)+coord_flip()## + 
   ## stat_summary(fun=median, geom="point", shape=23, size=1.8, stroke=0.9) ###ylim(-1,1)+
   p
}    


fn <- paste(outdir2, "3.2_corr_reorder.xlsx", sep="")
df_cl <- read.xlsx(fn)


###
### The correlation of score from motif and TF-genes

## df2 <- data.frame(x=1:48, rr=cc0$cor)
## ###
## p0 <- ggplot(df2, aes(x=x, y=rr))+
##    geom_point(size=1.5, shape=1, color="#dd1c77")+
##    guides(color=guide_legend(override.aes=list(size=2)))+ 
##    xlab("Canonical variable")+ylab("Correlation (TF motif-genes)")+
##    theme_bw()+
##    theme(axis.title=element_text(size=9),
##          axis.text=element_text(size=9))

## figfn <- paste(outdir2, "Figure3.1_corr_score.scatter.png", sep="")
## ggsave(figfn, p0, width=420, height=420, units="px", dpi=120)





###
### plot data for motif

plotDF <- motif[,1:5]%>%rownames_to_column(var="gene")%>%
    inner_join(df_cl, by="gene")%>%
    mutate(cluster2=paste("CL", cluster, sep=""))


###
###
p0 <- ggpairs(plotDF, columns=2:6, aes(color=cluster2),
              upper=list(continuous=wrap("cor", size=2)),
              diag=list(continuous=my_box),
              lower=list(continuous=wrap("points", size=0.8)))+
    theme_bw()+
    theme(strip.text=element_text(size=12))
    

figfn <- paste(outdir2, "Figure3.2_motif_score.pair.png", sep="")
ggsave(figfn, plot=p0, width=800, height=720, units="px", dpi=120)



###
### plot data for TF-regulated genes

plotDF2 <- rna[,1:5]%>%rownames_to_column(var="gene")%>%
    inner_join(df_cl, by="gene")%>%
    mutate(cluster2=paste("CL", cluster, sep=""))


###
###
p2 <- ggpairs(plotDF2, columns=2:6, aes(color=cluster2),
              upper=list(continuous=wrap("cor", size=2)),
              diag=list(continuous=my_box),
              lower=list(continuous=wrap("points", size=0.8)))+
    theme_bw()+
    theme(strip.text=element_text(size=12))
    

figfn <- paste(outdir2, "Figure3.2_TF-gene_score.pair.png", sep="")
ggsave(figfn, plot=p2, width=780, height=720, units="px", dpi=120)


###
###
plotDF3 <- map_dfr(1:5, function(i){
   ##
   df0 <- data.frame(x=motif[,i], y=rna[,i])
   df0$gr <- paste("Score-", i, sep="")
   df0
})
 
p3 <- ggplot(plotDF3, aes(x, y))+
   geom_point(size=1.5, shape=1, color="#dd1c77")+
   xlab("The score of TF activity")+ylab("The score of TF genes")+ 
   facet_wrap(~gr, nrow=3)+
   theme_bw()+
   theme(axis.title=element_text(size=10),
         axis.text=element_text(size=10))
##
figfn <- paste(outdir2, "Figure3.2_score.xy.png", sep="")
ggsave(figfn, p3, width=420, height=540, units="px", dpi=120)






##############################
### Heatmap score from cca ###
##############################

## motif <- cc2$xscores
## rna <- cc2$yscores

## m1 <- motif
## colnames(m1) <- paste("S", 1:48, sep="")
## m2 <- rna
## colnames(m2) <- paste("S", 1:48, sep="")

## mat <- cbind(m1, m2)



## #### setting colors
## y <- as.vector(mat)
## ##y0 <- y[abs(y)<2] ### quantile(abs(y), probs=0.99)
## mybreak <- c(min(y), seq(-3, 3, length.out=98), max(y)) 
## mycol <- colorRamp2(mybreak, colorRampPalette(rev(brewer.pal(n=7,name="RdBu")))(100))

## quantile(abs(y), probs=c(0.9, 0.95, 0.99))

## ### main plots
## p2 <- Heatmap(mat, col=mycol,
##    cluster_rows=T, row_km=6, cluster_columns=F,
##    column_split=rep(c("TF activity", "TF genes"), each=48),
##    show_row_dend=T, row_dend_width=grid::unit(2, "cm"),   
##    show_row_names=T, row_names_gp=gpar(fontsize=6),
##    show_column_names=T, column_names_gp=gpar(fontsize=6),
##    ## show_column_names=F, column_names_gp=gpar(fontsize=6),
##    ###top_annotation=col_ha,
##    heatmap_legend_param=list(title="score",
##       title_gp=gpar(fontsize=9),
##       at=seq(-4.5, 4.5, by=1.5), 
##       labels_gp=gpar(fontsize=9),
##       grid_width=grid::unit(0.45, "cm"),
##       legend_height=grid::unit(7.8, "cm")))
 
## ###
## figfn <- paste(outdir2, "Figure3.3_score_comb2.heatmap.png", sep="")
## ## ggsave(figfn, plot=p2, device="png", width=680, height=700, units="px", dpi=300) ## ggsave 
## png(figfn, width=1600, height=1200,res=120)
## set.seed(0)
## p2 <- draw(p2)
## dev.off()



########################
### loading from cca ###
########################

## m1 <- t(cc2$corr.X.xscores)
## m2 <- t(cc2$corr.Y.yscores)

## rownames(m1) <- paste("S", 1:48, sep="")
## rownames(m2) <- paste("S", 1:48, sep="")


## #### setting colors
## y <- as.vector(m1)
## ##y0 <- y[abs(y)<2] ### quantile(abs(y), probs=0.99)
## ## mybreak <- c(min(y), seq(-1, 1, length.out=98), max(y))
## mybreak <- seq(-1, 1, length.out=100)
## mycol <- colorRamp2(mybreak, colorRampPalette(rev(brewer.pal(n=7,name="RdBu")))(100))

## ## quantile(abs(y), probs=c(0.9, 0.95, 0.99))


## ###
## ### annotation columns
## col0 <- c("1"="#8c510a", "2"="#d8b365")
## col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
##   "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
##    "6_Monocyte"="#984ea3", "7_dnT"="black")
## col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
##        "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3")


## cvt <- str_split(colnames(m1), "_", simplify=T)
## cvt2 <- data.frame(comb=colnames(m1),
##    MCls=paste(cvt[,1], cvt[,2], sep="_"), contrast=cvt[,3])
## cvt2 <- cvt2%>%arrange(contrast)


## df_col <- cvt2%>%dplyr::select(MCls, contrast) ##%>%mutate(Feature=as.character(Feature))

## col_ha <- HeatmapAnnotation(df=df_col, col=list(MCls=col1, contrast=col2),
##     annotation_name_gp=gpar(fontsize=10),
##     annotation_legend_param=list(       
##   MCls=list(title_gp=gpar(fontsize=9), labels_gp=gpar(fontsize=9), title="celltype",
##                 grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm")),
##   contrast=list(title_gp=gpar(fontsize=9), labels_gp=gpar(fontsize=9), title="contrast",
##                 grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm"))))

## ### plot data
## m1a <- m1[,cvt2$comb]


## ### main plots
## p1 <- Heatmap(m1a, col=mycol,
##    cluster_rows=F, cluster_columns=F,
##    show_column_dend=F, column_dend_height=grid::unit(2, "cm"),   
##    show_row_names=T, row_names_gp=gpar(fontsize=7),
##    show_column_names=T, column_names_gp=gpar(fontsize=6),
##    top_annotation=col_ha,
##    heatmap_legend_param=list(title="loading",
##       title_gp=gpar(fontsize=9),
##       at=seq(-1, 1, by=0.5), 
##       labels_gp=gpar(fontsize=9),
##       grid_width=grid::unit(0.45, "cm"),
##       legend_height=grid::unit(7.8, "cm")))
 
## ###
## figfn <- paste(outdir2, "Figure3.4_loading_motif.heatmap.png", sep="")
## ## ggsave(figfn, plot=p2, device="png", width=680, height=700, units="px", dpi=300) ## ggsave 
## png(figfn, width=720, height=720,res=120)
## set.seed(0)
## p1 <- draw(p1)
## dev.off()



## ###
## ### loading of RNA

## #### setting colors
## y <- as.vector(m2)
## ##y0 <- y[abs(y)<2] ### quantile(abs(y), probs=0.99)
## ## mybreak <- c(min(y), seq(-1, 1, length.out=98), max(y))
## mybreak <- seq(-1, 1, length.out=100)
## mycol <- colorRamp2(mybreak, colorRampPalette(rev(brewer.pal(n=7,name="RdBu")))(100))

## ## quantile(abs(y), probs=c(0.9, 0.95, 0.99))


## ###
## ### annotation columns
## col0 <- c("1"="#8c510a", "2"="#d8b365")
## col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
##   "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
##    "6_Monocyte"="#984ea3", "7_dnT"="black")
## col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
##        "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3")


## cvt <- str_split(colnames(m2), "_", simplify=T)
## cvt2 <- data.frame(comb=colnames(m2),
##    MCls=paste(cvt[,1], cvt[,2], sep="_"), contrast=cvt[,3])
## cvt2 <- cvt2%>%arrange(contrast)


## df_col <- cvt2%>%dplyr::select(MCls, contrast) ##%>%mutate(Feature=as.character(Feature))

## col_ha <- HeatmapAnnotation(df=df_col, col=list(MCls=col1, contrast=col2),
##     annotation_name_gp=gpar(fontsize=10),
##     annotation_legend_param=list(       
##   MCls=list(title_gp=gpar(fontsize=9), labels_gp=gpar(fontsize=9), title="celltype",
##                 grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm")),
##   contrast=list(title_gp=gpar(fontsize=9), labels_gp=gpar(fontsize=9), title="contrast",
##                 grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm"))))



## ### plot data
## m2a <- m2[,cvt2$comb]

## ### main plots
## p2 <- Heatmap(m2a, col=mycol,
##    cluster_rows=F, cluster_columns=F,
##    show_column_dend=F, column_dend_height=grid::unit(2, "cm"),   
##    show_row_names=T, row_names_gp=gpar(fontsize=7),
##    show_column_names=T, column_names_gp=gpar(fontsize=6),
##    top_annotation=col_ha,
##    heatmap_legend_param=list(title="loading",
##       title_gp=gpar(fontsize=9),
##       at=seq(-1, 1, by=0.5), 
##       labels_gp=gpar(fontsize=9),
##       grid_width=grid::unit(0.45, "cm"),
##       legend_height=grid::unit(7.8, "cm")))
 
## ###
## figfn <- paste(outdir2, "Figure3.4_loading_TF-gene.heatmap.png", sep="")
## ## ggsave(figfn, plot=p2, device="png", width=680, height=700, units="px", dpi=300) ## ggsave 
## png(figfn, width=720, height=720,res=120)
## set.seed(0)
## p2 <- draw(p2)
## dev.off()





################################
### scatter plots of loading ###
################################

### custom function 
my_box <- function(data, mapping,...){
   ##
   p <- ggplot(data=data, mapping=mapping)+
       geom_boxplot(outlier.shape=NA, lwd=0.25)+coord_flip()
   ## geom_jitter(data=data, mapping=mapping, width=0.2, size=0.5, shape=23)+coord_flip()## + 
   ## stat_summary(fun=median, geom="point", shape=23, size=1.8, stroke=0.9) ###ylim(-1,1)+
   p
}


m1 <- as.data.frame(cc2$corr.X.xscores)
m2 <- as.data.frame(cc2$corr.Y.yscores)

colnames(m1) <- paste("S", 1:48, sep="")
colnames(m2) <- paste("S", 1:48, sep="")


col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
  "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
   "6_Monocyte"="#984ea3", "7_dnT"="black")
col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
       "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3")


###
### plot data for motif

plotDF <- m1[,1:5]%>%as.data.frame()%>%rownames_to_column(var="comb")
x <- str_split(plotDF$comb, "_", simplify=T)
plotDF <- plotDF%>%mutate(MCls=paste(x[,1], x[,2], sep="_"), contrast=x[,3])


###
###
p0 <- ggpairs(plotDF, columns=2:6, aes(color=contrast),
              upper="blank", ###list(continuous=wrap("cor", size=2)),
              diag=list(continuous=my_box),
              lower=list(continuous=wrap("points", size=0.8)))+
    scale_color_manual(values=col2)+
    theme_bw()+
    theme(strip.text=element_text(size=12))
    

figfn <- paste(outdir2, "Figure3.3_motif_loading.pair.png", sep="")
ggsave(figfn, plot=p0, width=800, height=720, units="px", dpi=120)



###
### plot data for TF-regulated genes

plotDF2 <- m2[,1:5]%>%as.data.frame()%>%rownames_to_column(var="comb")
x <- str_split(plotDF2$comb, "_", simplify=T)
plotDF2 <- plotDF2%>%mutate(MCls=paste(x[,1], x[,2], sep="_"), contrast=x[,3])


###
###
p2 <- ggpairs(plotDF2, columns=2:6, aes(color=contrast),
              upper="blank",   ###list(continuous=wrap("cor", size=2)),
              diag=list(continuous=my_box),
              lower=list(continuous=wrap("points", size=0.8)))+
    scale_color_manual(values=col2)+
    theme_bw()+
    theme(strip.text=element_text(size=12))
    

figfn <- paste(outdir2, "Figure3.3_TF-gene_loading.pair.png", sep="")
ggsave(figfn, plot=p2, width=780, height=720, units="px", dpi=120)




#########################################################
#### cluster CCA for TF-activity and TF-genes
#############################################################


rm(list=ls())
outdir2 <- "./4_motif_plots.outs/2_all_response/cluster_CCA/"
if ( !file.exists(outdir2) ) dir.create(outdir2, showWarnings=F, recursive=F)



### get covariance-variance matrix
getVar <- function(X, Y, df_cl){
    
   ### clusters    
   cls <- unique(df_cl$cluster)
   mm <- sum(table(df_cl$cluster)^2)
    
   ## variance, var_xx
   var_xx <- lapply(cls, function(i){
       ##
       rnSel <- df_cl%>%dplyr::filter(cluster==i)%>%pull(gene)
       n0 <- length(rnSel)
       x <- X[rnSel,]
       x <- scale(x)
       xx <- crossprod(x)*n0
       xx
   })
   var_xx <- Reduce("+", var_xx)
   var_xx <- var_xx/mm 

   ## variance, var_yy
   var_yy <- lapply(cls, function(i){
      ##
      rnSel <- df_cl%>%filter(cluster==i)%>%pull(gene)
      n0 <- length(rnSel)
      y <- Y[rnSel,]
      y <- scale(y) 
      yy <- crossprod(y)*n0
      yy
   })
   var_yy <- Reduce("+", var_yy)
   var_yy <- var_yy/mm

   ### 
   ### covariance, var_xy
   var_xy <- lapply(cls, function(i){
      ###
      rnSel <- df_cl%>%filter(cluster==i)%>%pull(gene)
      n0 <- length(rnSel)
      x <- X[rnSel,]
      x <- scale(x) 
      y <- Y[rnSel,]
      y <- scale(y) 
      np <- ncol(x) 
      ###
       
      xy <- matrix(0, np, np)   
      for (i in nrow(x)){
         ##
         for ( j in nrow(y)){
         ##
           x0 <- matrix(as.numeric(x[i,]), np, 1)
           y0 <- matrix(as.numeric(y[j,]), 1, np)
           xy <- xy+x0%*%y0
         }    
      }       
      xy
    })
    var_xy <- Reduce("+", var_xy)
    var_xy <- var_xy/mm

    ### return results
    results <- list(var_xx=var_xx, var_yy=var_yy, var_xy=var_xy)
    results         
}    


###
###
getCCA <- function(Var_list){
   ###
   xx <- Var_list$var_xx
   yy <- Var_list$var_yy
   xy <- Var_list$var_xy
    
   ### matrix
   sx <- sqrtm(solve(xx))
   sy <- sqrtm(solve(yy))
   T <- sx%*%xy%*%sy
   M <- tcrossprod(T) 

   results <- eigen(M)  
   uu <- results$vectors
   lam <- results$values
   ii <- lam>0  
   ### weights for x
   wa <- sx%*%uu[,ii]

   ### weights for y 
   wb <- sy%*%sy%*%t(xy)%*%wa
   lam_inv <- 1/sqrt(lam[ii])
   wb <- sweep(wb, 2, lam_inv, "*")
    
   ###
   results$wa <- wa
   results$wb <- wb
   results
}    
    




### cluster infor
fn <- "./4_motif_plots.outs/2_all_response/3.2_corr_reorder.xlsx"
df_cl <- read.xlsx(fn)
## df_cl2 <- df_cl%>%filter(gene%in%rnSel)



###
fn <- "./4_motif_plots.outs/2_all_response/3_TF_activity_Gene.mat.rds"
x <- read_rds(fn)
rnSel <- rownames(x)

comb2 <- gsub("^[12]_", "", colnames(x))
colnames(x) <- comb2




###
x1 <- x[,1:48]
x1 <- scale(x1)

###
x2 <- x[,49:96]
x2 <- scale(x2)
 

### cca object
Var <- getVar(x1, x2, df_cl)
cca <- getCCA(Var)
opfn <- paste(outdir2, "3_cca_TF-gene.rds", sep="")
write_rds(cca, file=opfn)

wa <- cca$wa
xscore <- x1%*%wa
colnames(xscore) <- paste("S", 1:ncol(xscore), sep="")
##
wb <- cca$wb
yscore <- x2%*%wb
colnames(yscore) <- paste("S", 1:ncol(yscore), sep="")



### custom function 
my_box <- function(data, mapping,...){
   ##
   p <- ggplot(data=data, mapping=mapping)+
       geom_boxplot(outlier.shape=1, lwd=0.25)+coord_flip()
   ## geom_jitter(data=data, mapping=mapping, width=0.2, size=0.5, shape=23)+coord_flip()## + 
   ## stat_summary(fun=median, geom="point", shape=23, size=1.8, stroke=0.9) ###ylim(-1,1)+
   p
}    



###
### plot data for motif
 
plotDF <- xscore[,1:5]%>%as.data.frame()%>%rownames_to_column(var="gene")%>%
    inner_join(df_cl, by="gene")%>%
    mutate(cluster2=paste("CL", cluster, sep=""))


###
###
p0 <- ggpairs(plotDF, columns=2:6, aes(color=cluster2),
              upper=list(continuous=wrap("cor", size=2)),
              diag=list(continuous=my_box),
              lower=list(continuous=wrap("points", size=0.8)))+
    theme_bw()+
    theme(strip.text=element_text(size=12))
    

figfn <- paste(outdir2, "Figure3.2_cluster_motif_score.pair.png", sep="")
ggsave(figfn, plot=p0, width=800, height=720, units="px", dpi=120)



###
### plot data for TF-regulated genes

plotDF2 <- yscore[,1:5]%>%as.data.frame()%>%rownames_to_column(var="gene")%>%
    inner_join(df_cl, by="gene")%>%
    mutate(cluster2=paste("CL", cluster, sep=""))


###
###
p2 <- ggpairs(plotDF2, columns=2:6, aes(color=cluster2),
              upper=list(continuous=wrap("cor", size=2)),
              diag=list(continuous=my_box),
              lower=list(continuous=wrap("points", size=0.8)))+
    theme_bw()+
    theme(strip.text=element_text(size=12))
    

figfn <- paste(outdir2, "Figure3.2_cluster_TF-gene_score.pair.png", sep="")
ggsave(figfn, plot=p2, width=780, height=720, units="px", dpi=120)



###
###
plotDF3 <- map_dfr(1:5, function(i){
   ##
   df0 <- data.frame(x=xscore[,i], y=yscore[,i])
   df0$gr <- paste("Score-", i, sep="")
   df0
})
 
p3 <- ggplot(plotDF3, aes(x, y))+
   geom_point(size=1.5, shape=1, color="#dd1c77")+
   xlab("The score of TF activity")+ylab("The score of TF genes")+
   geom_abline(color="grey")+ 
   facet_wrap(~gr, nrow=3)+
   theme_bw()+
   theme(axis.title=element_text(size=10),
         axis.text=element_text(size=10))
##
figfn <- paste(outdir2, "Figure3.2_cluster_score.xy.png", sep="")
ggsave(figfn, p3, width=420, height=540, units="px", dpi=120)


###############################################################
### correlation between score and original variables
################################################################
 

### custom function 
my_box <- function(data, mapping,...){
   ##
   p <- ggplot(data=data, mapping=mapping)+
       geom_boxplot(outlier.shape=1, lwd=0.25)+coord_flip()
   ## geom_jitter(data=data, mapping=mapping, width=0.2, size=0.5, shape=23)+coord_flip()## + 
   ## stat_summary(fun=median, geom="point", shape=23, size=1.8, stroke=0.9) ###ylim(-1,1)+
   p
}


m1 <- as.data.frame(t(cor(xscore, x1)))
m2 <- as.data.frame(t(cor(yscore, x2)))

ns <- ncol(m1)


colnames(m1) <- paste("S", 1:ns, sep="")
colnames(m2) <- paste("S", 1:ns, sep="")


col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
  "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
   "6_Monocyte"="#984ea3", "7_dnT"="black")
col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
       "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3")


###
### plot data for motif

plotDF <- m1[,1:5]%>%rownames_to_column(var="comb")
x <- str_split(plotDF$comb, "_", simplify=TRUE)
plotDF <- plotDF%>%mutate(MCls=paste(x[,1], x[,2], sep="_"), contrast=x[,3])


###
###
p0 <- ggpairs(plotDF, columns=2:6, aes(color=contrast),
              upper="blank",   ###list(continuous=wrap("cor", size=2)),
              diag=list(continuous=my_box),
              lower=list(continuous=wrap("points", size=0.8)))+
    scale_color_manual(values=col2)+
    theme_bw()+
    theme(strip.text=element_text(size=12))
    

figfn <- paste(outdir2, "Figure3.3_cluster_motif_loading.pair.png", sep="")
ggsave(figfn, plot=p0, width=800, height=720, units="px", dpi=120)



###
### plot data for TF-regulated genes

plotDF2 <- m2[,1:5]%>%as.data.frame()%>%rownames_to_column(var="comb")
x <- str_split(plotDF2$comb, "_", simplify=TRUE)
plotDF2 <- plotDF2%>%mutate(MCls=paste(x[,1], x[,2], sep="_"), contrast=x[,3])


###
###
p2 <- ggpairs(plotDF2, columns=2:6, aes(color=contrast),
              upper="blank",   ##list(continuous=wrap("cor", size=2)),
              diag=list(continuous=my_box),
              lower=list(continuous=wrap("points", size=0.8)))+
    scale_color_manual(values=col2)+
    theme_bw()+
    theme(strip.text=element_text(size=12))
    

figfn <- paste(outdir2, "Figure3.3_cluster_TF-gene_loading.pair.png", sep="")
ggsave(figfn, plot=p2, width=780, height=720, units="px", dpi=120)




########################################################################
### compare the results between R function cc and my function getCCA ###
########################################################################



###
### cca analysis
cc0 <- cc(x1, x2)
cc2 <- comput(x1, x2, cc0)






var_xx <- var(x1,x1)
var_yy <- var(x2, x2)
var_xy <- var(x1, x2)

Var_list <- list(var_xx=var_xx, var_yy=var_yy, var_xy=var_xy)
cc_my <- getCCA(Var_list)

wa <- cc_my$wa
xscore <- x1%*%wa

##
wb <- cc_my$wb
yscore <- x2%*%wb


###
###
xscore_0 <- cc2$xscores
plotDF <- data.frame(x=xscore_0[,1], y=xscore[,1])
p0 <- ggplot(plotDF, aes(x=x, y=y))+
   geom_point(size=1.5, shape=1, color="#dd1c77")+
   xlab("xscore of cc fun")+ylab("xscore of my fun")+
   theme_bw()

figfn <- paste(outdir2, "Figure3_888_xscore.png", sep="")
ggsave(figfn, p0, width=320, height=320, units="px", dpi=120)


### compare
yscore_0 <- cc2$yscore
plotDF2 <- data.frame(x=yscore_0[,1], y=yscore[,1])
p2 <- ggplot(plotDF, aes(x=x, y=y))+
   geom_point(size=1.5, shape=1, color="#dd1c77")+
   xlab("yscore of cc fun")+ylab("yscore of my fun")+
   theme_bw()

figfn <- paste(outdir2, "Figure3_888_yscore.png", sep="")
ggsave(figfn, p2, width=320, height=320, units="px", dpi=120)











##########################################################
### Correlation of TF activity and TF gene expression ####
##########################################################


rm(list=ls())

outdir2 <- "./4_motif_plots.outs/2_all_response/Examples_TF/"
if ( !file.exists(outdir2)) dir.create(outdir2, showWarnings=F, recursive=T)




fn <- "./4_motif_plots.outs/2_all_response/1.0_cl6_gene_reorder.xlsx"
df_cl <- read.xlsx(fn)
###
fn <- "./4_motif_plots.outs/2_all_response/3.2_corr.xlsx"
df_corr <- read.xlsx(fn)

df2 <- df_cl%>%inner_join(df_corr, by=c("gene"="motif_name"))%>%arrange(desc(rr))


gene0 <- "IRF7"
## gene1 <- "IRF7"


fn <- "./4_motif_plots.outs/3.2_motif_diff_fitInd_JW.rds"
res_motif <- read_rds(fn)%>%mutate(comb=paste(MCls, contrast, sep="_"), zscore=beta/stderr)
res2_motif <- res_motif%>%filter(motif_name==gene0)%>%
    dplyr::select(comb, motif_name, beta_TF=beta, zscore_TF=zscore)  



#### gene expression
### Differential results
fn <- "./sc_multiome_data/2_Differential/1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn_sampleID+new_treat.rds"
res_gene <- read_rds(fn)%>%mutate(comb=paste(MCls, contrast, sep="_"))
res2_gene <- res_gene%>%dplyr::filter(gene==gene0)%>%
    dplyr::select(comb,  beta_gene=estimate, zscore_gene=statistic)


res2_comb <- res2_motif%>%inner_join(res2_gene, by="comb")

 

##############################################
### box and violin plots of peaks and gene ###
##############################################



###
### setting colors
col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
  "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
   "6_Monocyte"="#984ea3", "7_dnT"="black")


col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
       "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3",
       "ctrl"="grey", "ctrl_etOH"="grey", "ctrl_water"="grey")

###
### differential gene expression data
fn <- "./sc_multiome_data/2_Differential/1_DiffRNA_2.outs/1_YtX.sel_clusters_0.1_cn.rds"
mat <- read_rds(fn)
depth <- colSums(mat)
mat_gene <- sweep(mat, 2, depth, "/")
mat_gene <- log2(mat_gene*1e+06+1)



###
### motif infor

fn <- "./4_motif_plots.outs/3.2_motif_diff_fitInd_JW.rds"
DF_motif <- read_rds(fn)%>%dplyr::select(motif_ID=gene, motif_name)%>%distinct(motif_ID, .keep_all=T)

fn <- "./sc_multiome_data/3_motif/2_motif.activities.outs/2_motif.ave_macs2_0.1_cn_control.rds"
mat_motif <- read_rds(fn)




### 
getData <- function(mat, gene, comb){

###    
x <- str_split(colnames(mat), "_", simplify=T)
cvt <- data.frame(rn=colnames(mat), MCls=paste(x[,1], x[,2], sep="_"), treats=x[,3], ind=x[,4])
cvt <- cvt%>%mutate(treat2=ifelse(treats%in%c("etOH", "water"), "control", treats),
                    comb=paste(MCls, treats, sep="_"),
                    comb2=paste(MCls, treat2, sep="_"))
    
###    
split0 <-  unlist(str_split(comb, "_"))
oneMCl <- paste(split0[1], split0[2], sep="_")

###    
### treat from comb column    
comb2 <- paste(rep(oneMCl, 3), c("etOH", "water", "control"), sep="_")         
DF <- map_dfr(c(comb, comb2[1:2]), function(ii){    
###
split0 <-  unlist(str_split(ii, "_"))
##
trt <- split0[3]
trt2 <- ifelse(grepl("etOH|water", trt), paste0("ctrl_", trt), trt)    
##    
rnSel <- cvt%>%dplyr::filter(comb==ii)%>%pull(rn) 
if ( length(rnSel)>0 ){    
   y2 <- mat[gene, rnSel]
   names(y2) <- gsub(".*_", "", names(y2))    
   df0 <- data.frame(y=y2, ind=names(y2), treats=trt2)
   rownames(df0) <- NULL
}else{
   df0 <- NULL
}
df0
}) ### End map_dfr

##
##control from comb2 column
ii <- comb2[3]
split0 <-  unlist(str_split(ii, "_"))
trt <- "ctrl"    
rnSel2 <- cvt%>%dplyr::filter(comb2==ii)%>%pull(rn)    
if ( length(rnSel2)>0 ){    
   y2 <- mat[gene, rnSel2]
   names(y2) <- gsub(".*_", "", names(y2))
   y2 <- tapply(y2, names(y2), mean, na.rm=T)
   df0 <- data.frame(y=y2, ind=names(y2), treats=trt)
}else{
  df0 <- NULL
}  

### combine
DF <- rbind(DF, df0)
DF
}



###
### gene expression 

gene0 <- "ELF1"
ii <- "3_TEM_nicotine"
    
plotDF <- getData(mat_gene, gene0, comb=ii)    
p1 <- ggplot(plotDF, aes(x=factor(treats), y=y, color=treats))+
   geom_boxplot(outlier.shape=NA, lwd=0.3)+
   geom_jitter(width=0.2, size=1)+
   stat_summary(fun=median, geom="point", shape=23, size=5, stroke=0.9)+ 
   scale_y_continuous("Normalized expression",
      expand=expansion(mult=c(0.2, 0.2)))+
   scale_color_manual("", values=col2)+ 
   theme_bw()+
   theme(## axis.text.x=element_text(angle=-90, hjust=0, size=10),
         axis.text=element_text(size=10),
         axis.title.x=element_blank(),
         axis.title.y=element_text(size=12),
         legend.position="none")
  
figfn <- paste(outdir2, "Figure1.1_", gene0, "_", ii, "_gene.box.png", sep="")
ggsave(figfn, p1,  width=480, height=360, units="px", dpi=120)





###
### TF motifs 
 
gene0 <- "ELF1"
m0_id <- DF_motif%>%filter(motif_name==gene0)%>%pull(motif_ID)
ii <- "3_TEM_nicotine"
  
###
plotDF2 <- getData(mat_motif, m0_id, comb=ii)
p2 <- ggplot(plotDF2, aes(x=factor(treats), y=y, color=treats))+
   geom_boxplot(outlier.shape=NA, lwd=0.3)+
   geom_jitter(width=0.2, size=1)+
   stat_summary(fun=median, geom="point", shape=23, size=5, stroke=0.9)+ 
   scale_y_continuous("TF activity",
      expand=expansion(mult=c(0.2, 0.2)))+
   scale_color_manual("", values=col2)+ 
   theme_bw()+
   theme(## axis.text.x=element_text(angle=-90, hjust=0, size=10),
         axis.text=element_text(size=10),
         axis.title.x=element_blank(),
         axis.title.y=element_text(size=12),
         legend.position="none")
  
figfn <- paste(outdir2, "Figure1.2_", gene0, "_", ii, "_TF.box.png", sep="")
ggsave(figfn, p2,  width=420, height=360, units="px", dpi=120)




### scatter plots
trt <- gsub(".*_", "", ii)

## gene
df0 <- plotDF%>%filter(treats==trt)%>%dplyr::select(ind, y2=y)
df1 <- plotDF%>%filter(treats=="ctrl")%>%dplyr::select(ind, y0=y)

df2 <- df0%>%full_join(df1, by="ind")%>%mutate(y_diff=y2-y0)%>%
   dplyr::select(ind, y_diff_gene=y_diff) 


### motif
df0 <- plotDF2%>%filter(treats==trt)%>%dplyr::select(ind, y2=y)
df1 <- plotDF2%>%filter(treats=="ctrl")%>%dplyr::select(ind, y0=y)

df3 <- df0%>%full_join(df1, by="ind")%>%mutate(y_diff=y2-y0)%>%
   dplyr::select(ind, y_diff_motif=y_diff) 


df_comb <- df2%>%full_join(df3, by="ind")%>%drop_na(y_diff_gene, y_diff_motif)

corr <- cor.test(df_comb$y_diff_gene, df_comb$y_diff_motif, method="spearman")
rr <- round(as.numeric(corr$estimate), digits=3)
eq <- deparse(bquote(italic(rho)==.(rr)~",NS")) 


p3 <- ggplot(df_comb, aes(y_diff_gene, y_diff_motif))+
   geom_point(size=3, color="#dd1c77", shape=1)+
   geom_smooth(se=F, method=lm)+
   annotate("text", label=eq, x=0, y=0.55, size=3, parse=T)+ 
   xlab("Relative Expression to ctrl")+
   ylab("Relative Motif to ctrl")+
   theme_bw()+
   theme(axis.title=element_text(size=10),
         axis.text=element_text(size=10))
##
figfn <- paste(outdir2, "Figure1.3_", gene0, "_", ii, ".scatter.png", sep="")
ggsave(figfn, p3,  width=380, height=380, units="px", dpi=120)





###
### another example 



###
### gene expression 

gene0 <- "IRF7" 
ii <- "2_NKcell_caffeine"
    
plotDF <- getData(mat_gene, gene0, comb=ii)    
p1 <- ggplot(plotDF, aes(x=factor(treats), y=y, color=treats))+
   geom_boxplot(outlier.shape=NA, lwd=0.3)+
   geom_jitter(width=0.2, size=1)+
   stat_summary(fun=median, geom="point", shape=23, size=5, stroke=0.9)+ 
   scale_y_continuous("Normalized expression",
      expand=expansion(mult=c(0.2, 0.2)))+
   scale_color_manual("", values=col2)+ 
   theme_bw()+
   theme(## axis.text.x=element_text(angle=-90, hjust=0, size=10),
         axis.text=element_text(size=10),
         axis.title.x=element_blank(),
         axis.title.y=element_text(size=12),
         legend.position="none")
  
figfn <- paste(outdir2, "Figure2.1_", gene0, "_", ii, "_gene.box.png", sep="")
ggsave(figfn, p1,  width=480, height=360, units="px", dpi=120)



###
### TF motifs 
 
gene0 <- "IRF7"
m0_id <- DF_motif%>%filter(motif_name==gene0)%>%pull(motif_ID)
ii <- "2_NKcell_caffeine"
  
###
plotDF2 <- getData(mat_motif, m0_id, comb=ii)
p2 <- ggplot(plotDF2, aes(x=factor(treats), y=y, color=treats))+
   geom_boxplot(outlier.shape=NA, lwd=0.3)+
   geom_jitter(width=0.2, size=1)+
   stat_summary(fun=median, geom="point", shape=23, size=5, stroke=0.9)+ 
   scale_y_continuous("TF activity",
      expand=expansion(mult=c(0.2, 0.2)))+
   scale_color_manual("", values=col2)+ 
   theme_bw()+
   theme(## axis.text.x=element_text(angle=-90, hjust=0, size=10),
         axis.text=element_text(size=10),
         axis.title.x=element_blank(),
         axis.title.y=element_text(size=12),
         legend.position="none")
  
figfn <- paste(outdir2, "Figure2.2_", gene0, "_", ii, "_TF.box.png", sep="")
ggsave(figfn, p2,  width=420, height=360, units="px", dpi=120)



###
### scatter plots

trt <- gsub(".*_", "", ii)

## gene
df0 <- plotDF%>%filter(treats==trt)%>%dplyr::select(ind, y2=y)
df1 <- plotDF%>%filter(treats=="ctrl")%>%dplyr::select(ind, y0=y)

df2 <- df0%>%full_join(df1, by="ind")%>%mutate(y_diff=y2-y0)%>%
   dplyr::select(ind, y_diff_gene=y_diff) 


### motif
df0 <- plotDF2%>%filter(treats==trt)%>%dplyr::select(ind, y2=y)
df1 <- plotDF2%>%filter(treats=="ctrl")%>%dplyr::select(ind, y0=y)

df3 <- df0%>%full_join(df1, by="ind")%>%mutate(y_diff=y2-y0)%>%
   dplyr::select(ind, y_diff_motif=y_diff) 


df_comb <- df2%>%full_join(df3, by="ind")%>%drop_na(y_diff_gene, y_diff_motif)

corr <- cor.test(df_comb$y_diff_gene, df_comb$y_diff_motif, method="pearson")
rr <- round(as.numeric(corr$estimate), digits=3)
eq <- deparse(bquote(italic(rho)==.(rr)~",NS")) 
 

p3 <- ggplot(df_comb, aes(y_diff_gene, y_diff_motif))+
   geom_point(size=3, color="#dd1c77", shape=1)+
   geom_smooth(se=F, method=lm)+
   annotate("text", label=eq, x=0, y=0.55, size=3, parse=T)+ 
   xlab("Relative Expression to ctrl")+
   ylab("Relative Motif to ctrl")+
   theme_bw()+
   theme(axis.title=element_text(size=10),
         axis.text=element_text(size=10))
##
figfn <- paste(outdir2, "Figure2.3_", gene0, "_", ii, ".scatter.png", sep="")
ggsave(figfn, p3,  width=380, height=380, units="px", dpi=120)
