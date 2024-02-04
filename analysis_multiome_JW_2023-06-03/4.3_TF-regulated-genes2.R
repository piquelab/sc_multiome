##
library(Matrix)
library(tidyverse)
library(data.table)
library(expm)
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
library("rainbow")  ##, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages")
library("fds") ##, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages")
library(fda)   ##, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages")
library(CCA) ##, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages")

library("GGally") 


rm(list=ls())





####################################################################
### Heatmap-1, correlation, TF motif activity vs gene expression ###
####################################################################

### input data
## res <- x%>%dplyr::select(comb, MCls, contrast, motif_name, beta_x, beta_y=beta_rna)%>%drop_na(beta_x, beta_y)


## ###
## ### correlation mat
## dfcorr <- res%>%group_by(MCls, contrast)%>%
##     summarize(rr=cor(beta_x, beta_y, method="spearman"), .groups="drop")%>%as.data.frame()
## dfmat <- dfcorr%>%pivot_wider(id_cols=MCls, names_from=contrast, values_from=rr)%>%as.data.frame()

## mat <- dfmat%>%column_to_rownames(var="MCls")
## mat <- as.matrix(mat)



## ###
## ### significance
## sig <- res%>%group_by(MCls, contrast)%>%
##     summarize(pval=cor.test(beta_x, beta_y, method="spearman")$p.value, .groups="drop")%>%as.data.frame()
## sig <- sig%>%mutate(FDR=p.adjust(pval, method="BH"), is_sig=ifelse(FDR<0.05, 1, 0))

## sig_mat <- sig%>%pivot_wider(id_cols=MCls, names_from=contrast, values_from=is_sig)%>%as.data.frame()
 
## sig_mat <- sig_mat%>%column_to_rownames(var="MCls")
## sig_mat <- as.matrix(sig_mat)


## mat2 <- mat*sig_mat
## mat2[mat2==0] <- NA


## ### setting color
## ## mybreak <- seq(0, 1, length.out=20)
## ## col0 <- brewer.pal(n=9,name="Reds")
## ## mycol <- colorRamp2(mybreak, colorRampPalette(col0)(20))

## mycol <- colorRamp2(seq(-1, 1, length.out=50), colorRampPalette(rev(brewer.pal(n=7,name="RdBu")))(50))


## p <- Heatmap(mat2, name="SCC", na_col="grey90", 
##     col=mycol, cluster_rows=F, cluster_columns=F,
##     row_names_gp=gpar(fontsize=10), column_names_gp=gpar(fontsize=10), ##column_names_rot=-45,
##     heatmap_legend_param=list(at=seq(-1, 1, by=0.5),
##         grid_width=grid::unit(0.38, "cm"), legend_height=grid::unit(4, "cm"),
##         title_gp=gpar(fontsize=8), labels_gp=gpar(fontsize=8)),
##     cell_fun=function(j, i, x, y, width, height, fill){
##        grid.text(round(mat[i,j],digits=3), x, y, gp=gpar(fontsize=9))
##     })

## figfn <- paste(outdir2, "Figure4.0_corr_DE_Motif_heatmap.png", sep="")
## png(figfn, height=500, width=520, res=120)
## print(p)
## dev.off()



## ###
## ### TF activity vs TF-regulated-genes, z-score

## fn <- "./4_motif_plots.outs/2_all_response/1.0_cl6_gene_reorder.xlsx"
## df_cl <- read.xlsx(fn)
## motifSel <- unique(df_cl$gene)

## ### input data
## x2 <- x%>%dplyr::filter(motif_name%in%motifSel)
## res <- x2%>%dplyr::select(comb, MCls, contrast, motif_name, x=zscore_x, y=zscore_rna)%>%drop_na(x, y)


## ###
## ### correlation mat
## dfcorr <- res%>%group_by(MCls, contrast)%>%
##     summarize(rr=cor(x, y, method="spearman"), .groups="drop")%>%as.data.frame()
## dfmat <- dfcorr%>%pivot_wider(id_cols=MCls, names_from=contrast, values_from=rr)%>%as.data.frame()

## mat <- dfmat%>%column_to_rownames(var="MCls")
## mat <- as.matrix(mat)



## ###
## ### significance
## sig <- res%>%group_by(MCls, contrast)%>%
##     summarize(pval=cor.test(x, y, method="spearman")$p.value, .groups="drop")%>%as.data.frame()
## sig <- sig%>%mutate(FDR=p.adjust(pval, method="BH"), is_sig=ifelse(FDR<0.05, 1, 0))

## sig_mat <- sig%>%pivot_wider(id_cols=MCls, names_from=contrast, values_from=is_sig)%>%as.data.frame()
 
## sig_mat <- sig_mat%>%column_to_rownames(var="MCls")
## sig_mat <- as.matrix(sig_mat)


## mat2 <- mat*sig_mat
## mat2[mat2==0|mat2<0] <- NA


## ### setting color
## mybreak <- seq(0, 1, length.out=20)
## col0 <- brewer.pal(n=9,name="Reds")
## mycol <- colorRamp2(mybreak, colorRampPalette(col0)(20))


## ## mycol <- colorRamp2(seq(0, 1, length.out=12), colorRampPalette(c("white", "red"))(12))

## p <- Heatmap(mat2, name="SCC", na_col="grey90", 
##     col=mycol, cluster_rows=F, cluster_columns=F,
##     row_names_gp=gpar(fontsize=10), column_names_gp=gpar(fontsize=10), ##column_names_rot=-45,
##     heatmap_legend_param=list(at=seq(0, 1, by=0.25),
##         grid_width=grid::unit(0.38, "cm"), legend_height=grid::unit(4, "cm"),
##         title_gp=gpar(fontsize=8), labels_gp=gpar(fontsize=8)),
##     cell_fun=function(j, i, x, y, width, height, fill){
##        grid.text(round(mat[i,j],digits=3), x, y, gp=gpar(fontsize=9))
##     })

## figfn <- paste(outdir2, "Figure4.0_corr3_zscore.heatmap.png", sep="")
## png(figfn, height=500, width=520, res=120)
## print(p)
## dev.off()





#####################################################
### TF activity and TF-regulated gene expression ####
#####################################################

## rm(list=ls())

## outdir2 <- "./4_motif_plots.outs/2_all_response/"

## ###
## fn <- paste(outdir2, "1.0_cl6_gene_reorder.xlsx", sep="")
## df_cl <- read.xlsx(fn)
## resp_motif <- unique(df_cl$gene)

## ##
## fn <- paste(outdir2, "4_LFC.TF_LFC.DE.rds", sep="")
## res <- read_rds(fn)

## x2 <- res%>%filter(comb=="4_Bcell_caffeine")%>%dplyr::select(motif_name, x1=beta_x, x2=beta_rna) 
## df2 <- 

## res2 <- res%>%group_by(comb, motif_name)%>%slice_max(order_by=abs(zscore_x), n=1)%>%as.data.frame()
 
## dfmat <- res2%>%pivot_wider(id_cols=motif_name, names_from=comb, values_from=zscore_rna, values_fill=NA)
## mat <- dfmat%>%column_to_rownames(var="motif_name")%>%as.matrix()


## mat <- mat[resp_motif,]
## iSel <- rowSums(is.na(mat))==0
## mat <- mat[iSel,]

## ###
## ### get colnames and re-order by treats
## rn <- colnames(mat)
## x <- str_split(rn, "_", simplify=T)
## cvt <- data.frame(comb=rn, MCls=paste(x[,1], x[,2], sep="_"), contrast=x[,3])
## cvt <- cvt%>%arrange(contrast)

## mat2 <- mat[,cvt$comb]


## ###
## ### color for heatmap value
## y <- as.vector(mat2)
## ## y0 <- y[abs(y)<2]
## quantile(abs(y), probs=0.99)

## mybreak <- c(min(y,na.rm=T), seq(-6, 6, length.out=98), max(y,na.rm=T))

## ## quantile(abs(y), probs=c(0.9,0.95,0.99))
## ## range(y)

## mycol <- colorRamp2(mybreak, colorRampPalette(rev(brewer.pal(n=7,name="RdBu")))(100))

 
## ###
## ### annotation columns
## col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
##   "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
##    "6_Monocyte"="#984ea3", "7_dnT"="black")
## col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
##        "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3")


## x <- str_split(colnames(mat2), "_", simplify=T)
## df_col <- data.frame(celltype=paste(x[,1], x[,2], sep="_"), contrast=x[,3])
## col_ha <- HeatmapAnnotation(df=df_col, col=list(celltype=col1, contrast=col2),
##     annotation_name_gp=gpar(fontsize=10),
##     annotation_legend_param=list(
##   celltype=list(title_gp=gpar(fontsize=9), labels_gp=gpar(fontsize=9), title="celltype",
##                 grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm")),
##   contrast=list(title_gp=gpar(fontsize=9), labels_gp=gpar(fontsize=9), title="contrast",
##                 grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm"))))




## ###
## ### row annotation
## df_cl2 <- df_cl%>%filter(gene%in%rownames(mat))
## identical(df_cl2$gene, rownames(mat))


## ## df_row <- df_cl2%>%dplyr::select(cluster)
## ## row_ha <- rowAnnotation(cl=anno_block(gp=gpar(fill=2:7), labels=c("CL2", "CL1", "CL3", "CL4", "CL5", "CL6")))
 
## ## table(DF$cluster)
## ##mat[is.na(mat)] <- 0
## p0 <- Heatmap(mat2, col=mycol, 
##    cluster_rows=T, 
##    cluster_columns=F,
##    ## row_split=c(rep("1",32), rep("2", 8), rep("3", 21), rep("4", 14), rep("5",12), rep("6", 20)),
##    ##row_split=c(rep("1", 9), rep("2", 17), rep("3", 38)), 
##    show_row_names=T, row_names_gp=gpar(fontsize=5),
##    show_column_names=T, column_names_gp=gpar(fontsize=7),
##    show_row_dend=F, show_column_dend=F,
##    row_title=NULL,
##    ##column_names_rot=-45,   
##    top_annotation=col_ha,
##    ## left_annotation=row_ha,
##    heatmap_legend_param=list(title=bquote(~italic(Z)~"-score"),
##       title_gp=gpar(fontsize=9),
##       at=seq(-9, 9, by=3), 
##       labels_gp=gpar(fontsize=9),
##       grid_width=grid::unit(0.45, "cm"),
##       legend_height=grid::unit(7.8, "cm")))

## ###
## figfn <- paste(outdir2, "Figure4.1_TF-regulated-gene.heatmap.png", sep="")
## png(figfn, width=1000, height=1050,res=120)
## set.seed(0)
## p0 <- draw(p0)
## dev.off()


 

#########################################################################################################
### generate TF-regulated-genes data as a matrix form used for plots and secondary analysis 
#########################################################################################################


### ./4_motif_plots.outs/TF-regulated-genes2/
### 0_TF_peaks.R generate TF-peaks
### 1_TF-regulated-genes.R generate beta and z-score of TF-regulated-genes for all response motifs
### 4.0_TF_peaks.rds, TF-peaks for response motifs
### 4_response_TF-regulated-genes.rds for all response motifs
### 4_response_TF-regulated-genes.mat.rds for the top 10 response motifs (102)


rm(list=ls())
outdir2 <- "./4_motif_plots.outs/TF-regulated-genes2/"

###
### top 10 motifs

fn <- "./4_motif_plots.outs/3.2_motif_diff_fitInd_JW.rds"
res <- read_rds(fn)%>%
    mutate(comb=paste(MCls, contrast, sep="_"), zscore=beta/stderr,
                      is_sig=ifelse(qval<0.1, 1, 0))

### response motifs
fn <- "./4_motif_plots.outs/response_motif_correct/2.2_response_motif_th0.1.txt"
resp_motif <- read.table(fn, header=TRUE)$motif_name%>%unique()


### top 10
res2 <- res%>%filter(motif_name%in%resp_motif)
top10 <- res2%>%group_by(comb)%>%slice_max(order_by=abs(zscore), n=10)%>%pull(motif_name)%>%unique()

 


###
### TF activity and TF-regulated genes
fn <- paste(outdir2, "4_response_TF-regulated-genes.rds", sep="")
df2 <- read_rds(fn)
df2 <- df2%>%group_by(comb, motif_name)%>%slice_max(order_by=abs(zscore_x), n=1)%>%as.data.frame()


###
x1 <- df2%>%dplyr::select(comb, motif_name, zscore_x)%>%
   pivot_wider(id_cols=motif_name, names_from=comb, values_from=zscore_x, values_fill=NA)%>%
   column_to_rownames(var="motif_name")%>%
   as.matrix()

###
x2 <- df2%>%dplyr::select(comb, motif_name, zscore_rna)%>%
   pivot_wider(id_cols=motif_name, names_from=comb, values_from=zscore_rna, values_fill=NA)%>%
   column_to_rownames(var="motif_name")%>%
   as.matrix()

identical(colnames(x1), colnames(x2))
identical(rownames(x1), rownames(x2))


###
### rename x1
comb2 <- paste0("1_", colnames(x1), sep="")
colnames(x1) <- comb2
x1 <- x1[top10,]

###
### rename x2
comb2 <- paste0("2_", colnames(x2), sep="")
colnames(x2) <- comb2
x2 <- x2[top10,]

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

opfn <- paste(outdir2, "4_response_TF-regulated-genes.mat.rds", sep="")
write_rds(mat2, file=opfn)



######################################################################################
### correlation values between TF activity and TF-regulated-genes 
######################################################################################

rm(list=ls())
outdir2 <- "./4_motif_plots.outs/TF-regulated-genes2/"

###
fn <- paste(outdir2, "4_response_TF-regulated-genes.rds", sep="")
x <- read_rds(fn)
x <- x%>%group_by(comb, motif_name)%>%slice_max(order_by=abs(zscore_x), n=1)%>%as.data.frame()

### beta correlation

## motif changes-beta
df1 <- x%>%dplyr::select(comb, motif_name, b=beta_x)
mat1 <- df1%>%pivot_wider(id_cols=motif_name, names_from=comb, values_from=b, values_fill=NA)%>%
    column_to_rownames(var="motif_name")
id1 <- rownames(mat1)
comb1 <- names(mat1)



### TF gene expression changes
df2 <- x%>%dplyr::select(comb, motif_name, b=beta_rna)
mat2 <- df2%>%pivot_wider(id_cols=motif_name, names_from=comb, values_from=b, values_fill=NA)%>%
    column_to_rownames(var="motif_name")
id2 <- rownames(mat2)
comb2 <- names(mat2)

identical(id1, id2)
identical(comb1, comb2)



id_sel <- rowSums(is.na(cbind(mat1, mat2)))==0


###
###
mat1 <- mat1[id_sel,]
mat2 <- mat2[id_sel,]
nm <- nrow(mat1)
id <- rownames(mat1)

dfcorr2 <- map_dfr(1:nm, function(i){
   ##
   id0 <- id[i] 
   b1 <- as.numeric(mat1[i,])
   b2 <- as.numeric(mat2[i,])
   rr <- cor(b1, b2, method="spearman")
   df0 <- data.frame(motif_name=id0, "rr_beta"=rr)
   df0 
})


###################################
### z-score correlation 
#########################################


## motif changes-zscore
df1 <- x%>%dplyr::select(comb, motif_name, b=zscore_x)
mat1 <- df1%>%pivot_wider(id_cols=motif_name, names_from=comb, values_from=b, values_fill=NA)%>%
    column_to_rownames(var="motif_name")
id1 <- rownames(mat1)
comb1 <- names(mat1)



### TF-regulated-genes expression changes
df2 <- x%>%dplyr::select(comb, motif_name, b=zscore_rna)
mat2 <- df2%>%pivot_wider(id_cols=motif_name, names_from=comb, values_from=b, values_fill=NA)%>%
    column_to_rownames(var="motif_name")
id2 <- rownames(mat2)
comb2 <- names(mat2)

identical(id1, id2)
identical(comb1, comb2)


id_sel <- rowSums(is.na(cbind(mat1, mat2)))==0



###
###
mat1 <- mat1[id_sel,]
mat2 <- mat2[id_sel,]
nm <- nrow(mat1)
id <- rownames(mat1)

dfcorr3 <- map_dfr(1:nm, function(i){
   ##
   id0 <- id[i] 
   b1 <- as.numeric(mat1[i,])
   b2 <- as.numeric(mat2[i,])
   rr <- cor(b1, b2, method="spearman")
   df0 <- data.frame(motif_name=id0, "rr_zscore"=rr)
   df0 
})


dfcomb <- dfcorr2%>%left_join(dfcorr3, by="motif_name")%>%arrange(desc(rr_zscore))
opfn <- paste(outdir2, "4.2_corr.xlsx", sep="")
write.xlsx(dfcomb, file=opfn)



########################################
### heatmap of correlation 
#########################################


rm(list=ls())
outdir2 <- "./4_motif_plots.outs/TF-regulated-genes2/"

###
fn <- paste(outdir2, "4_response_TF-regulated-genes.rds", sep="")
x <- read_rds(fn)
x <- x%>%group_by(comb, motif_name)%>%slice_max(order_by=abs(zscore_x), n=1)%>%as.data.frame()

### motif select
fn <- paste(outdir2, "4.2_corr.xlsx", sep="")
motifSel <- read.xlsx(fn)%>%pull(motif_name)


### input data
x2 <- x%>%dplyr::filter(motif_name%in%motifSel)
res <- x2%>%dplyr::select(comb, MCls, contrast, motif_name, x=zscore_x, y=zscore_rna)


###
### correlation mat
dfcorr <- res%>%group_by(MCls, contrast)%>%
    summarize(rr=cor(x, y, method="spearman"), .groups="drop")%>%as.data.frame()
dfmat <- dfcorr%>%pivot_wider(id_cols=MCls, names_from=contrast, values_from=rr)%>%as.data.frame()

mat <- dfmat%>%column_to_rownames(var="MCls")%>%as.matrix()



###
### significance
sig <- res%>%group_by(MCls, contrast)%>%
    summarize(pval=cor.test(x, y, method="spearman")$p.value, .groups="drop")%>%as.data.frame()
sig <- sig%>%mutate(FDR=p.adjust(pval, method="BH"), is_sig=ifelse(FDR<0.05, 1, 0))

sig_mat <- sig%>%pivot_wider(id_cols=MCls, names_from=contrast, values_from=is_sig)%>%as.data.frame() 
sig_mat <- sig_mat%>%column_to_rownames(var="MCls")%>%as.matrix()


mat2 <- mat*sig_mat
mat2[mat2==0] <- NA


### setting color
## mybreak <- seq(0, 1, length.out=20)
## col0 <- brewer.pal(n=9,name="Reds")
## mycol <- colorRamp2(mybreak, colorRampPalette(col0)(20))

mycol <- colorRamp2(seq(-1, 1, length.out=50), colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(50))


## mycol <- colorRamp2(seq(0, 1, length.out=12), colorRampPalette(c("white", "red"))(12))

p <- Heatmap(mat2, name="SCC", na_col="grey90", 
    col=mycol, cluster_rows=F, cluster_columns=F,
    row_names_gp=gpar(fontsize=10), column_names_gp=gpar(fontsize=10), ##column_names_rot=-45,
    heatmap_legend_param=list(at=seq(-1, 1, by=0.5),
        grid_width=grid::unit(0.38, "cm"), legend_height=grid::unit(4, "cm"),
        title_gp=gpar(fontsize=8), labels_gp=gpar(fontsize=8)),
    cell_fun=function(j, i, x, y, width, height, fill){
       grid.text(round(mat[i,j],digits=3), x, y, gp=gpar(fontsize=9))
    })

figfn <- paste(outdir2, "summary_results/Figure0_corr_zscore.heatmap.png", sep="")
png(figfn, height=500, width=520, res=120)
print(p)
dev.off()




###
### old analysis 

## rm(list=ls())
## outdir2 <- "./4_motif_plots.outs/2_all_response/"

## ###
## fn <- paste(outdir2, "4_LFC.TF_LFC.DE.rds", sep="")
## x <- read_rds(fn)
## x <- x%>%group_by(comb, motif_name)%>%slice_max(order_by=abs(zscore_x), n=1)%>%as.data.frame()

## mat2 <- x%>%pivot_wider(id_cols=motif_name, names_from=comb, values_from=zscore_rna, values_fill=NA)%>%
##     column_to_rownames(var="motif_name")%>%as.matrix()
## motif <- rownames(mat2)
## idSel <- rowSums(is.na(mat2))==0
## motifSel <- motif[idSel]    


## x2 <- x%>%dplyr::filter(motif_name%in%motifSel)
## res <- x2%>%dplyr::select(comb, MCls, contrast, motif_name, x=zscore_x, y=zscore_rna)


## ###
## ### correlation mat
## dfcorr <- res%>%group_by(MCls, contrast)%>%
##     summarize(rr=cor(x, y, method="spearman"), .groups="drop")%>%as.data.frame()
## dfmat <- dfcorr%>%pivot_wider(id_cols=MCls, names_from=contrast, values_from=rr)%>%as.data.frame()

## mat <- dfmat%>%column_to_rownames(var="MCls")%>%as.matrix()



## ###
## ### significance
## sig <- res%>%group_by(MCls, contrast)%>%
##     summarize(pval=cor.test(x, y, method="spearman")$p.value, .groups="drop")%>%as.data.frame()
## sig <- sig%>%mutate(FDR=p.adjust(pval, method="BH"), is_sig=ifelse(FDR<0.05, 1, 0))

## sig_mat <- sig%>%pivot_wider(id_cols=MCls, names_from=contrast, values_from=is_sig)%>%as.data.frame() 
## sig_mat <- sig_mat%>%column_to_rownames(var="MCls")%>%as.matrix()


## mat2 <- mat*sig_mat
## mat2[mat2==0] <- NA


## ### setting color
## ## mybreak <- seq(0, 1, length.out=20)
## ## col0 <- brewer.pal(n=9,name="Reds")
## ## mycol <- colorRamp2(mybreak, colorRampPalette(col0)(20))

## mycol <- colorRamp2(seq(-1, 1, length.out=50), colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(50))


## ## mycol <- colorRamp2(seq(0, 1, length.out=12), colorRampPalette(c("white", "red"))(12))

## p <- Heatmap(mat2, name="SCC", na_col="grey90", 
##     col=mycol, cluster_rows=F, cluster_columns=F,
##     row_names_gp=gpar(fontsize=10), column_names_gp=gpar(fontsize=10), ##column_names_rot=-45,
##     heatmap_legend_param=list(at=seq(-1, 1, by=0.5),
##         grid_width=grid::unit(0.38, "cm"), legend_height=grid::unit(4, "cm"),
##         title_gp=gpar(fontsize=8), labels_gp=gpar(fontsize=8)),
##     cell_fun=function(j, i, x, y, width, height, fill){
##        grid.text(round(mat[i,j],digits=3), x, y, gp=gpar(fontsize=9))
##     })

## figfn <- paste(outdir2, "Figure4.0_corr_zscore.heatmap.png", sep="")
## png(figfn, height=500, width=520, res=120)
## print(p)
## dev.off()






#############################
### heatmap plots 
##################################

 
rm(list=ls())

outdir2 <- "./4_motif_plots.outs/TF-regulated-genes2/"
if ( !file.exists(outdir2) ) dir.create(outdir2, showWarnings=F, recursive=T)

fn <- paste(outdir2, "4_response_TF-regulated-genes.mat.rds", sep="")
mat <- read_rds(fn)

###
### setting color for heatmap value
y <- as.vector(mat)
##y0 <- y[abs(y)<2] ### quantile(abs(y), probs=0.99)
mybreak <- c(min(y), seq(-6, 6, length.out=98), max(y))
## mybreak <- c(min(y), seq(-1.5, 1.5, length.out=98), max(y))
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


df_col <- x2%>%dplyr::select(MCls, contrast) ##%>%mutate(Feature=as.character(Feature))

col_ha <- HeatmapAnnotation(df=df_col, col=list(MCls=col1, contrast=col2),
    annotation_name_gp=gpar(fontsize=10),
    annotation_legend_param=list(       
  MCls=list(title_gp=gpar(fontsize=9), labels_gp=gpar(fontsize=9), title="celltype",
                grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm")),
  contrast=list(title_gp=gpar(fontsize=9), labels_gp=gpar(fontsize=9), title="contrast",
                grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm"))))

mat2 <- mat[,x2$comb]

comb <- gsub("^[12]_", "", colnames(mat2))
colnames(mat2) <- comb


## show_legend=c(F,F))
### main plots
p2 <- Heatmap(mat2, col=mycol,
   cluster_rows=T, cluster_columns=F, row_km=6, column_split=rep(c("TF activity", "TF regulated genes"), each=48),
   show_row_names=T, row_names_gp=gpar(fontsize=5),
   show_column_names=T, column_names_gp=gpar(fontsize=8),
   show_row_dend=T, row_dend_width=grid::unit(2, "cm"),
   show_column_dend=F,
   ##column_names_rot=-45,   
   top_annotation=col_ha,
   heatmap_legend_param=list(title=bquote(~italic(Z)~"-score"),
      title_gp=gpar(fontsize=9),
      at=seq(-9, 9, by=3), 
      labels_gp=gpar(fontsize=9),
      grid_width=grid::unit(0.45, "cm"),
      legend_height=grid::unit(7.8, "cm")))
 
###
figfn <- paste(outdir2, "summary_results/Figure1_top10_comb2.heatmap.png", sep="")
## ggsave(figfn, plot=p2, device="png", width=680, height=700, units="px", dpi=300) ## ggsave 
png(figfn, width=1600, height=1200,res=120)
set.seed(0)
p2 <- draw(p2)
dev.off()


###
### TF clusters
hmap <- Heatmap(mat2, cluster_rows=T, row_km=6, cluster_columns=F)
set.seed(0)
hmap <- draw(hmap)
cl <- row_order(hmap)
lengths(cl)
 
geneSel <- rownames(mat2)
DF_cl <- NULL
for (i in names(cl)){
   cl_tmp <- data.frame(cluster=i, motif_name=geneSel[cl[[i]]])
   DF_cl <- rbind(DF_cl, cl_tmp)
}

## newCL <- c("3"="1", "2"="2", "1"="3")
## DF_cl <- DF_cl%>%mutate(cluster2=newCL[cluster])
## DF2 <- DF_cl%>%arrange(cluster2)

fn <- paste(outdir2, "4.2_corr.xlsx", sep="")
DFcorr <- read.xlsx(fn)
DFcomb <- DF_cl%>%left_join(DFcorr, by="motif_name")

opfn <- paste(outdir2, "4.2_corr_reorder.xlsx", sep="")
write.xlsx(DFcomb, file=opfn)

 


###############################
### Boxplots of correlation ###
###############################

rm(list=ls())

outdir2 <- "./4_motif_plots.outs/TF-regulated-genes2/"
if ( !file.exists(outdir2) ) dir.create(outdir2, showWarnings=F, recursive=T)


fn <- paste(outdir2, "4.2_corr_reorder.xlsx", sep="")
df2 <- read.xlsx(fn)%>%mutate(cluster2=paste("CL", cluster, sep=""))
###

p <- ggplot(df2, aes(x=cluster2, y=rr_zscore, color=cluster2))+
   geom_boxplot(outlier.shape=NA, lwd=0.25)+
   geom_jitter(width=0.2, size=0.5)+ 
   stat_summary(fun=median, geom="point", shape=23, size=1.8, stroke=0.9)+ylim(-1,1)+
   ggtitle("SCC of z-score between TFs activity and TF-regulated genes")+  
   theme_bw()+
   theme(legend.position="none",
         axis.title=element_blank(),
         axis.text=element_text(size=10),
         plot.title=element_text(size=9, hjust=0.5))
 
figfn <- paste(outdir2, "summary_results/Figure1.2_corr.box.png", sep="")
ggsave(figfn, p, width=520, height=350, units="px", dpi=120)



###
### old
## rm(list=ls())

## outdir2 <- "./4_motif_plots.outs/2_all_response/"
## if ( !file.exists(outdir2) ) dir.create(outdir2, showWarnings=F, recursive=T)

## fn <- paste(outdir2, "4_TF_activity_TF-regulated-gene.mat.rds", sep="")
## x <- read_rds(fn)
## x1 <- x[,1:48]
## x2 <- x[,49:96]
## motifs <- rownames(x)

## dfcorr <- map_dfr(motifs, function(ii){
##     ##
##     x1 <- as.numeric(x[ii, 1:48])
##     x2 <- as.numeric(x[ii, 49:96])
##     df0 <- data.frame(motif_name=ii, rr_zscore=cor(x1, x2, method="spearman"))
##     df0
## })    

## fn <- paste(outdir2, "4_comb2_cl6_motif_reorder.xlsx", sep="")
## df_cl <- read.xlsx(fn)%>%rename("motif_name"="gene")

## dfcomb <- df_cl%>%left_join(dfcorr, by="motif_name")
## opfn <- paste(outdir2, "4.2_corr_zscore.xlsx", sep="")
## write.xlsx(dfcomb, file=opfn)


## ####
## ####
 
## rm(list=ls())

## outdir2 <- "./4_motif_plots.outs/2_all_response/"
## fn <- paste(outdir2, "4.2_corr_zscore.xlsx", sep="") 
## df2 <- read.xlsx(fn)%>%mutate(cluster2=paste("CL", cluster, sep=""))


## p <- ggplot(df2, aes(x=cluster2, y=rr_zscore, color=cluster2))+
##    geom_boxplot(outlier.shape=NA, lwd=0.25)+
##    geom_jitter(width=0.2, size=0.5)+ 
##    stat_summary(fun=median, geom="point", shape=23, size=1.8, stroke=0.9)+ylim(-1,1)+
##    ggtitle("SCC of z-score between TFs activity and TF-regulated genes")+  
##    theme_bw()+
##    theme(legend.position="none",
##          axis.title=element_blank(),
##          axis.text=element_text(size=10),
##          plot.title=element_text(size=9, hjust=0.5))
 
## figfn <- paste(outdir2, "Figure4.2_corr.box.png", sep="")
## ggsave(figfn, p, width=520, height=350, units="px", dpi=120)s






######################
#### CCA analysis ####
######################

 
rm(list=ls())
outdir2 <- "./4_motif_plots.outs/2_all_response/"


##
fn <- paste(outdir2, "1.0_cl6_gene_reorder.xlsx", sep="")
df_cl <- read.xlsx(fn)
motifSel <- unique(df_cl$gene)

###
fn <- paste(outdir2, "4_LFC.TF_LFC.DE.rds", sep="")
x <- read_rds(fn)
x <- x%>%group_by(comb, motif_name)%>%slice_max(order_by=abs(zscore_x), n=1)%>%as.data.frame()


###
### motif changes
df1 <- x%>%dplyr::select(comb, motif_name, beta_x)
mat1 <- df1%>%pivot_wider(id_cols=motif_name, names_from=comb, values_from=beta_x, values_fill=NA)%>%
    column_to_rownames(var="motif_name")
##
mat1 <- mat1[motifSel,]
id1 <- rownames(mat1)
comb1 <- names(mat1)

id1_sel <- rowSums(is.na(mat1))==0


###
### TF regulated gene expression changes
df2 <- x%>%dplyr::select(comb, motif_name, beta_rna)
mat2 <- df2%>%pivot_wider(id_cols=motif_name, names_from=comb, values_from=beta_rna, values_fill=NA)%>%
    column_to_rownames(var="motif_name")
###
mat2 <- mat2[motifSel,]
id2 <- rownames(mat2)
comb2 <- names(mat2)

id2_sel <- rowSums(is.na(mat2))==0

id_sel <- id1_sel&id2_sel

identical(id1, id2)
identical(comb1, comb2)


####################
### cca analysis ###
####################



rm(list=ls())
outdir2 <- "./4_motif_plots.outs/TF-regulated-genes2/"


###
fn <- paste(outdir2, "4_TF_activity_TF-regulated-gene.mat.rds", sep="")
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
names(motif) <- paste("s", 1:48, sep="")
###
rna <- as.data.frame(cc2$yscores)
names(rna) <- paste("s", 1:48, sep="")


identical(rownames(motif), rownames(rna))




###
### custom function 
my_box <- function(data, mapping,...){
   ##
   p <- ggplot(data=data, mapping=mapping)+
       geom_boxplot(outlier.shape=NA, lwd=0.25)+coord_flip()
   ## geom_jitter(data=data, mapping=mapping, width=0.2, size=0.5, shape=23)+coord_flip()## + 
   ## stat_summary(fun=median, geom="point", shape=23, size=1.8, stroke=0.9) ###ylim(-1,1)+
   p
}    


fn <- paste(outdir2, "1.0_cl6_gene_reorder.xlsx", sep="")
df_cl <- read.xlsx(fn)




###
### scatter plots of correlation 
df2 <- data.frame(x=1:48, rr=cc0$cor)
###
p0 <- ggplot(df2, aes(x=x, y=rr))+
   geom_point(size=1.5, shape=1, color="#dd1c77")+
   guides(color=guide_legend(override.aes=list(size=2)))+ 
   xlab("Canonical variable")+ylab("Correlation (TF motif-regulated genes)")+
   theme_bw()+
   theme(axis.title=element_text(size=9),
         axis.text=element_text(size=9))
 
figfn <- paste(outdir2, "Figure4.3_corr_score.scatter.png", sep="")
ggsave(figfn, p0, width=420, height=420, units="px", dpi=120)



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
    

figfn <- paste(outdir2, "Figure4.4_motif_score.pair.png", sep="")
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
    

figfn <- paste(outdir2, "Figure4.4_TF-regulated-gene_score.pair.png", sep="")
ggsave(figfn, plot=p2, width=780, height=720, units="px", dpi=120)

 



## ###
## ###
## p <- ggplot(DFcomb, aes(x=cluster2, y=motif_1, color=cluster2))+
##    geom_boxplot(outlier.shape=NA, lwd=0.25)+
##    geom_jitter(width=0.2, size=0.5)+ 
##    stat_summary(fun=median, geom="point", shape=23, size=1.8, stroke=0.9)+ ###ylim(-1,1)+
##    ggtitle("1st score of TF motif")+  
##    theme_bw()+
##    theme(legend.position="none",
##          axis.title=element_blank(),
##          axis.text=element_text(size=10),
##          plot.title=element_text(size=10,hjust=0.5))

## figfn <- paste(outdir2, "Figure4.3_motif1_box.png", sep="")
## ggsave(figfn, p, width=450, height=350, units="px", dpi=120)
 

## ###
## ###
## p2 <- ggplot(DFcomb, aes(x=cluster2, y=rna_1, color=cluster2))+
##    geom_boxplot(outlier.shape=NA, lwd=0.25)+
##    geom_jitter(width=0.2, size=0.5)+ 
##    stat_summary(fun=median, geom="point", shape=23, size=1.8, stroke=0.9)+ ###ylim(-1,1)+
##    ggtitle("1st score of TF-regulated genes")+  
##    theme_bw()+
##    theme(legend.position="none",
##          axis.title=element_blank(),
##          axis.text=element_text(size=10),
##          plot.title=element_text(size=9,hjust=0.5))

## figfn <- paste(outdir2, "Figure4.4_gene1_box.png", sep="")
## ggsave(figfn, p2, width=450, height=350, units="px", dpi=120)


## ###
## ###
## p3 <- ggplot(DFcomb, aes(x=motif_1, y=rna_1, color=cluster2))+
##    geom_point(size=1)+
##    guides(color=guide_legend(override.aes=list(size=2)))+ 
##    xlab("1st score of TF motif")+ylab("1st score of TF-regulated gene")+
##    theme_bw()+
##    theme(axis.title=element_text(size=10),
##          axis.text=element_text(size=10),
##          legend.key.size=grid::unit(0.4, "cm"),
##          legend.title=element_blank())

## figfn <- paste(outdir2, "Figure4.5_1score_scatter.png", sep="")
## ggsave(figfn, p3, width=450, height=320, units="px", dpi=120)




 




#################################
### heatmap show CCA loadings ###
#################################


rm(list=ls())
outdir2 <- "./4_motif_plots.outs/2_all_response/"


###
fn <- paste(outdir2, "4_TF_activity_TF-regulated-gene.mat.rds", sep="")
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



##############################
### Heatmap score from cca ###
##############################

motif <- cc2$xscores
rna <- cc2$yscores

m1 <- motif
colnames(m1) <- paste("S", 1:48, sep="")
m2 <- rna
colnames(m2) <- paste("S", 1:48, sep="")

mat <- cbind(m1, m2)



#### setting colors
y <- as.vector(mat)
##y0 <- y[abs(y)<2] ### quantile(abs(y), probs=0.99)
mybreak <- c(min(y), seq(-3, 3, length.out=98), max(y)) 
mycol <- colorRamp2(mybreak, colorRampPalette(rev(brewer.pal(n=7,name="RdBu")))(100))

quantile(abs(y), probs=c(0.9, 0.95, 0.99))

### main plots
p2 <- Heatmap(mat, col=mycol,
   cluster_rows=T, row_km=6, cluster_columns=F,
   column_split=rep(c("TF activity", "TF-regulated genes"), each=48),
   show_row_dend=T, row_dend_width=grid::unit(2, "cm"),   
   show_row_names=T, row_names_gp=gpar(fontsize=6),
   show_column_names=T, column_names_gp=gpar(fontsize=6),
   ## show_column_names=F, column_names_gp=gpar(fontsize=6),
   ###top_annotation=col_ha,
   heatmap_legend_param=list(title="score",
      title_gp=gpar(fontsize=9),
      at=seq(-4.5, 4.5, by=1.5), 
      labels_gp=gpar(fontsize=9),
      grid_width=grid::unit(0.45, "cm"),
      legend_height=grid::unit(7.8, "cm")))
 
###
figfn <- paste(outdir2, "Figure4.5_score_comb2.heatmap.png", sep="")
## ggsave(figfn, plot=p2, device="png", width=680, height=700, units="px", dpi=300) ## ggsave 
png(figfn, width=1600, height=1200,res=120)
set.seed(0)
p2 <- draw(p2)
dev.off()



##########################
#### loading  from cca ###
##########################


m1 <- t(cc2$corr.X.xscores)
m2 <- t(cc2$corr.Y.yscores)

rownames(m1) <- paste("S", 1:48, sep="")
rownames(m2) <- paste("S", 1:48, sep="")


#### setting colors
y <- as.vector(m1)
##y0 <- y[abs(y)<2] ### quantile(abs(y), probs=0.99)
## mybreak <- c(min(y), seq(-1, 1, length.out=98), max(y))
mybreak <- seq(-1, 1, length.out=100)
mycol <- colorRamp2(mybreak, colorRampPalette(rev(brewer.pal(n=7,name="RdBu")))(100))

## quantile(abs(y), probs=c(0.9, 0.95, 0.99))


###
### annotation columns
col0 <- c("1"="#8c510a", "2"="#d8b365")
col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
  "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
   "6_Monocyte"="#984ea3", "7_dnT"="black")
col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
       "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3")


cvt <- str_split(colnames(m1), "_", simplify=T)
cvt2 <- data.frame(comb=colnames(m1),
   MCls=paste(cvt[,1], cvt[,2], sep="_"), contrast=cvt[,3])
cvt2 <- cvt2%>%arrange(contrast)


df_col <- cvt2%>%dplyr::select(MCls, contrast) ##%>%mutate(Feature=as.character(Feature))

col_ha <- HeatmapAnnotation(df=df_col, col=list(MCls=col1, contrast=col2),
    annotation_name_gp=gpar(fontsize=10),
    annotation_legend_param=list(       
  MCls=list(title_gp=gpar(fontsize=9), labels_gp=gpar(fontsize=9), title="celltype",
                grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm")),
  contrast=list(title_gp=gpar(fontsize=9), labels_gp=gpar(fontsize=9), title="contrast",
                grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm"))))

### plot data
m1a <- m1[,cvt2$comb]


### main plots
p1 <- Heatmap(m1a, col=mycol,
   cluster_rows=F, cluster_columns=F,
   show_column_dend=F, column_dend_height=grid::unit(2, "cm"),   
   show_row_names=T, row_names_gp=gpar(fontsize=7),
   show_column_names=T, column_names_gp=gpar(fontsize=6),
   top_annotation=col_ha,
   heatmap_legend_param=list(title="loading",
      title_gp=gpar(fontsize=9),
      at=seq(-1, 1, by=0.5), 
      labels_gp=gpar(fontsize=9),
      grid_width=grid::unit(0.45, "cm"),
      legend_height=grid::unit(7.8, "cm")))
 
###
figfn <- paste(outdir2, "Figure4.6_loading_motif.heatmap.png", sep="")
## ggsave(figfn, plot=p2, device="png", width=680, height=700, units="px", dpi=300) ## ggsave 
png(figfn, width=720, height=720,res=120)
set.seed(0)
p1 <- draw(p1)
dev.off()



###
### loading of RNA

#### setting colors
y <- as.vector(m2)
##y0 <- y[abs(y)<2] ### quantile(abs(y), probs=0.99)
## mybreak <- c(min(y), seq(-1, 1, length.out=98), max(y))
mybreak <- seq(-1, 1, length.out=100)
mycol <- colorRamp2(mybreak, colorRampPalette(rev(brewer.pal(n=7,name="RdBu")))(100))

## quantile(abs(y), probs=c(0.9, 0.95, 0.99))


###
### annotation columns
col0 <- c("1"="#8c510a", "2"="#d8b365")
col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
  "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
   "6_Monocyte"="#984ea3", "7_dnT"="black")
col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
       "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3")


cvt <- str_split(colnames(m2), "_", simplify=T)
cvt2 <- data.frame(comb=colnames(m2),
   MCls=paste(cvt[,1], cvt[,2], sep="_"), contrast=cvt[,3])
cvt2 <- cvt2%>%arrange(contrast)


df_col <- cvt2%>%dplyr::select(MCls, contrast) ##%>%mutate(Feature=as.character(Feature))

col_ha <- HeatmapAnnotation(df=df_col, col=list(MCls=col1, contrast=col2),
    annotation_name_gp=gpar(fontsize=10),
    annotation_legend_param=list(       
  MCls=list(title_gp=gpar(fontsize=9), labels_gp=gpar(fontsize=9), title="celltype",
                grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm")),
  contrast=list(title_gp=gpar(fontsize=9), labels_gp=gpar(fontsize=9), title="contrast",
                grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm"))))



### plot data
m2a <- m2[,cvt2$comb]

### main plots
p2 <- Heatmap(m2a, col=mycol,
   cluster_rows=F, cluster_columns=F,
   show_column_dend=F, column_dend_height=grid::unit(2, "cm"),   
   show_row_names=T, row_names_gp=gpar(fontsize=7),
   show_column_names=T, column_names_gp=gpar(fontsize=6),
   top_annotation=col_ha,
   heatmap_legend_param=list(title="loading",
      title_gp=gpar(fontsize=9),
      at=seq(-1, 1, by=0.5), 
      labels_gp=gpar(fontsize=9),
      grid_width=grid::unit(0.45, "cm"),
      legend_height=grid::unit(7.8, "cm")))
 
###
figfn <- paste(outdir2, "Figure4.6_loading_TF-regulated-gene.heatmap.png", sep="")
## ggsave(figfn, plot=p2, device="png", width=680, height=700, units="px", dpi=300) ## ggsave 
png(figfn, width=720, height=720,res=120)
set.seed(0)
p2 <- draw(p2)
dev.off()



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
              upper=list(continuous=wrap("cor", size=2)),
              diag=list(continuous=my_box),
              lower=list(continuous=wrap("points", size=0.8)))+
    scale_color_manual(values=col2)+
    theme_bw()+
    theme(strip.text=element_text(size=12))
    

figfn <- paste(outdir2, "Figure4.7_motif_loading.pair.png", sep="")
ggsave(figfn, plot=p0, width=800, height=720, units="px", dpi=120)



###
### plot data for TF-regulated genes

plotDF2 <- m2[,1:5]%>%as.data.frame()%>%rownames_to_column(var="comb")
x <- str_split(plotDF2$comb, "_", simplify=T)
plotDF2 <- plotDF2%>%mutate(MCls=paste(x[,1], x[,2], sep="_"), contrast=x[,3])


###
###
p2 <- ggpairs(plotDF2, columns=2:6, aes(color=contrast),
              upper=list(continuous=wrap("cor", size=2)),
              diag=list(continuous=my_box),
              lower=list(continuous=wrap("points", size=0.8)))+
    scale_color_manual(values=col2)+
    theme_bw()+
    theme(strip.text=element_text(size=12))
    

figfn <- paste(outdir2, "Figure4.7_TF-regulated-gene_loading.pair.png", sep="")
ggsave(figfn, plot=p2, width=780, height=720, units="px", dpi=120)





###################
### cluster CCA ###
###################

 
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
fn <- "./4_motif_plots.outs/2_all_response/4_comb2_cl6_motif_reorder.xlsx"
df_cl <- read.xlsx(fn)
## df_cl2 <- df_cl%>%filter(gene%in%rnSel)



###
fn <- "./4_motif_plots.outs/2_all_response/4_TF_activity_TF-regulated-gene.mat.rds"
x <- read_rds(fn)
rnSel <- rownames(x)

comb2 <- gsub("^[12]_", "", colnames(x))
colnames(x) <- comb2




###
x1 <- x[,1:48]
## x1 <- scale(x1)

###
x2 <- x[,49:96]
## x2 <- scale(x2)


### cca object
Var <- getVar(x1, x2, df_cl)
cca <- getCCA(Var)
## opfn <- paste(outdir2, "4_cca.rds", sep="")
## write_rds(cca, file=opfn)

wa <- cca$wa
xscore <- x1%*%wa
##
wb <- cca$wb
yscore <- x2%*%wb




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
    

figfn <- paste(outdir2, "Figure4.4_motif_score.pair.png", sep="")
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
    

figfn <- paste(outdir2, "Figure4.4_TF-regulated-gene_score.pair.png", sep="")
ggsave(figfn, plot=p2, width=780, height=720, units="px", dpi=120)



########################################
### correlation of score and x and y ###
########################################


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
              upper=list(continuous=wrap("cor", size=2)),
              diag=list(continuous=my_box),
              lower=list(continuous=wrap("points", size=0.8)))+
    scale_color_manual(values=col2)+
    theme_bw()+
    theme(strip.text=element_text(size=12))
    

figfn <- paste(outdir2, "Figure4.7_motif_loading.pair.png", sep="")
ggsave(figfn, plot=p0, width=800, height=720, units="px", dpi=120)



###
### plot data for TF-regulated genes

plotDF2 <- m2[,1:5]%>%as.data.frame()%>%rownames_to_column(var="comb")
x <- str_split(plotDF2$comb, "_", simplify=TRUE)
plotDF2 <- plotDF2%>%mutate(MCls=paste(x[,1], x[,2], sep="_"), contrast=x[,3])


###
###
p2 <- ggpairs(plotDF2, columns=2:6, aes(color=contrast),
              upper=list(continuous=wrap("cor", size=2)),
              diag=list(continuous=my_box),
              lower=list(continuous=wrap("points", size=0.8)))+
    scale_color_manual(values=col2)+
    theme_bw()+
    theme(strip.text=element_text(size=12))
    

figfn <- paste(outdir2, "Figure4.7_TF-regulated-gene_loading.pair.png", sep="")
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
cc_mine <- getCCA(Var_list)

wa <- cc_mine$wa
xscore <- x1%*%wa

##
wb <- cc_mine$wb
yscore <- x2%*%wb


###
###
xscore_0 <- cc2$xscores
plotDF <- data.frame(x=xscore_0[,1], y=xscore[,1])
p0 <- ggplot(plotDF, aes(x=x, y=y))+
   geom_point(size=1.5, shape=1, color="#dd1c77")+
   xlab("xscore of cc fun")+ylab("xscore of my fun")+
   theme_bw()

figfn <- paste(outdir2, "Figure888_xscore.png", sep="")
ggsave(figfn, p0, width=320, height=320, units="px", dpi=120)


### compare
yscore_0 <- cc2$yscore
plotDF2 <- data.frame(x=yscore_0[,1], y=yscore[,1])
p2 <- ggplot(plotDF, aes(x=x, y=y))+
   geom_point(size=1.5, shape=1, color="#dd1c77")+
   xlab("yscore of cc fun")+ylab("yscore of my fun")+
   theme_bw()

figfn <- paste(outdir2, "Figure888_yscore.png", sep="")
ggsave(figfn, p2, width=320, height=320, units="px", dpi=120)

   

               






########################################
### Heatmap of coefficients from cca ###
########################################
## coef_motif <- cc0$xcoef
## coef_rna <- cc0$ycoef

## m1 <- t(coef_motif)
## m2 <- t(coef_rna)



## #### setting colors
## y <- as.vector(m1)
## ##y0 <- y[abs(y)<2] ### quantile(abs(y), probs=0.99)
## mybreak <- c(min(y), seq(-1, 1, length.out=98), max(y)) 
## mycol <- colorRamp2(mybreak, colorRampPalette(rev(brewer.pal(n=7,name="RdBu")))(100))

## quantile(abs(y), probs=c(0.9, 0.95, 0.99))


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
##    show_row_names=F, 
##    show_column_names=T, column_names_gp=gpar(fontsize=6),
##    top_annotation=col_ha,
##    heatmap_legend_param=list(title="coef",
##       title_gp=gpar(fontsize=9),
##       at=seq(-1.5, 1.5, by=0.5), 
##       labels_gp=gpar(fontsize=9),
##       grid_width=grid::unit(0.45, "cm"),
##       legend_height=grid::unit(7.8, "cm")))
 
## ###
## figfn <- paste(outdir2, "Figure4.8_coef_motif.heatmap.png", sep="")
## ## ggsave(figfn, plot=p2, device="png", width=680, height=700, units="px", dpi=300) ## ggsave 
## png(figfn, width=650, height=720,res=120)
## set.seed(0)
## p1 <- draw(p1)
## dev.off()


## ###
## ### m2, TF-regulated genes

## #### setting colors
## y <- as.vector(m2)
## ##y0 <- y[abs(y)<2] ### quantile(abs(y), probs=0.99)
## mybreak <- c(min(y), seq(-4, 4, length.out=98), max(y)) 
## mycol <- colorRamp2(mybreak, colorRampPalette(rev(brewer.pal(n=7,name="RdBu")))(100))

## quantile(abs(y), probs=c(0.9, 0.95, 0.99))


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
##    show_row_names=F, 
##    show_column_names=T, column_names_gp=gpar(fontsize=6),
##    top_annotation=col_ha,
##    heatmap_legend_param=list(title="coef",
##       title_gp=gpar(fontsize=9),
##       at=seq(-6, 6, by=2), 
##       labels_gp=gpar(fontsize=9),
##       grid_width=grid::unit(0.45, "cm"),
##       legend_height=grid::unit(7.8, "cm")))
 
## ###
## figfn <- paste(outdir2, "Figure4.8_coef_TF-regulated-gene.heatmap.png", sep="")
## ## ggsave(figfn, plot=p2, device="png", width=680, height=700, units="px", dpi=300) ## ggsave 
## png(figfn, width=650, height=720,res=120)
## set.seed(0)
## p2 <- draw(p2)
## dev.off()


















######################
#### DEs vs DP-TFs ###
######################


## rm(list=ls())

## outdir2 <- "./4_motif_plots.outs/2_all_response/"

## fn <- paste(outdir2, "5_LFC.DP.TF_LFC.DE.TF.rds", sep="")
## x <- read_rds(fn)



## res <- x%>%dplyr::select(comb, MCls, contrast, beta_x=beta_peak, beta_y=beta_rna)%>%drop_na(beta_x, beta_y)


## col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
##   "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
##    "6_Monocyte"="#984ea3", "7_dnT"="black")
## col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
##        "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3")

## corr <- cor.test(res$beta_y, res$beta_x, method="spearman")
## rr <- round(as.numeric(corr$estimate), digits=3)
## eq <- deparse(bquote(italic(rho)==.(rr)~"***"))
  
     
## p <- ggplot(res, aes(x=beta_x, y=beta_y))+
##    geom_point(aes(color=factor(MCls)), size=0.6)+
##    annotate("text", label=eq, x=0.5, y=5, size=3, parse=T)+ 
##    scale_color_manual(values=col1, guide=guide_legend(override.aes=list(size=1.5)))+
##    geom_smooth(method="lm", formula=y~x, linewidth=0.5, se=F, color="grey50")+ 
##    xlab("Response changes of DARs with TFs")+ ##xlim(-3, 3)+ylim(-3, 9)+ 
##    ylab("Response changes of TF-regulated genes")+ 
##    theme_bw()+
##    theme(legend.key.size=grid::unit(0.4, "cm"),
##          legend.title=element_blank(),
##          legend.text=element_text(size=8))
 
## ###
## figfn <- paste(outdir2, "Figure6.0_scatter_DE_DP.Motif.png", sep="")
## ggsave(figfn, p, width=500, height=380, units="px", dpi=120)






## ########################################
## ### facet by cell type and treatment ###
## ########################################

## ### annotation text

## feq <- function(x){
##   r <- round(as.numeric(x$estimate),digits=3)
##   p <- x$p.value
##   if(p<0.001) symb <- "***"
##   if(p>=0.001 & p<0.01) symb <- "**"
##   if (p>=0.01 & p<0.05) symb <- "*"
##   if(p>0.05) symb <- "NS"
  
##   eq <- bquote(italic(rho)==.(r)~.(symb))
##   eq 
## }


## ###    
## ### annotation text 
## anno_df1 <- res%>%group_by(MCls, contrast)%>%
##    nest()%>%
##    mutate(corr=map(data, ~cor.test((.x)$beta_x, (.x)$beta_y, method="spearman")),
          
##           eq=map(corr,feq),
##           rr=map_dbl(corr,~(.x)$estimate), xpos=-0.5, ypos=0)%>%
##    dplyr::select(-data,-corr)


## ####
## p2 <- ggplot(res, aes(beta_x, beta_y, color=MCls))+
##    geom_point(size=0.6)+
##    scale_color_manual(values=col1)+    
##    geom_text(data=anno_df1, aes(x=xpos, y=ypos, label=eq), color="black", size=2.5, parse=T)+
##    facet_grid(contrast~MCls, scales="free")+
##    xlab("Response changes of DARs with TF")+ ##xlim(-3, 3)+
##    ylab("Response changes of TF-regulated genes")+##ylim(-3, 9)+
##    geom_smooth(method="lm", formula=y~x, linewidth=0.5, se=F, color="grey50")+    
##    theme_bw()+
##    theme(legend.position="none") 
    
 
## ### output 
## figfn <- paste(outdir2, "Figure6.2_facet_scatter_DE_DP.Motif.png", sep="")
## ggsave(figfn, p2, width=1100, height=720, units="px", dpi=120)




########################
### show some motifs ###
########################

## rm(list=ls())

## outdir2 <- "./4_motif_plots.outs/Examples/"
## if ( !file.exists(outdir2) ) dir.create(outdir2, showWarnings=F, recursive=T)


## ## fn <- "./4_motif_plots.outs/3.2_motif_diff_fitInd_JW.rds"
## ## res <- read_rds(fn)%>%mutate(comb=paste(MCls, contrast, sep="_"))


## ###
## ### Function-get the data
## getData <- function(gene, mat, cvt){
##    ##
##    ##gene <- "MA0098.3" 
##    cvt$y <- as.numeric(mat[gene,])
##    MCls <- sort(unique(cvt$MCls))

##    cvt <- cvt%>%group_by(MCls, treats, ind)%>%summarise(y=mean(y, na.rm=T), .groups="drop")%>%ungroup()   
##    ##
##    cvt2 <- map_dfr(MCls, function(oneMCl){
##       ##
##       x <- cvt%>%filter(MCls==oneMCl)
##       x0 <- x%>%filter(treats=="control")
##       y0 <- x0$y
##       names(y0) <- x0$ind
##       ### 
##       x <- x%>%mutate("y0"=y0[as.character(ind)])
##       x
##    })
##    cvt2 <- cvt2%>%mutate(y2=y-y0)
##    cvt2
## }    


## ###
## ### setting colors
## col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
##   "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
##    "6_Monocyte"="#984ea3", "7_dnT"="black")
## col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
##        "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3")



## ###
## ### cell-types and values
## MCls_sel <- c("0_CD4Naive", "5_CD8Naive", "1_TCM", "3_TEM", "2_NKcell", "4_Bcell", "6_Monocyte")
## MCls_val <- 1:length(MCls_sel)
## names(MCls_val) <- MCls_sel



## #######################
## ### motifs activity ###
## #######################


## ###
## ### diff motifs
## fn <- "./4_motif_plots.outs/response_motif_correct/2.2_response_motif_th0.1.txt"
## df_motif <- read.table(fn, header=T)


## ###
## fn <- "./4.0_motif_compare.outs/2_motif.ave_macs2_0.1_cn_control.rds"
## mat <- read_rds(fn)

## x <- str_split(colnames(mat), "_", simplify=T)
## cvt <- data.frame(MCls=paste(x[,1], x[,2], sep="_"), treats=x[,3], ind=x[,4])%>%
##     mutate(treat2=ifelse(treats%in%c("etOH", "water"), "control", treats))




## motif0 <- "ETS1"
## motif_ID <- "MA0098.3"
## df_motif%>%filter(motif_name==motif0)%>%distinct(motif_ID, .keep_all=T)



## ###
## ### generate plot data
## plotDF <- getData(motif_ID, mat, cvt)
## plotDF2 <- plotDF%>%filter(treats%in%c("caffeine", "nicotine"), MCls%in%MCls_sel)
## plotDF2 <- plotDF2%>%
##     mutate("MCls_val"=as.numeric(MCls_val[MCls]), MCl2=fct_reorder(MCls, MCls_val))


## ###
## ###
## comb_sigs <- df_motif%>%filter(motif_name==motif0)%>%pull(comb) 

## sigs_df <- plotDF2%>%group_by(MCls, treats)%>%
##     summarize(ymax=max(y2), .groups="drop")%>%
##     mutate(MCls_value=as.numeric(MCls_val[MCls]),
##            MCl2=fct_reorder(MCls, MCls_value), comb=paste(MCls, treats, sep="_"),
##            sigs=ifelse(comb%in%comb_sigs, "Sig", "Not"))

## p0 <- ggplot(plotDF2, aes(x=factor(MCl2), y=y2, fill=MCls))+
##    ##geom_violin(position="dodge", width=0.8)+
##    geom_boxplot(outlier.size=0.8)+    
##    geom_text(data=sigs_df, aes(x=MCl2, y=ymax+0.05, label=sigs), size=3)+
##    scale_y_continuous(paste("Motif activity changes(", motif0, ")", sep=""),
##                       expand=expansion(mult=c(0.2, 0.2)))+
##    scale_fill_manual("", values=col1)+ 
##    facet_wrap(~treats)+
##    theme_bw()+
##    theme(## axis.text.x=element_text(angle=-90, hjust=0, size=10),
##          axis.text.x=element_text(angle=45, hjust=1, size=8),
##          axis.title.x=element_blank(),
##          axis.title.y=element_text(size=9),
##          axis.text.y=element_text(size=9),
##          strip.text=element_text(size=10), 
##          ##plot.title=element_text(hjust=0.5, size=12),
##          legend.position="none")
 
## figfn <- paste(outdir2, "Figure1.1_", motif0, "_TF_box.png", sep="")
## ggsave(figfn, p0, device="png", width=520, height=380, units="px", dpi=120)


## ###
## ### gene expression
## fn <- "./sc_multiome_data/2_Differential/1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn_sampleID+new_treat.rds"
## res <- read_rds(fn)%>%as.data.frame()%>%
##     mutate(comb=paste(MCls, contrast, sep="_"), is_sig=ifelse(p.adjusted<0.1&abs(estimate)>0.5, 1, 0))   

## fn <- "./sc_multiome_data/2_Differential/1_DiffRNA_2.outs/1_YtX.sel_clusters_0.1_cn.rds"
## mat <- read_rds(fn)
## rnz <- colSums(mat)
## mat2 <- sweep(mat, 2, rnz, "/")
## mat2 <- log2(mat2*1e+06+1)

## x <- str_split(colnames(mat2), "_", simplify=T)
## cvt <- data.frame(MCls=paste(x[,1], x[,2], sep="_"), treats=x[,3], ind=x[,4])
## cvt <- cvt%>%mutate(treats=ifelse(treats%in%c("etOH", "water"), "control", treats))




## gene0 <- "ETS1"
## plotDF <- getData(gene0, mat2, cvt)
## plotDF2 <- plotDF%>%filter(treats%in%c("caffeine", "nicotine"), MCls%in%MCls_sel)
## plotDF2 <- plotDF2%>%
##     mutate("MCls_val"=as.numeric(MCls_val[MCls]), MCl2=fct_reorder(MCls, MCls_val))%>%drop_na(y2)


## ###
## ###
## comb_sigs <- res%>%filter(is_sig==1, gene==gene0)%>%pull(comb) 

## sigs_df <- plotDF2%>%group_by(MCls, treats)%>%
##     summarize(ymax=max(y2), .groups="drop")%>%
##     mutate(MCls_value=as.numeric(MCls_val[MCls]),
##            MCl2=fct_reorder(MCls, MCls_value), comb=paste(MCls, treats, sep="_"),
##            sigs=ifelse(comb%in%comb_sigs, "Sig", "Not"))
 
## p2 <- ggplot(plotDF2, aes(x=factor(MCl2), y=y2, fill=MCls))+
##    ##geom_violin(position="dodge", width=0.8)+
##    geom_boxplot(outlier.size=0.8)+    
##    geom_text(data=sigs_df, aes(x=MCl2, y=ymax+0.05, label=sigs), size=3)+
##    scale_y_continuous(bquote("Gene expression changes("~italic(.(gene0))~")"),
##                       expand=expansion(mult=c(0.2, 0.2)))+
##    scale_fill_manual("", values=col1)+ 
##    facet_wrap(~treats)+
##    theme_bw()+
##    theme(## axis.text.x=element_text(angle=-90, hjust=0, size=10),
##          axis.text.x=element_text(angle=45, hjust=1, size=8),
##          axis.title.x=element_blank(),
##          axis.title.y=element_text(size=9),
##          axis.text.y=element_text(size=9),
##          strip.text=element_text(size=10), 
##          ##plot.title=element_text(hjust=0.5, size=12),
##          legend.position="none")
  
## figfn <- paste(outdir2, "Figure1.2_", gene0, "_gene_box.png", sep="")
## ggsave(figfn, p2, device="png", width=520, height=380, units="px", dpi=120)




###
### MTF1 genes  
###
### diff motifs
## fn <- "./4_motif_plots.outs/response_motif_correct/2.2_response_motif_th0.1.txt"
## df_motif <- read.table(fn, header=T)


## ###
## fn <- "./4.0_motif_compare.outs/2_motif.ave_macs2_0.1_cn_control.rds"
## mat <- read_rds(fn)

## x <- str_split(colnames(mat), "_", simplify=T)
## cvt <- data.frame(MCls=paste(x[,1], x[,2], sep="_"), treats=x[,3], ind=x[,4])%>%
##     mutate(treat2=ifelse(treats%in%c("etOH", "water"), "control", treats))


 

## motif0 <- "MTF1"
## motif_ID <- "MA0863.1"
## df_motif%>%filter(motif_name==motif0)

###
### generate plot data
## plotDF <- getData(motif_ID, mat, cvt)
## plotDF2 <- plotDF%>%filter(treats%in%c("zinc"), MCls%in%MCls_sel)
## plotDF2 <- plotDF2%>%
##     mutate("MCls_val"=as.numeric(MCls_val[MCls]), MCl2=fct_reorder(MCls, MCls_val))

## ###
## ###
## comb_sigs <- df_motif%>%filter(motif_name==motif0)%>%pull(comb) 

## sigs_df <- plotDF2%>%group_by(MCls, treats)%>%
##     summarize(ymax=max(y2), .groups="drop")%>%
##     mutate(MCls_value=as.numeric(MCls_val[MCls]),
##            MCl2=fct_reorder(MCls, MCls_value), comb=paste(MCls, treats, sep="_"),
##            sigs=ifelse(comb%in%comb_sigs, "Sig", "Not"))

## p0 <- ggplot(plotDF2, aes(x=factor(MCl2), y=y2, fill=MCls))+
##    ##geom_violin(position="dodge", width=0.8)+
##    geom_boxplot(outlier.size=0.8)+    
##    geom_text(data=sigs_df, aes(x=MCl2, y=ymax+0.05, label=sigs), size=3)+
##    scale_y_continuous(paste("Motif activity changes(", motif0, ")", sep=""),
##                       expand=expansion(mult=c(0.2, 0.2)))+
##    scale_fill_manual("", values=col1)+ 
##    ##facet_wrap(~treats)+
##    theme_bw()+
##    theme(## axis.text.x=element_text(angle=-90, hjust=0, size=10),
##          axis.text.x=element_text(angle=45, hjust=1, size=8),
##          axis.title.x=element_blank(),
##          axis.title.y=element_text(size=9),
##          axis.text.y=element_text(size=9),
##          ##strip.text=element_text(size=10), 
##          ##plot.title=element_text(hjust=0.5, size=12),
##          legend.position="none")
 
## figfn <- paste(outdir2, "Figure2.1_", motif0, "_TF_box.png", sep="")
## ggsave(figfn, p0, device="png", width=380, height=380, units="px", dpi=120)



###
### genes expression 


###
### gene expression
## fn <- "./sc_multiome_data/2_Differential/1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn_sampleID+new_treat.rds"
## res <- read_rds(fn)%>%as.data.frame()%>%
##     mutate(comb=paste(MCls, contrast, sep="_"), is_sig=ifelse(p.adjusted<0.1&abs(estimate)>0.5, 1, 0))   

## fn <- "./sc_multiome_data/2_Differential/1_DiffRNA_2.outs/1_YtX.sel_clusters_0.1_cn.rds"
## mat <- read_rds(fn)
## rnz <- colSums(mat)
## mat2 <- sweep(mat, 2, rnz, "/")
## mat2 <- log2(mat2*1e+06+1)

## x <- str_split(colnames(mat2), "_", simplify=T)
## cvt <- data.frame(MCls=paste(x[,1], x[,2], sep="_"), treats=x[,3], ind=x[,4])
## cvt <- cvt%>%mutate(treats=ifelse(treats%in%c("etOH", "water"), "control", treats))



## gene0 <- "MTF1"
## plotDF <- getData(gene0, mat2, cvt)
## plotDF2 <- plotDF%>%filter(treats%in%c("zinc"), MCls%in%MCls_sel)
## plotDF2 <- plotDF2%>%
##     mutate("MCls_val"=as.numeric(MCls_val[MCls]), MCl2=fct_reorder(MCls, MCls_val))%>%drop_na(y2)


###
###
## comb_sigs <-

##     res%>%filter(gene==gene0, contrast=="zinc") ##%>%pull(comb) 

## sigs_df <- plotDF2%>%group_by(MCls, treats)%>%
##     summarize(ymax=max(y2), .groups="drop")%>%
##     mutate(MCls_value=as.numeric(MCls_val[MCls]),
##            MCl2=fct_reorder(MCls, MCls_value), comb=paste(MCls, treats, sep="_"),
##            sigs=ifelse(comb%in%comb_sigs, "Sig", "Not"))
 
## p2 <- ggplot(plotDF2, aes(x=factor(MCl2), y=y2, fill=MCls))+
##    ##geom_violin(position="dodge", width=0.8)+
##    geom_boxplot(outlier.size=0.8)+    
##    ## geom_text(data=sigs_df, aes(x=MCl2, y=ymax+0.05, label=sigs), size=3)+
##    scale_y_continuous(bquote("Gene expression changes("~italic(.(gene0))~")"),
##                       expand=expansion(mult=c(0.2, 0.2)))+
##    scale_fill_manual("", values=col1)+ 
##    facet_wrap(~treats)+
##    theme_bw()+
##    theme(## axis.text.x=element_text(angle=-90, hjust=0, size=10),
##          axis.text.x=element_text(angle=45, hjust=1, size=8),
##          axis.title.x=element_blank(),
##          axis.title.y=element_text(size=9),
##          axis.text.y=element_text(size=9),
##          strip.text=element_text(size=10), 
##          ##plot.title=element_text(hjust=0.5, size=12),
##          legend.position="none")
  
## figfn <- paste(outdir2, "Figure2.2_", gene0, "_gene_box.png", sep="")
## ggsave(figfn, p2, device="png", width=380, height=380, units="px", dpi=120)



## fn <- "./sc_multiome_data/1_processing/3_Clustering.outs/1_seurat.cluster.combined.mitofilt.rds"
## sc0 <- read_rds(fn)


## fn <- "./sc_multiome_data/3_motif/2_motif.activities.outs/1_scATAC.motifActivities.rds"
## sc2 <- read_rds(fn)
 
## sc2[["RNA"]] <- sc0[["RNA"]]
## sc2[["SCT"]] <- sc0[["SCT"]]
## sc2[["pca"]] <- sc0[["pca"]]
## sc2[["lsi"]] <- sc0[["lsi"]]
## sc2[["umap"]] <- sc0[["umap"]]

## identical(Cells(sc0), Cells(sc2))

## ### add column in meta data
## x <- sc2@meta.data
## x2 <- x%>%mutate(treat2=ifelse(treats%in%c("etOH", "water"), "control", treats))
## sc2 <- AddMetaData(sc2, metadata=x2)

## opfn <- "./sc_multiome_data_rs/1_seurat_combined_all_data.rds"
## write_rds(sc2, file=opfn)


## ###
## ### plots

## outdir2 <- "./4_motif_plots.outs/Examples/"
## if ( !file.exists(outdir2) ) dir.create(outdir2, showWarnings=F, recursive=T)



## ### examples motifs
## fn <- "./4_motif_plots.outs/response_motif_correct/2.2_response_motif_th0.1.txt"
## df_motif <- read.table(fn, header=T)

## x <- df_motif%>%filter(motif_name=="ETS1")%>%distinct(motif_ID, .keep_all=T)

## sc <- sc2
## DefaultAssay(sc) <- "chromvar"

## motif0 <- "MA0098.3"
## gene <- "ETS1"


## yval <- FetchData(sc, vars=motif0)

## sc_tr <- subset(sc, treat2%in%c("caffeine", "nicotine"))
## sc_ctr <- subset(sc, treat2=="control")
## DefaultAssay(sc_tr) <- "chromvar"
## DefaultAssay(sc_ctr) <- "chromvar"

## ### TF
## p0 <- FeaturePlot(sc_ctr, feature=motif0)+
##     scale_colour_gradient(gene, low="lightgrey", high="blue", limits=c(-5.5, 5.5))+
##     ggtitle("Caffeine and nicotine")+
##     theme_bw()+
##     theme(legend.title=element_text(size=8),
##           legend.key.size=grid::unit(0.3, "cm"),
##           legend.text=element_text(size=8),
##           plot.title=element_text(hjust=0.5, size=10),
##           axis.title=element_text(size=9),
##           axis.text=element_text(size=9))
    
## p2 <- FeaturePlot(sc_tr, features=motif0)+
##     scale_colour_gradient(gene, low="lightgrey", high="blue", limits=c(-5.5, 5.5))+
##     ggtitle("Control")+
##     theme_bw()+
##     theme(legend.title=element_text(size=8),
##           legend.key.size=grid::unit(0.3, "cm"),
##           legend.text=element_text(size=8),
##           plot.title=element_text(hjust=0.5, size=10),
##           axis.title=element_text(size=9),
##           axis.text=element_text(size=9))

## comb <- plot_grid(p0, p2)
 
## figfn <- paste(outdir2, "Figure4.1_ETS1_TFs_umap.png", sep="")
## ggsave(figfn, comb, device="png", units="px", height=320, width=700, dpi=120)



## ## gene expression
## DefaultAssay(sc) <- "RNA"
## yval <- FetchData(sc, vars=gene)
## DefaultAssay(sc_tr) <- "RNA"
## DefaultAssay(sc_ctr) <- "RNA"
## p0 <- FeaturePlot(sc_ctr, feature=gene)+
##     scale_colour_gradient(gene, low="lightgrey", high="blue", limits=c(0,21))+
##     ggtitle("Caffeine and nicotine")+
##     theme_bw()+
##     theme(legend.title=element_text(size=8),
##           legend.key.size=grid::unit(0.3, "cm"),
##           legend.text=element_text(size=8),
##           plot.title=element_text(hjust=0.5, size=12),
##           axis.title=element_text(size=9),
##           axis.text=element_text(size=9))
    
## p2 <- FeaturePlot(sc_tr, features=gene)+
##     scale_colour_gradient(gene, low="lightgrey", high="blue", limits=c(0,21))+
##     ggtitle("Control")+
##     theme_bw()+
##     theme(legend.title=element_text(size=8),
##           legend.key.size=grid::unit(0.3, "cm"),
##           legend.text=element_text(size=8),
##           plot.title=element_text(hjust=0.5, size=10),
##           axis.title=element_text(size=9),
##           axis.text=element_text(size=9))

## comb <- plot_grid(p0, p2)

## figfn <- paste(outdir2, "Figure4.1_ETS1_gene_umap.png", sep="")
## ggsave(figfn, comb, device="png", units="px", height=320, width=700, dpi=120)




###
### summary response motif
## fn <- "./4_motif_plots.outs/3.2_motif_diff_fitInd_JW.rds"
## res <- read_rds(fn)%>%mutate(comb=paste(MCls, contrast, sep="_"), zscore=beta/stderr)
## comb <- sort(unique(res$comb))

## summ <- map_dfr(comb, function(ii){
##    ###
##    cat(ii, "\n")    
##    x <- res%>%filter(comb==ii)
##    th0 <- quantile(abs(x$beta), probs=0.9)
##    motif <- x%>%filter(qval<0.1, abs(beta)>th0)%>%pull(gene)%>%unique()
##    ###
##    summ2 <- data.frame(comb=ii, nmotif=length(motif), th=round(th0, 3))
##    summ2
## })

## opfn <- paste(outdir, "0.1_summary_response_motif.xlsx", sep="")
## write.xlsx(summ, file=opfn, overwrite=T)


## ## union of response motifs
## resp_motif <- lapply(comb, function(ii){
##    ###
##    cat(ii, "\n")    
##    x <- res%>%filter(comb==ii)
##    th0 <- quantile(abs(x$beta), probs=0.9)
##    motif <- x%>%filter(qval<0.1, abs(beta)>th0)%>%pull(gene)%>%unique()
##    motif
## })
## resp_motif <- unique(unlist(resp_motif))


## ###
## fn <- "./4_motif_plots.outs/response_motif_correct/2.2_response_motif_th0.1.txt"
## df_motif <- read.table(fn, header=T)

## motif_zinc <- df_motif$motif_ID



