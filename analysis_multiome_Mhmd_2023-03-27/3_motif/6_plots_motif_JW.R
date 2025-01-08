###
library(tidyverse)
##
library(cowplot)
library(grid)
library(gridExtra)
library(ggExtra)
library(RColorBrewer)
library(ggsci)
library(viridis)
library(ComplexHeatmap)
library(circlize)
library(openxlsx)



## motifLs <- unique(read.table("./Response_motif/2.2_response_motif_th0.1.txt", header=T)$motif_ID)
## motifLs <- unique(read.table("./Response_motif/2.3_response_motif_th0.2.txt", header=T)$motif_ID)

###
###
fn <- "./Response_motif/1_diff_motif.results.rds"
res <- read_rds(fn)
res2 <- res%>%filter(!grepl("8_MAIT|water", comb))%>%arrange(comb)

###
df <- res2%>%pivot_wider(id_cols=gene, names_from=comb, values_from=beta)
df2 <- df[,-1]
mat <- cor(df2)

###
comb <- sort(unique(res2$comb))
ncomb <- length(comb)

### significance 
is_sig <- matrix(0, ncomb, ncomb)
for (i in 1:ncomb){
    ##
    for (j in 1:ncomb){
        x <- df2[, i, drop=T]
        y <- df2[, j, drop=T]
        pval <- cor.test(x,y)$p.value
        is_sig[i,j] <- ifelse(pval<0.05, 1, NA)
   }
}

mat2 <- mat*is_sig


###
### treatment together 
treat2 <- rep(c("caffeine", "nicotine", "vitA", "vitD", "vitE", "zinc"), each=8)
MCl2 <- rep(c("0_CD4Naive", "1_TCM", "3_TEM", "5_CD8Naive", "7_dnT", "2_NKcell", "4_Bcell", "6_Monocyte"), times=6)
comb2 <- paste(MCl2, treat2, sep="_")

mat2 <- mat2[comb2, comb2]

###
### set colors
mycol <- colorRamp2(seq(-1, 1, length.out=20), colorRampPalette(rev(brewer.pal(n=7,name="RdBu")))(20))

### annotation
cvt <- str_split(rownames(mat2), "_", simplify=T)
anno_df <- data.frame(celltype=paste(cvt[,1], cvt[,2], sep="_"), treats=cvt[,3])

column_ha <- HeatmapAnnotation(df=anno_df,
   col=list(celltype=c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56", "3_TEM"="blue",
                   "4_Bcell"="#4daf4a", "5_CD8Naive"="green", "6_Monocyte"="#984ea3", "7_dnT"="black"),
            treats=c("caffeine"="red", "nicotine"="tan", "vitA"="tan4", "vitD"="seagreen4",
                         "vitE"="salmon3", "zinc"="maroon3")),
   annotation_legend_param=list(
      celltype=list(labels_gp=gpar(fontsize=10), title_gp=gpar(fontsize=10),
                grid_width=grid::unit(0.5, "cm"), grid_height=grid::unit(0.6, "cm") ),
      treats=list(labels_gp=gpar(fontsize=10), title_gp=gpar(fontsize=10),
                  grid_width=grid::unit(0.5, "cm"), grid_height=grid::unit(0.6, "cm"))),
   annotation_name_gp=gpar(fontsize=12))




### Heatmap
p <- Heatmap(mat2, name="PCC", na_col="grey90", 
    col=mycol, cluster_rows=F, cluster_columns=F,
    show_row_dend=F, show_column_dend=F,
    top_annotation=column_ha,
    show_row_names=T, show_column_names=T,
    row_names_gp=gpar(fontsize=9),
    column_names_gp=gpar(fontsize=9), column_names_rot=-90,
    heatmap_legend_param=list(at=seq(-1, 1, by=0.25),
        grid_width=grid::unit(0.5, "cm"), legend_height=grid::unit(10, "cm"),
        title_gp=gpar(fontsize=10), labels_gp=gpar(fontsize=9)),
    ## cell_fun=function(j, i, x, y, width, height, fill){
    ##    grid.text(round(mat[i,j],digits=3), x, y, gp=gpar(fontsize=10))
    ## },
    use_raster=F, raster_device="png")


###
###
## figfn <- "./Response_motif/Figure1.2_th0.1_corr_heat.png"
## figfn <- "./Response_motif/Figure1.3_th0.2_corr_heat.png"
figfn <- "./Response_motif/Figure1_corr_heat.png"
png(figfn, height=900, width=1050, res=120)
set.seed(0)
p2 <- draw(p)
dev.off()




##################################
### cell type cluster based on ###
##################################

treat2 <- rep(c("caffeine", "nicotine", "vitA", "vitD", "vitE", "zinc"), times=8)
MCl2 <- rep(c("0_CD4Naive", "1_TCM", "3_TEM", "5_CD8Naive", "7_dnT", "2_NKcell", "4_Bcell", "6_Monocyte"), each=6)
comb2 <- paste(MCl2, treat2, sep="_")

mat2 <- mat2[comb2, comb2]

###
### set colors
mycol <- colorRamp2(seq(-1, 1, length.out=20), colorRampPalette(rev(brewer.pal(n=7,name="RdBu")))(20))

### annotation
cvt <- str_split(rownames(mat2), "_", simplify=T)
anno_df <- data.frame(celltype=paste(cvt[,1], cvt[,2], sep="_"), treats=cvt[,3])

column_ha <- HeatmapAnnotation(df=anno_df,
   col=list(celltype=c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56", "3_TEM"="blue",
                   "4_Bcell"="#4daf4a", "5_CD8Naive"="green", "6_Monocyte"="#984ea3", "7_dnT"="black"),
            treats=c("caffeine"="red", "nicotine"="tan", "vitA"="tan4", "vitD"="seagreen4",
                         "vitE"="salmon3", "zinc"="maroon3")),
   annotation_legend_param=list(
      celltype=list(labels_gp=gpar(fontsize=10), title_gp=gpar(fontsize=10),
                grid_width=grid::unit(0.5, "cm"), grid_height=grid::unit(0.6, "cm") ),
      treats=list(labels_gp=gpar(fontsize=10), title_gp=gpar(fontsize=10),
                  grid_width=grid::unit(0.5, "cm"), grid_height=grid::unit(0.6, "cm"))),
   annotation_name_gp=gpar(fontsize=12))




### Heatmap
p <- Heatmap(mat2, name="PCC", na_col="grey90", 
    col=mycol, cluster_rows=F, cluster_columns=F,
    show_row_dend=F, show_column_dend=F,
    top_annotation=column_ha,
    show_row_names=T, show_column_names=T,
    row_names_gp=gpar(fontsize=9),    
    column_names_gp=gpar(fontsize=9), column_names_rot=-90,
    heatmap_legend_param=list(at=seq(-1, 1, by=0.25),
        grid_width=grid::unit(0.5, "cm"), legend_height=grid::unit(10, "cm"),
        title_gp=gpar(fontsize=10), labels_gp=gpar(fontsize=9)),
    ## cell_fun=function(j, i, x, y, width, height, fill){
    ##    grid.text(round(mat[i,j],digits=3), x, y, gp=gpar(fontsize=10))
    ## },
    use_raster=F, raster_device="png")


###
###
## figfn <- "./Response_motif/Figure1.2_th0.1_corr_heat.png"
figfn <- "./Response_motif/Figure2_corr_heat.png"
png(figfn, height=900, width=1050, res=120)
set.seed(0)
p2 <- draw(p)
dev.off()



###################################################
### overlap of response motif across cell types ###
###################################################

fn <- "./Response_motif/2.2_response_motif_th0.1.txt"
## fn <- "./Response_motif/2.3_response_motif_th0.2.txt"
res <- read.table(fn, header=T)
res2 <- res%>%filter(!grepl("8_MAIT|water", comb))%>%arrange(comb)
###
###

comb_ls <- sort(unique(res2$comb))
ncomb <- length(comb_ls) 
mat_olap <- matrix(0, ncomb, ncomb)
for (i in 1:ncomb){
   ## 
   motif1 <- res2%>%filter(comb==comb_ls[i])%>%pull(motif_ID)
   ##
   for (j in i:ncomb){
      ###
      motif2 <- res2%>%dplyr::filter(comb==comb_ls[j])%>%pull(motif_ID)
      nolap <- length(intersect(motif1, motif2)) 
      mat_olap[i,j] <- nolap
      mat_olap[j,i] <- nolap
   }    
}
colnames(mat_olap) <- comb_ls
rownames(mat_olap) <- comb_ls

treat2 <- rep(c("caffeine", "nicotine", "vitA", "vitD", "vitE", "zinc"), each=8)
MCl2 <- rep(c("0_CD4Naive", "1_TCM", "3_TEM", "5_CD8Naive", "7_dnT", "2_NKcell", "4_Bcell", "6_Monocyte"), times=6)
comb2 <- paste(MCl2, treat2, sep="_")
comb2 <- comb2[comb2%in%comb_ls]

mat2 <- mat_olap[comb2, comb2]

###
### set colors
olap_min <- min(mat2)
olap_max <- max(mat2)
mycol <- colorRamp2(seq(olap_min, olap_max, length.out=20), colorRampPalette(brewer.pal(n=7,name="BuGn"))(20))

### annotation
cvt <- str_split(rownames(mat2), "_", simplify=T)
anno_df <- data.frame(celltype=paste(cvt[,1], cvt[,2], sep="_"), treats=cvt[,3])

column_ha <- HeatmapAnnotation(df=anno_df,
   col=list(celltype=c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56", "3_TEM"="blue",
                   "4_Bcell"="#4daf4a", "5_CD8Naive"="green", "6_Monocyte"="#984ea3", "7_dnT"="black"),
            treats=c("caffeine"="red", "nicotine"="tan", "vitA"="tan4", "vitD"="seagreen4",
                         "vitE"="salmon3", "zinc"="maroon3")),
   annotation_legend_param=list(
      celltype=list(labels_gp=gpar(fontsize=10), title_gp=gpar(fontsize=10),
                grid_width=grid::unit(0.5, "cm"), grid_height=grid::unit(0.6, "cm") ),
      treats=list(labels_gp=gpar(fontsize=10), title_gp=gpar(fontsize=10),
                  grid_width=grid::unit(0.5, "cm"), grid_height=grid::unit(0.6, "cm"))),
   annotation_name_gp=gpar(fontsize=12))



 
### Heatmap
p <- Heatmap(mat2, name="olap", na_col="grey90", 
    col=mycol, cluster_rows=F, cluster_columns=F,
    show_row_dend=F, show_column_dend=F,
    top_annotation=column_ha,
    show_row_names=T, show_column_names=T,
     row_names_gp=gpar(fontsize=9),
    column_names_gp=gpar(fontsize=9), column_names_rot=-90,
    heatmap_legend_param=list(at=seq(0, olap_max, by=20),
        grid_width=grid::unit(0.5, "cm"), legend_height=grid::unit(10, "cm"),
        title_gp=gpar(fontsize=10), labels_gp=gpar(fontsize=9)),
    ## cell_fun=function(j, i, x, y, width, height, fill){
    ##    grid.text(round(mat[i,j],digits=3), x, y, gp=gpar(fontsize=10))
    ## },
    use_raster=F, raster_device="png")


###
###
figfn <- "./Response_motif/Figure3_olap_heat.png"
## figfn <- "./Response_motif/Figure3.3_th0.2_olap_heat.png"
png(figfn, height=900, width=1050, res=120)
set.seed(0)
p2 <- draw(p)
dev.off()



#####################
### Jaccard index ###
#####################

fn <- "./Response_motif/2.2_response_motif_th0.1.txt"
## fn <- "./Response_motif/2.3_response_motif_th0.2.txt"
res <- read.table(fn, header=T)
res2 <- res%>%filter(!grepl("8_MAIT|water", comb))%>%arrange(comb)
###
###

comb_ls <- sort(unique(res2$comb))
ncomb <- length(comb_ls) 
mat_olap <- matrix(0, ncomb, ncomb)
for (i in 1:ncomb){
   ## 
   motif1 <- res2%>%filter(comb==comb_ls[i])%>%pull(motif_ID)
   ##
   for (j in i:ncomb){
      ###
      motif2 <- res2%>%dplyr::filter(comb==comb_ls[j])%>%pull(motif_ID)
      nolap <- length(intersect(motif1, motif2))
      n_union <- length(union(motif1, motif2))
       
      mat_olap[i,j] <- nolap/n_union
      mat_olap[j,i] <- nolap/n_union
   }    
}
colnames(mat_olap) <- comb_ls
rownames(mat_olap) <- comb_ls

treat2 <- rep(c("caffeine", "nicotine", "vitA", "vitD", "vitE", "zinc"), each=8)
MCl2 <- rep(c("0_CD4Naive", "1_TCM", "3_TEM", "5_CD8Naive", "7_dnT", "2_NKcell", "4_Bcell", "6_Monocyte"), times=6)
comb2 <- paste(MCl2, treat2, sep="_")
comb2 <- comb2[comb2%in%comb_ls]

mat2 <- mat_olap[comb2, comb2]

###
### set colors
olap_min <- min(mat2)
olap_max <- max(mat2)
mycol <- colorRamp2(seq(0, 1, length.out=20), colorRampPalette(brewer.pal(n=7,name="BuGn"))(20))

### annotation
cvt <- str_split(rownames(mat2), "_", simplify=T)
anno_df <- data.frame(celltype=paste(cvt[,1], cvt[,2], sep="_"), treats=cvt[,3])

column_ha <- HeatmapAnnotation(df=anno_df,
   col=list(celltype=c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56", "3_TEM"="blue",
                   "4_Bcell"="#4daf4a", "5_CD8Naive"="green", "6_Monocyte"="#984ea3", "7_dnT"="black"),
            treats=c("caffeine"="red", "nicotine"="tan", "vitA"="tan4", "vitD"="seagreen4",
                         "vitE"="salmon3", "zinc"="maroon3")),
   annotation_legend_param=list(
      celltype=list(labels_gp=gpar(fontsize=10), title_gp=gpar(fontsize=10),
                grid_width=grid::unit(0.5, "cm"), grid_height=grid::unit(0.6, "cm") ),
      treats=list(labels_gp=gpar(fontsize=10), title_gp=gpar(fontsize=10),
                  grid_width=grid::unit(0.5, "cm"), grid_height=grid::unit(0.6, "cm"))),
   annotation_name_gp=gpar(fontsize=12))



 
### Heatmap
p <- Heatmap(mat2, name="Jaccard", na_col="grey90", 
    col=mycol, cluster_rows=F, cluster_columns=F,
    show_row_dend=F, show_column_dend=F,
    top_annotation=column_ha,
    show_row_names=T, show_column_names=T,
  row_names_gp=gpar(fontsize=9),   
    column_names_gp=gpar(fontsize=9), column_names_rot=-90,
    heatmap_legend_param=list(at=seq(0, 1, by=0.2),
        grid_width=grid::unit(0.5, "cm"), legend_height=grid::unit(10, "cm"),
        title_gp=gpar(fontsize=10), labels_gp=gpar(fontsize=9)),
    ## cell_fun=function(j, i, x, y, width, height, fill){
    ##    grid.text(round(mat[i,j],digits=3), x, y, gp=gpar(fontsize=10))
    ## },
    use_raster=F, raster_device="png")


###
###
figfn <- "./Response_motif/Figure4_jaccard_heat.png"
## figfn <- "./Response_motif/Figure4.3_th0.2_jaccard_heat.png"
png(figfn, height=900, width=1050, res=120)
set.seed(0)
p2 <- draw(p)
dev.off()



