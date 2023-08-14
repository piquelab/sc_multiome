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

rm(list=ls())

outdir <- "./4_motif_plots.outs/"
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)



###
### We generated respone motif analysis using Mhmd's script, ./4.0_motif_compare.outs/3.2_motif_diff_fitInd_JW.rds
### Then cp this file to the folder 4_motif_plots.outs for further analysis
#### The folder, ./4_motif_plots.outs/response_motif_correct contain response motifs using the above file
### Aug-8-2023, By JW



###
### motif name
fn <- "./sc_multiome_data/3_motif/2_motif.activities.outs/1_scATAC.motifActivities.rds"
atac <- read_rds(fn)
x <- Motifs(atac)
motifs <- unlist(x@motif.names)
## motifs_name <- unlist(unname(motifs))
 

###
### response motifs
fn <- "./sc_multiome_data/3_motif/2_motif.activities.outs/3_motif.diff.results_macs2_0.1_cn_indfilt_incsampleID_all_cols_control.rds"
res <- read_rds(fn)%>%mutate(comb=paste(MCls, contrast, sep="_"), zscore=beta/stderr)
comb <- sort(unique(res$comb))


### motif matrix
dfmat <- res%>%pivot_wider(id_cols=gene, names_from=comb, values_from=zscore)
mat <- dfmat%>%column_to_rownames(var="gene")%>%as.matrix()


###
### get colnames and re-order by treats
rn <- colnames(mat)
x <- str_split(rn, "_", simplify=T)
cvt <- data.frame(comb=rn, MCls=paste(x[,1], x[,2], sep="_"), contrast=x[,3])
cvt <- cvt%>%arrange(contrast)
mat <- mat[,cvt$comb]


####
#### select motifs

### response motif
top_motif <- lapply(comb, function(ii){
    ##
    x <- res%>%filter(comb==ii)
    th0 <- quantile(abs(x$beta), probs=0.9)
    x2 <- x%>%filter(qval<0.1, abs(beta)>th0)
    ## if ( nrow(x2)>0){
    ##    x2%>%slice_max 
    motif <- x2 %>%pull(gene)%>%unique()
    motif
})
top_motif <- unique(unlist(top_motif))


###
### Final select motif
vm <- apply(mat, 1, var)
vm_top <- sort(vm[top_motif], T)
sel_motif1 <- names(vm_top)[1:80]

## mu_m <- apply(abs(mat), 1, mean)
## mu_top <- sort(mu_m[top_motif],T)
## sel_motif2 <- names(mu_top)[1:80]

## sel_motif <- union(sel_motif1, sel_motif2)

sel_motif <- sel_motif1

###
### plot data
mat2 <- mat[sel_motif,]
rownames(mat2) <- unname(motifs[sel_motif])


### color for heatmap value
y <- as.vector(mat2)
##y0 <- y[abs(y)<2] ### quantile(abs(y), probs=0.99)
mybreak <- c(min(y,na.rm=T), seq(-6, 6, length.out=98), max(y,na.rm=T))

quantile(abs(y), probs=c(0.9,0.95,0.99))
range(y)

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
                grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm"))),
  simple_anno_size=unit(0.5, "cm"))


###
### main heatmap

p1 <- Heatmap(mat2, col=mycol,
   cluster_rows=T, cluster_columns=F,
   show_row_names=T, row_names_gp=gpar(fontsize=6),
   show_column_names=T, column_names_gp=gpar(fontsize=6),
   show_row_dend=F, show_column_dend=F,
   ##column_names_rot=-45,   
   top_annotation=col_ha,
   heatmap_legend_param=list(title=bquote(~italic(Z)~"-score"),
      title_gp=gpar(fontsize=9),
      at=seq(-9, 9, by=3), 
      labels_gp=gpar(fontsize=9),
      grid_width=grid::unit(0.45, "cm"),
      legend_height=grid::unit(7.8, "cm")))
 
###
figfn <- paste(outdir, "Figure1.1_motif.heatmap.png", sep="")
## ggsave(figfn, plot=p2, device="png", width=680, height=700, units="px", dpi=300) ## ggsave 
png(figfn, width=900, height=1000,res=120)
set.seed(0)
p1 <- draw(p1)
dev.off()





### cluster gene and re-order gene
hmap <- Heatmap(mat2, cluster_rows=T, cluster_columns=F)
set.seed(0)
hmap <- draw(hmap)
cl <- row_order(hmap)

geneSel <- rownames(mat2)
DF_cl <- NULL
for (i in 1:length(cl)){
   cl_tmp <- data.frame(cluster=i, gene=geneSel[cl[[i]]])
   DF_cl <- rbind(DF_cl, cl_tmp)
}
opfn <- paste(outdir, "1_gene_reorder.xlsx", sep="")
write.xlsx(DF_cl, file=opfn, overwrite=T)



###
### summary response motif
fn <- "./sc_multiome_data/3_motif/2_motif.activities.outs/3_motif.diff.results_macs2_0.1_cn_indfilt_incsampleID_all_cols_control.rds"
res <- read_rds(fn)%>%mutate(comb=paste(MCls, contrast, sep="_"), zscore=beta/stderr)
comb <- sort(unique(res$comb))

summ <- map_dfr(comb, function(ii){
   ###
   cat(ii, "\n")    
   x <- res%>%filter(comb==ii)
   th0 <- quantile(abs(x$beta), probs=0.9)
   motif <- x%>%filter(qval<0.1, abs(beta)>th0)%>%pull(gene)%>%unique()
   ###
   summ2 <- data.frame(comb=ii, nmotif=length(motif), th=round(th0, 3))
   summ2
})

opfn <- paste(outdir, "0.1_summary_response_motif.xlsx", sep="")
write.xlsx(summ, file=opfn, overwrite=T)




## fn <- "./sc_multiome_data/3_motif/2_motif.activities.outs/3_motif.diff.results_macs2_0.1_cn.rds"
## res2 <- read_rds(fn)%>%mutate(comb=paste(MCls, contrast, sep="_"))

## ## compare
## df_olap <- map_dfr(comb, function(ii){
##    ###
##    cat(ii, "\n")
    
##    x <- res%>%filter(comb==ii)
##    th0 <- quantile(abs(x$beta), probs=0.9)
##    motif <- x%>%filter(qval<0.1, abs(beta)>th0)%>%pull(gene)%>%unique()

##    ##
##    x2 <- res2%>%filter(comb==ii)
##    th0 <- quantile(abs(x2$beta), probs=0.9)
##    motif2 <- x2%>%filter(qval<0.1, abs(beta)>th0)%>%pull(gene)%>%unique()    

##    summ <- data.frame(comb=ii, nmotif_new=length(motif), nmotif_old=length(motif2),
##                       nolap=length(intersect(motif, motif2)))
##    summ
## })

## opfn <- paste(outdir, "0_compare_new_old.xlsx", sep="")
## write.xlsx(df_olap, file=opfn, overwrite=T)




### motif activities
## fn <- "./sc_multiome_data/3_motif/2_motif.activities.outs/2_motif.ave_macs2_0.1_cn_control.rds"
## mat_motif <- read_rds(fn)
## rn <- rownames(mat_motif)
## rn2 <- motifs[rn]

### motif gene expression


## fn <- "./sc_multiome_data/2_Differential/1_DiffRNA_2.outs/1_YtX.sel_clusters_0.1_cn.rds"
## mat_RNA <- read_rds(fn)
## sum(rn2%in%rownames(mat_RNA))
