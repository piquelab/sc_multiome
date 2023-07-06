##
library(Matrix)
library(tidyverse)
library(data.table)
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

outdir <- "./2_diff_plots.outs/"
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)
 


######################
### summary of DEG ###
######################

###
### heatmap of DEGs
## Mohammed excel files
## "../2_Differential/1.2_DiffPeak.outs/2.0_DESeq.results_0.1_cn_sampleID+new_treat.rds
## ../2_Differential/1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn.rds" from 

###
### DE results

fn <- "./sc_multiome_data/2_Differential/1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn_sampleID+new_treat.rds"
res <- read_rds(fn)%>%as.data.frame()%>%mutate(comb=paste(MCls, contrast, sep="_"))

res2 <- res%>%filter(p.adjusted<0.1, abs(estimate)>0.5)%>%
   mutate(direction=ifelse(estimate>0, "Up", "Down")) # sign(estimate)))

summ <- res2%>%group_by(MCls, contrast)%>%summarize(ngene=n(), .groups="drop")
mat <- summ%>%pivot_wider(id_cols=MCls, names_from=contrast, values_from=ngene)

### data for heatmap
mat2 <- as.matrix(mat[, -1])
rownames(mat2) <- mat$MCls
MCls_sort <- sort(rownames(mat2), decreasing=T)
mat2 <- mat2[MCls_sort,]

### setting colors
olap_min <- min(mat2)
olap_max <- max(mat2)
mycol <- colorRamp2(seq(olap_min, olap_max, length.out=20), colorRampPalette(brewer.pal(n=7,name="YlGnBu"))(20))

 
p <- Heatmap(mat2, name="#DEGs", col=mycol, cluster_rows=F, cluster_columns=F,             
    row_names_gp=gpar(fontsize=10), row_names_side="left",      
    column_names_gp=gpar(fontsize=10), column_names_rot=45,
    heatmap_legend_param=list(at=round(seq(olap_min, olap_max, length.out=8),0),
       grid_width=grid::unit(0.4, "cm"), legend_height=grid::unit(6, "cm"),
       title_gp=gpar(fontsize=9), labels_gp=gpar(fontsize=9), direction="horizontal"),
    cell_fun=function(j, i, x, y, width, height, fill){
        grid.text(mat2[i,j], x, y, gp=gpar(fontsize=9))
    })    

figfn <- paste(outdir, "Figure1.1_DEGs.heatmap.png", sep="")
png(figfn, width=500, height=520, res=120)
draw(p, heatmap_legend_side="bottom")
dev.off()



###
### bar plots for top annotation 
summDF <- res2%>%group_by(contrast, direction)%>%summarize(ny=n(),.groups="drop")

## col=list(celltype=c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56", "3_TEM"="blue",
##                    "4_Bcell"="#4daf4a", "5_CD8Naive"="green", "6_Monocyte"="#984ea3", "7_dnT"="black"),
##             treats=)

p2 <- ggplot(summDF, aes(x=contrast, y=ny, fill=contrast, color=contrast, alpha=direction))+
    geom_bar(stat="identity")+
    scale_fill_manual(values=c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
       "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3"), guide="none")+
    scale_color_manual(values=c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
       "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3"), guide="none")+
    scale_alpha_manual(values=c("Down"=0.5, "Up"=1),
         guide=guide_legend(override.aes=list(fill="red", size=1.5)))+
    ylab("DEGs")+                                                     
    theme_bw()+
    theme(legend.position=c(0.8, 0.8),
          legend.title=element_blank(),
          legend.key=element_blank(),
          legend.box.background=element_blank(),
          legend.background=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_text(size=10),
          axis.text.x=element_blank(),
          axis.text.y=element_text(size=10),
          axis.ticks.x=element_blank())

figfn <- paste(outdir, "Figure1.2_DEGs.bar.png", sep="")
ggsave(figfn, plot=p2, device="png", width=550, height=300, units="px", dpi=120) ## ggsave 


## figfn <- paste(outdir, "Figure1.2_bar_DEG.png", sep="")
## png(figfn, width=550, height=380, res=120)
## p2
## dev.off()


###
### bar plots for row annotation, group_by, cell type

summDF2 <- res2%>%group_by(MCls, direction)%>%summarize(ny=n(),.groups="drop")

p3 <- ggplot(summDF2, aes(x=MCls, y=ny, fill=MCls, color=MCls, alpha=direction))+
    geom_bar(stat="identity")+
    coord_flip()+
    ## scale_y_reverse()+ y axis in the other side
    scale_fill_manual(values=c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
        "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
        "6_Monocyte"="#984ea3", "7_dnT"="black"), guide="none")+
    scale_color_manual(values=c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
        "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
        "6_Monocyte"="#984ea3", "7_dnT"="black"), guide="none")+
    scale_alpha_manual(values=c("Down"=0.5, "Up"=1),
         guide=guide_legend(override.aes=list(fill="#ffaa00", size=1.5)))+
    ylab("DEGs")+                                                     
    theme_bw()+
    theme(legend.position="none",
          legend.title=element_blank(),
          legend.key=element_blank(),
          legend.box.background=element_blank(),
          legend.background=element_blank(),
          axis.title.x=element_text(size=10),
          axis.title.y=element_blank(),
          axis.text.x=element_text(size=8),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
 
figfn <- paste(outdir, "Figure1.3_DEGs.bar.png", sep="")
ggsave(figfn, plot=p3, device="png", width=250, height=420, units="px", dpi=120) ## ggsave 

## figfn <- paste(outdir, "Figure1.3_bar_DEG.png", sep="")
## png(figfn, width=320, height=450, res=120)
## p3
## dev.off()


###############
### heatmap ###
###############


fn <- "./sc_multiome_data/2_Differential/1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn_sampleID+new_treat.rds"
res <- read_rds(fn)%>%as.data.frame()%>%mutate(comb=paste(MCls, contrast, sep="_"))

sigs <- res%>%filter(p.adjusted<0.1, abs(estimate)>0.5)%>%pull(gene)%>%unique()

res2 <- res%>%filter(gene%in%sigs)%>%dplyr::select(comb, gene, estimate, statistic)

###
### data format #DEG*combination
### Orginal data is a matrix of 7241*48 and after filtering missing value, it is 7117*48
mat <- res2%>%pivot_wider(id_cols=gene, names_from=comb, values_from=estimate)
mat2 <- as.matrix(mat[,-1])
rnz <- rowSums(is.na(mat2))
mat2 <- mat2[rnz==0,] ### 7117*48

### correlation data for plot
corr_mat <- cor(mat2, method="spearman") 
mycol <- colorRamp2(seq(-1, 1, length.out=100), colorRampPalette(rev(brewer.pal(n=7,name="RdBu")))(100))


###
### annotation columns
col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
  "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
   "6_Monocyte"="#984ea3", "7_dnT"="black")
col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
       "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3")

x <- str_split(colnames(corr_mat), "_", simplify=T)
df_col <- data.frame(celltype=paste(x[,1], x[,2], sep="_"), contrast=x[,3])
col_ha <- HeatmapAnnotation(df=df_col, col=list(celltype=col1, contrast=col2),
    annotation_name_gp=gpar(fontsize=10),
    annotation_legend_param=list(
  celltype=list(title_gp=gpar(fontsize=10), labels_gp=gpar(fontsize=8), title="celltype",
                grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm")),
  contrast=list(title_gp=gpar(fontsize=10), labels_gp=gpar(fontsize=8), title="contrast",
                grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm")))) 


###
### main heatmap  
p4 <- Heatmap(corr_mat, col=mycol,
   cluster_rows=F, cluster_columns=F,
   show_row_names=F,
   show_column_names=T, column_names_gp=gpar(fontsize=6),
   ##column_names_rot=-45,   
   top_annotation=col_ha,
   heatmap_legend_param=list(title="rho",
      title_gp=gpar(fontsize=10),
      at=seq(-1, 1, by=0.25), 
      labels_gp=gpar(fontsize=8),
      grid_width=grid::unit(0.45, "cm"),
      legend_height=grid::unit(7.8, "cm")), 
   use_raster=T, raster_device="png")
 
###
figfn <- paste(outdir, "Figure1.4_corr_LFC_spearman.heatmap.png", sep="")
## ggsave(figfn, plot=p2, device="png", width=680, height=700, units="px", dpi=300) ## ggsave 
png(figfn, width=700, height=780,res=120)
p4 <- draw(p4)
dev.off()



#################################
### boxplots show correlation ###
#################################

###
### generate pairwise correlation matrix
fn <- "./sc_multiome_data/2_Differential/1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn_sampleID+new_treat.rds"
res <- read_rds(fn)%>%as.data.frame()%>%mutate(comb=paste(MCls, contrast, sep="_"))
sigs <- res%>%filter(p.adjusted<0.1, abs(estimate)>0.5)%>%pull(gene)%>%unique()
res2 <- res%>%filter(gene%in%sigs)%>%dplyr::select(comb, gene, estimate, statistic)
 
###
### data format #DEG*combination 
### Orginal data is a matrix of 7241*48 and after filtering missing value, it is 7117*48
mat <- res2%>%pivot_wider(id_cols=gene, names_from=comb, values_from=estimate)
mat2 <- as.matrix(mat[,-1])
rnz <- rowSums(is.na(mat2))
mat2 <- mat2[rnz==0,] ### 7117*48

### correlation data for plot
## corr_mat <- cor(mat2, method="spearman") 
## opfn <- paste(outdir, "1.x_DEG_corr.rds", sep="")
## write_rds(corr_mat, file=opfn)

###
### setting colors
col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
  "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
   "6_Monocyte"="#984ea3", "7_dnT"="black")
col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
       "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3")



###
### plot data group by cell types
rn<- rownames(corr_mat)
cvt <- str_split(rownames(corr_mat), "_", simplify=T)
cvt2 <- data.frame(rn=rownames(corr_mat), MCls=paste(cvt[,1], cvt[,2], sep="_"), treats=cvt[,3])

###
MCls <- unique(cvt2$MCls)
DFcorr <- map_dfr(MCls, function(ii){
   ###
   rn2 <- cvt2%>%filter(MCls==ii)%>%pull(rn) 
   corr2 <- corr_mat[rn2, rn2]
   corr2 <- corr2%>%as.data.frame()%>%rownames_to_column(var="comb") 
   ## 
   DF2 <- corr2%>%pivot_longer(!comb, names_to="comb2", values_to="rr")  ##%>%filter(!rr==1)
   DF2$MCls <- ii
   DF2 
})
 

###
treats <- unique(cvt2$treats)
DFcorr2 <- map_dfr(treats, function(ii){
   ###
   rn2 <- cvt2%>%filter(treats==ii)%>%pull(rn) 
   corr2 <- corr_mat[rn2, rn2]
   corr2 <- corr2%>%as.data.frame()%>%rownames_to_column(var="comb") 
   ## 
   DF2 <- corr2%>%pivot_longer(!comb, names_to="comb2", values_to="rr") ##%>%filter(!rr==1)
   DF2$treats <- ii
   DF2 
})


###
### boxplot group by cell types
DFcorr <- DFcorr%>%mutate(rr=round(rr,3))%>%filter(rr<1)
p5_1 <- ggplot(DFcorr, aes(x=MCls, y=rr, fill=MCls))+
   ##geom_violin()+ 
   geom_boxplot(outlier.shape=NA, lwd=0.25)+
   stat_summary(fun=median, color="grey", geom="point", shape=23, size=2.5, stroke=0.9)+
   geom_jitter(width=0.2, size=0.1)+     
   scale_fill_manual(values=col1)+
   ylab("Spearman correlation")+ylim(0, 0.8)+
   theme_bw()+
   theme(legend.position="none",
         axis.title.x=element_blank(),
         axis.title.y=element_text(size=9),
         axis.text.x=element_text(angle=45, hjust=1, size=9),
         axis.text.y=element_text(size=9))                         


### boxplots group by treatments
DFcorr2 <- DFcorr2%>%mutate(rr=round(rr,3))%>%filter(rr<1)
p5_2 <- ggplot(DFcorr2, aes(x=treats, y=rr, fill=treats))+
   ##geom_violin()+ 
   geom_boxplot(outlier.shape=NA, lwd=0.25)+
   stat_summary(fun=median, color="grey", geom="point", shape=23, size=2.5, stroke=0.9)+
   geom_jitter(width=0.2, size=0.1)+ 
   scale_fill_manual(values=col2)+
   ylab("Spearman correlation")+ylim(0,0.8)+
   theme_bw()+
   theme(legend.position="none",
         axis.title.x=element_blank(),
         axis.title.y=element_text(size=9),
         axis.text.x=element_text(angle=45, hjust=1, size=9),
         axis.text.y=element_text(size=9))                         


 
comb_plots <- plot_grid(p5_1, p5_2, nrow=2, align="v", rel_heights=c(1,0.9))
figfn <- paste(outdir, "Figure1.5_comb.box.png", sep="")
ggsave(figfn, plot=comb_plots, device="png", width=320, height=560, units="px", dpi=120) ## ggsave 



###
### pair treatments within the same cell types, each box has 8 points
mat <- read_rds("./2_diff_plots.outs/1.x_DEG_corr.rds")
x <- str_split(colnames(mat), "_", simplify=T)
df <- data.frame(rn=colnames(mat), MCls=paste(x[,1], x[,2], sep="_"), treats=x[,3])

treat <- sort(unique(df$treats))
MCls <- sort(unique(df$MCls))
ntr <- length(treat)
comb <- NULL
for (i in 1:(ntr-1)){
   ##
   for (j in (i+1):ntr){
       ii <- paste(treat[i], treat[j], sep="_")
       comb <- c(comb, ii)
   }
}

       
DFcorr3 <- map_dfr(comb, function(ii){
   ##
   tr1 <- gsub("_.*", "", ii)
   tr2 <- gsub(".*_", "", ii) 
   df <- map_dfr(MCls, function(oneMCl){
      ##
      comb1 <- paste(oneMCl, tr1, sep="_")
      comb2 <- paste(oneMCl, tr2, sep="_")
      df2 <- data.frame(comb=ii, MCls=oneMCl, rr=mat[comb1, comb2])
      df2
   })
   df
})    
 
p5_3 <- ggplot(DFcorr3, aes(x=comb, y=as.numeric(rr)))+ ##, fill=factor(MCls), group=comb))+
   ## geom_violin()+
   geom_boxplot(outlier.shape=NA, lwd=0.25)+
   stat_summary(fun=median, color="grey", geom="point", shape=23, size=1.8, stroke=0.9)+
   geom_jitter(width=0.2, size=0.5, aes(color=MCls))+ 
   scale_color_manual(values=col1)+
   ylab("Spearman correlation")+ylim(0,0.8)+
   theme_bw()+
   theme(legend.position="none",
         axis.title.x=element_blank(),
         axis.title.y=element_text(size=9),
         axis.text.x=element_text(angle=45, hjust=1, size=8),
         axis.text.y=element_text(size=9))                         
  
figfn <- paste(outdir, "Figure1.6_pairtreats.box.png", sep="")
ggsave(figfn, plot=p5_3, device="png", width=460, height=380, units="px", dpi=120)




#################################################
### summary Differentially accessible regions ###
#################################################

rm(list=ls())

outdir <- "./2_diff_plots.outs/"

###
### DE results

fn <- "./sc_multiome_data/2_Differential/1.2_DiffPeak.outs/2.0_DESeq.results_0.1_cn_sampleID+new_treat.rds"
res <- read_rds(fn)%>%as.data.frame()%>%mutate(comb=paste(MCls, contrast, sep="_"))

res2 <- res%>%filter(p.adjusted<0.1, abs(estimate)>0.5)%>%
   mutate(direction=ifelse(estimate>0, "Up", "Down")) # sign(estimate)))

sigs <- unique(res2$gene)

summ <- res2%>%group_by(MCls, contrast)%>%summarize(ngene=n(), .groups="drop")
mat <- summ%>%pivot_wider(id_cols=MCls, names_from=contrast, values_from=ngene)

### data for heatmap
mat2 <- as.matrix(mat[, -1])
rownames(mat2) <- mat$MCls
MCls_sort <- sort(rownames(mat2), decreasing=T)
mat2 <- mat2[MCls_sort,]

### setting colors
olap_min <- min(mat2)
olap_max <- max(mat2)
mycol <- colorRamp2(seq(olap_min, olap_max, length.out=20), colorRampPalette(brewer.pal(n=7,name="YlGnBu"))(20))

 
p <- Heatmap(mat2, name="#DARs", col=mycol, cluster_rows=F, cluster_columns=F,             
    row_names_gp=gpar(fontsize=10), row_names_side="left",      
    column_names_gp=gpar(fontsize=10), column_names_rot=45,
    heatmap_legend_param=list(at=round(seq(olap_min, olap_max, length.out=8),0),
       grid_width=grid::unit(0.4, "cm"), legend_height=grid::unit(6, "cm"),
       title_gp=gpar(fontsize=9), labels_gp=gpar(fontsize=9), direction="horizontal"),
    cell_fun=function(j, i, x, y, width, height, fill){
        grid.text(mat2[i,j], x, y, gp=gpar(fontsize=9))
    })    


figfn <- paste(outdir, "Figure2.1_DARs.heatmap.png", sep="")
png(figfn, width=500, height=520, res=120)
draw(p, heatmap_legend_side="bottom")
dev.off()


## ggsave(figfn, plot=p, device="png", width=580, height=620, units="px", dpi=300) ## ggsave 
## unlink(figfn)


###
### bar plots for top annotation 
summDF <- res2%>%group_by(contrast, direction)%>%summarize(ny=n(),.groups="drop")

## col=list(celltype=c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56", "3_TEM"="blue",
##                    "4_Bcell"="#4daf4a", "5_CD8Naive"="green", "6_Monocyte"="#984ea3", "7_dnT"="black"),
##             treats=)

p2 <- ggplot(summDF, aes(x=contrast, y=ny, fill=contrast, color=contrast, alpha=direction))+
    geom_bar(stat="identity")+
    scale_fill_manual(values=c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
       "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3"), guide="none")+
    scale_color_manual(values=c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
       "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3"), guide="none")+
    scale_alpha_manual(values=c("Down"=0.5, "Up"=1),
         guide=guide_legend(override.aes=list(fill="red", size=1.5)))+
    ylab("DARs")+                                                     
    theme_bw()+
    theme(legend.position=c(0.8, 0.8),
          legend.title=element_blank(),
          legend.key=element_blank(),
          legend.box.background=element_blank(),
          legend.background=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_text(size=10),
          axis.text.x=element_blank(),
          axis.text.y=element_text(size=10),
          axis.ticks.x=element_blank())
##
##
figfn <- paste(outdir, "Figure2.2_DARs.bar.png", sep="")
ggsave(figfn, plot=p2, device="png", width=550, height=300, units="px", dpi=120) ## ggsave 
##unlink(figfn)

## figfn <- paste(outdir, "Figure2.2_bar_DEG.png", sep="")
## png(figfn, width=550, height=380, res=120)
## p2
## dev.off()


###
### bar plots for row annotation, group_by, cell type

summDF2 <- res2%>%group_by(MCls, direction)%>%summarize(ny=n(),.groups="drop")


p3 <- ggplot(summDF2, aes(x=MCls, y=ny, fill=MCls, color=MCls, alpha=direction))+
    geom_bar(stat="identity")+
    coord_flip()+
    ##scale_y_reverse()+
    scale_fill_manual(values=c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
        "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
        "6_Monocyte"="#984ea3", "7_dnT"="black"), guide="none")+
    scale_color_manual(values=c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
        "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
        "6_Monocyte"="#984ea3", "7_dnT"="black"), guide="none")+
    scale_alpha_manual(values=c("Down"=0.5, "Up"=1),
         guide=guide_legend(override.aes=list(fill="#ffaa00", size=1.5)))+
    ylab("DARs")+                                                     
    theme_bw()+
    theme(legend.position="none",
          legend.title=element_blank(),
          legend.key=element_blank(),
          legend.box.background=element_blank(),
          legend.background=element_blank(),
          axis.title.x=element_text(size=10),
          axis.title.y=element_blank(),
          axis.text.x=element_text(size=10),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())

figfn <- paste(outdir, "Figure2.3_DARs.bar.png", sep="")
ggsave(figfn, plot=p3, device="png", width=250, height=420, units="px", dpi=120) ## ggsave 

## figfn <- paste(outdir, "Figure1.3_bar_DEG.png", sep="")
## png(figfn, width=320, height=450, res=120)
## p3
## dev.off()



######################################################
### Heatmap of correlation using spearman with LFC ### 
######################################################

## fn <- "./sc_multiome_data/2_Differential/1.2_DiffPeak.outs/2.0_DESeq.results_0.1_cn_sampleID+new_treat.rds"
## res <- read_rds(fn)%>%as.data.frame()%>%mutate(comb=paste(MCls, contrast, sep="_"))
## sigs <- res%>%filter(p.adjusted<0.1, abs(estimate)>0.5)%>%pull(gene)%>%unique()
## res2 <- res%>%filter(gene%in%sigs)%>%dplyr::select(comb, gene, estimate, statistic)


## ###
## mat <- res2%>%pivot_wider(id_cols=gene, names_from=comb, values_from=estimate)
## mat2 <- as.matrix(mat[,-1])
## rnz <- rowSums(is.na(mat2))
## mat2 <- mat2[rnz==0,] ### 17,445*48

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
##     show_legend=F,                        
##     annotation_name_gp=gpar(fontsize=10),
##     annotation_legend_param=list(
##   celltype=list(title_gp=gpar(fontsize=10), labels_gp=gpar(fontsize=8), title="celltype",
##                 grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm")),
##   contrast=list(title_gp=gpar(fontsize=10), labels_gp=gpar(fontsize=8), title="contrast",
##                 grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm")))) 


## ###
## ### main heatmap  
## p4 <- Heatmap(corr_mat, col=mycol,
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
## figfn <- paste(outdir, "Figure2.4_corr_LFC_spearman.heatmap.png", sep="")
## ## ggsave(figfn, plot=p2, device="png", width=680, height=700, units="px", dpi=300) ## ggsave 
## png(figfn, width=580, height=700,res=120)
## p4 <- draw(p4)
## dev.off()




###################################################################
#### boxplots show spearman correlation of LFC of union of DARs ###
###################################################################
## fn <- "./sc_multiome_data/2_Differential/1.2_DiffPeak.outs/2.0_DESeq.results_0.1_cn_sampleID+new_treat.rds"
## res <- read_rds(fn)%>%as.data.frame()%>%mutate(comb=paste(MCls, contrast, sep="_"))
## sigs <- res%>%filter(p.adjusted<0.1, abs(estimate)>0.5)%>%pull(gene)%>%unique()
## res2 <- res%>%filter(gene%in%sigs)%>%dplyr::select(comb, gene, estimate, statistic)



## mat <- res2%>%pivot_wider(id_cols=gene, names_from=comb, values_from=estimate)
## mat2 <- as.matrix(mat[,-1])
## rnz <- rowSums(is.na(mat2))
## mat2 <- mat2[rnz==0,] ### 17,445*48

## ### correlation data for plot
## corr_mat <- cor(mat2, method="spearman") 



## ###
## ### setting colors
## col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
##   "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
##    "6_Monocyte"="#984ea3", "7_dnT"="black")
## col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
##        "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3")



## ###
## ### plot data group by cell types
## rn<- rownames(corr_mat)
## cvt <- str_split(rownames(corr_mat), "_", simplify=T)
## cvt2 <- data.frame(rn=rownames(corr_mat), MCls=paste(cvt[,1], cvt[,2], sep="_"), treats=cvt[,3])

## ###
## MCls <- unique(cvt2$MCls)
## DFcorr <- map_dfr(MCls, function(ii){
##    ###
##    rn2 <- cvt2%>%filter(MCls==ii)%>%pull(rn) 
##    corr2 <- corr_mat[rn2, rn2]
##    corr2 <- corr2%>%as.data.frame()%>%rownames_to_column(var="comb") 
##    ## 
##    DF2 <- corr2%>%pivot_longer(!comb, names_to="comb2", values_to="rr")
##    DF2$MCls <- ii
##    DF2 
## })
 

## ###
## treats <- unique(cvt2$treats)
## DFcorr2 <- map_dfr(treats, function(ii){
##    ###
##    rn2 <- cvt2%>%filter(treats==ii)%>%pull(rn) 
##    corr2 <- corr_mat[rn2, rn2]
##    corr2 <- corr2%>%as.data.frame()%>%rownames_to_column(var="comb") 
##    ## 
##    DF2 <- corr2%>%pivot_longer(!comb, names_to="comb2", values_to="rr")
##    DF2$treats <- ii
##    DF2 
## })


## ###
## ### boxplot group by cell types
## DFcorr <- DFcorr%>%mutate(rr=round(rr,3))%>%filter(rr<1)
## p5_1 <- ggplot(DFcorr, aes(x=MCls, y=rr, fill=MCls))+
##    ##geom_violin()+ 
##    geom_boxplot(outlier.shape=NA, lwd=0.25)+
##    stat_summary(fun=median, color="grey", geom="point", shape=23, size=2.5, stroke=0.9)+
##    geom_jitter(width=0.2, size=0.1)+     
##    scale_fill_manual(values=col1)+
##    ylab("Spearman correlation")+ylim(0, 0.8)+   
##    theme_bw()+
##    theme(legend.position="none",
##          axis.title.x=element_blank(),
##          axis.title.y=element_text(size=9),
##          axis.text.x=element_text(angle=45, hjust=1, size=9),
##          axis.text.y=element_text(size=9))                         


## ### boxplots group by treatments
## DFcorr2 <- DFcorr2%>%mutate(rr=round(rr,3))%>%filter(rr<1)
## p5_2 <- ggplot(DFcorr2, aes(x=treats, y=rr, fill=treats))+
##    ##geom_violin()+ 
##    geom_boxplot(outlier.shape=NA, lwd=0.25)+
##    stat_summary(fun=median, color="grey", geom="point", shape=23, size=2.5, stroke=0.9)+
##    geom_jitter(width=0.2, size=0.1)+ 
##    scale_fill_manual(values=col2)+
##    ylab("Spearman correlation")+ylim(0, 0.8)+
##    theme_bw()+
##    theme(legend.position="none",
##          axis.title.x=element_blank(),
##          axis.title.y=element_text(size=9),
##          axis.text.x=element_text(angle=45, hjust=1, size=9),
##          axis.text.y=element_text(size=9))                         


## comb_plots <- plot_grid(p5_1, p5_2, nrow=2, align="v", rel_heights=c(1,0.9))
## figfn <- paste(outdir, "Figure2.5_comb.box.png", sep="")
## ggsave(figfn, plot=comb_plots, device="png", width=320, height=560, units="px", dpi=120) ## ggsave 


 

#################################################
#### another option to calculate correlation ####
#################################################
rm(list=ls())
outdir <- "./2_diff_plots.outs/"

fn <- "./sc_multiome_data/2_Differential/1.2_DiffPeak.outs/2.0_DESeq.results_0.1_cn_sampleID+new_treat.rds"
res <- read_rds(fn)%>%as.data.frame()%>%mutate(comb=paste(MCls, contrast, sep="_"))

## sig <- res%>%filter(p.adjusted<0.1, abs(estimate)>0.5)%>%pull(gene)%>%unique()

###
### peaks
MCls <- unique(res$MCls)
DFpeak <- map_dfr(MCls, function(ii){
   ###
   peakSel <- res%>%filter(p.adjusted<0.1, abs(estimate)>0.5, MCls==ii)%>%pull(gene)%>%unique() 
   DF <- data.frame(peak=peakSel)
   DF$MCls <- ii
   DF
})    

## summ <- DFpeak%>%group_by(MCls)%>%summarize(npeak=n(), .groups="drop")
## write.xlsx(summ, "./2_diff_plots.outs/2_summ.xlsx")


###
### correlation 
comb <- sort(unique(res$comb))
ncomb <- length(comb)
cvt <- str_split(comb, "_", simplify=T)
DFcomb <- data.frame(rn=comb, MCls=paste(cvt[,1], cvt[,2], sep="_"), treats=cvt[,3])


corr_mat <- matrix(0, ncomb, ncomb)
colnames(corr_mat) <- comb
rownames(corr_mat) <- comb
MCls <- sort(unique(DFcomb$MCls))
for (oneMCl in MCls){
   ###
   peakSel <- DFpeak%>%filter(MCls==oneMCl)%>%pull(peak)%>%unique()
   rn2 <- DFcomb%>%filter(MCls==oneMCl)%>%pull(rn)%>%unique()%>%sort()

   res2 <- res%>%filter(comb%in%rn2, gene%in%peakSel)%>%select(comb, gene, LFC=estimate)
   dfmat <- res2%>%pivot_wider(id_cols=gene, names_from=comb, values_from=LFC)
   rr_mat  <- cor(dfmat[,-1], method="spearman")
   corr_mat[rn2, rn2] <- rr_mat[rn2, rn2]     
   cat(oneMCl, "\n")
}    
     
opfn <- paste(outdir, "2.x_DP_corr.rds", sep="")
write_rds(corr_mat, file=opfn)




###
### Heatmap 
fn <- paste(outdir, "2.x_DP_corr.rds", sep="")
corr_mat <- read_rds(fn)




### setting colors
mycol <- colorRamp2(seq(-1, 1, length.out=100), colorRampPalette(rev(brewer.pal(n=7,name="RdBu")))(100))


###
### annotation columns
col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
  "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
   "6_Monocyte"="#984ea3", "7_dnT"="black")
col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
       "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3")

x <- str_split(colnames(corr_mat), "_", simplify=T)
df_col <- data.frame(celltype=paste(x[,1], x[,2], sep="_"), contrast=x[,3])
col_ha <- HeatmapAnnotation(df=df_col, col=list(celltype=col1, contrast=col2),
    show_legend=F,                        
    annotation_name_gp=gpar(fontsize=10),
    annotation_legend_param=list(
  celltype=list(title_gp=gpar(fontsize=10), labels_gp=gpar(fontsize=8), title="celltype",
                grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm")),
  contrast=list(title_gp=gpar(fontsize=10), labels_gp=gpar(fontsize=8), title="contrast",
                grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm")))) 


###
### main heatmap  
p2 <- Heatmap(corr_mat, col=mycol,
   cluster_rows=F, cluster_columns=F,
   show_row_names=F,
   show_column_names=T, column_names_gp=gpar(fontsize=6),
   ##column_names_rot=-45,   
   top_annotation=col_ha,
   heatmap_legend_param=list(title="rho",
      title_gp=gpar(fontsize=10),
      at=seq(-1, 1, by=0.25), 
      labels_gp=gpar(fontsize=8),
      grid_width=grid::unit(0.45, "cm"),
      legend_height=grid::unit(6, "cm")), 
   use_raster=T, raster_device="png")

###
figfn <- paste(outdir, "Figure2.4_corr_LFC_spearman.heatmap.png", sep="")
## ggsave(figfn, plot=p2, device="png", width=680, height=700, units="px", dpi=300) ## ggsave 
png(figfn, width=580, height=700,res=120)
p2 <- draw(p2)
dev.off()




###################################################################
#### boxplots show spearman correlation of LFC of union of DARs ###
###################################################################

rm(list=ls())

outdir <- "./2_diff_plots.outs/"

corr_mat <- read_rds("./2_diff_plots.outs/2.x_DP_corr.rds")
###
### setting colors
col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
  "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
   "6_Monocyte"="#984ea3", "7_dnT"="black")
col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
       "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3")



###
### plot data group by cell types
rn<- rownames(corr_mat)
cvt <- str_split(rownames(corr_mat), "_", simplify=T)
cvt2 <- data.frame(rn=rownames(corr_mat), MCls=paste(cvt[,1], cvt[,2], sep="_"), treats=cvt[,3])

###
MCls <- unique(cvt2$MCls)
DFcorr <- map_dfr(MCls, function(ii){
   ###
   rn2 <- cvt2%>%filter(MCls==ii)%>%pull(rn) 
   corr2 <- corr_mat[rn2, rn2]
   corr2 <- corr2%>%as.data.frame()%>%rownames_to_column(var="comb") 
   ## 
   DF2 <- corr2%>%pivot_longer(!comb, names_to="comb2", values_to="rr")%>%
       mutate(rr=round(rr,3))%>%filter(rr<1)
   DF2$MCls <- ii
   DF2 
})
 

###
## treats <- unique(cvt2$treats)
## DFcorr2 <- map_dfr(treats, function(ii){
##    ###
##    rn2 <- cvt2%>%filter(treats==ii)%>%pull(rn) 
##    corr2 <- corr_mat[rn2, rn2]
##    corr2 <- corr2%>%as.data.frame()%>%rownames_to_column(var="comb") 
##    ## 
##    DF2 <- corr2%>%pivot_longer(!comb, names_to="comb2", values_to="rr")
##    DF2$treats <- ii
##    DF2 
## })


 

###
### boxplot group by cell types
DFcorr <- DFcorr%>%mutate(rr=round(rr,3))%>%filter(rr<1)
p6_1 <- ggplot(DFcorr, aes(x=MCls, y=rr, fill=MCls))+
   ##geom_violin()+ 
   geom_boxplot(outlier.shape=NA, lwd=0.25)+
   stat_summary(fun=median, color="grey", geom="point", shape=23, size=2.5, stroke=0.9)+
   geom_jitter(width=0.2, size=0.1)+     
   scale_fill_manual(values=col1)+
   ylab("Spearman correlation")+ylim(-0.1,0.8)+
   theme_bw()+
   theme(legend.position="none",
         axis.title.x=element_blank(),
         axis.title.y=element_text(size=9),
         axis.text.x=element_text(angle=45, hjust=1, size=9),
         axis.text.y=element_text(size=9))                         

###
figfn <- paste(outdir, "Figure2.5_comb.box.png", sep="")
ggsave(figfn, plot=p6_1, device="png", width=350, height=320, units="px", dpi=120) ## ggsave 



###
### pair treatments within the same cell types, each box has 8 points
mat <- read_rds("./2_diff_plots.outs/2.x_DP_corr.rds")
x <- str_split(colnames(mat), "_", simplify=T)
df <- data.frame(rn=colnames(mat), MCls=paste(x[,1], x[,2], sep="_"), treats=x[,3])

treat <- sort(unique(df$treats))
MCls <- sort(unique(df$MCls))
ntr <- length(treat)
comb <- NULL
for (i in 1:(ntr-1)){
   ##
   for (j in (i+1):ntr){
       ii <- paste(treat[i], treat[j], sep="_")
       comb <- c(comb, ii)
   }
}

       
DFcorr3 <- map_dfr(comb, function(ii){
   ##
   tr1 <- gsub("_.*", "", ii)
   tr2 <- gsub(".*_", "", ii) 
   df <- map_dfr(MCls, function(oneMCl){
      ##
      comb1 <- paste(oneMCl, tr1, sep="_")
      comb2 <- paste(oneMCl, tr2, sep="_")
      df2 <- data.frame(comb=ii, MCls=oneMCl, rr=mat[comb1, comb2])
      df2
   })
   df
})    
 
p6_3 <- ggplot(DFcorr3, aes(x=comb, y=as.numeric(rr)))+ ##, fill=factor(MCls), group=comb))+
   ## geom_violin()+
   geom_boxplot(outlier.shape=NA, lwd=0.25)+
   stat_summary(fun=median, color="grey", geom="point", shape=23, size=1.8, stroke=0.9)+
   geom_jitter(width=0.2, size=0.5, aes(color=MCls))+ 
   scale_color_manual(values=col1)+
   ylab("Spearman correlation")+ylim(-0.1,0.8)+
   theme_bw()+
   theme(legend.position="none",
         axis.title.x=element_blank(),
         axis.title.y=element_text(size=9),
         axis.text.x=element_text(angle=45, hjust=1, size=7.5),
         axis.text.y=element_text(size=9))                         
  
figfn <- paste(outdir, "Figure2.6_pairtreats.box.png", sep="")
ggsave(figfn, plot=p6_3, device="png", width=460, height=380, units="px", dpi=120)





### boxplots group by treatments
## DFcorr2 <- DFcorr2%>%mutate(rr=round(rr,3))%>%filter(rr<1)
## p6_2 <- ggplot(DFcorr2, aes(x=treats, y=rr, fill=treats))+
##    ##geom_violin()+ 
##    geom_boxplot(outlier.shape=NA, lwd=0.25)+
##    stat_summary(fun=median, color="grey", geom="point", shape=23, size=2.5, stroke=0.9)+
##    geom_jitter(width=0.2, size=0.1)+ 
##    scale_fill_manual(values=col2)+
##    ylab("Spearman correlation")+ylim(-0.2,0.8)+
##    theme_bw()+
##    theme(legend.position="none",
##          axis.title.x=element_blank(),
##          axis.title.y=element_text(size=9),
##          axis.text.x=element_text(angle=45, hjust=1, size=9),
##          axis.text.y=element_text(size=9))                         


## comb_plots <- plot_grid(p6_1, p6_2, nrow=2, align="v", rel_heights=c(1,0.9))
## figfn <- paste(outdir, "Figure2.6.2_comb.box.png", sep="")
## ggsave(figfn, plot=comb_plots, device="png", width=320, height=560, units="px", dpi=120) ## ggsave 
