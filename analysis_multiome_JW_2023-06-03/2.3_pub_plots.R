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

outdir <- "./Plots_pub/2_diff_plots/"
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)
 
###
### script for doing plots in the manuscript  


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

###
### setting colors
col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
  "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
   "6_Monocyte"="#984ea3", "7_dnT"="black")
col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
       "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3")





fn <- "./sc_multiome_data/2_Differential/1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn_sampleID+new_treat.rds"
res <- read_rds(fn)%>%as.data.frame()%>%mutate(comb=paste(MCls, contrast, sep="_"))

res2 <- res%>%filter(p.adjusted<0.1, abs(estimate)>0.5)%>%
   mutate(direction=ifelse(estimate>0, "Up", "Down")) # sign(estimate)))

summ <- res2%>%group_by(MCls, contrast)%>%summarize(ngene=n(), .groups="drop")
mat <- summ%>%pivot_wider(id_cols=MCls, names_from=contrast, values_from=ngene)


res2%>%group_by(contrast)%>%summarize(ngene=length(unique(gene)), .groups="drop")%>%arrange(desc(ngene))
res2%>%group_by(MCls)%>%summarize(ngene=length(unique(gene)), .groups="drop")%>%arrange(desc(ngene))

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
    row_names_gp=gpar(fontsize=12), row_names_side="left",      
    column_names_gp=gpar(fontsize=12), column_names_rot=45,
    heatmap_legend_param=list(at=round(seq(olap_min, olap_max, length.out=5),0),
       grid_width=grid::unit(0.4, "cm"), legend_height=grid::unit(6, "cm"), title="",
       title_gp=gpar(fontsize=9), labels_gp=gpar(fontsize=12), direction="horizontal"),
    cell_fun=function(j, i, x, y, width, height, fill){
        grid.text(mat2[i,j], x, y, gp=gpar(fontsize=12))
    })    

figfn <- paste(outdir, "Figure1.1_DEGs.heatmap.png", sep="")
png(figfn, width=500, height=520, res=120)
draw(p, heatmap_legend_side="bottom")
dev.off()



###
### bar plots for top annotation for treatment

x <- res2%>%group_by(contrast, direction)%>%summarize(ny=length(unique(gene)), .groups="drop")%>%
    group_by(contrast)%>%summarize(ny=sum(ny),.groups="drop")

summDF <- res2%>%group_by(contrast)%>%summarize(ny=length(unique(gene)),.groups="drop")


p2 <- ggplot(summDF, aes(x=contrast, y=ny, fill=contrast))+ ## , color=contrast))+ ##, alpha=direction))+
    geom_bar(stat="identity")+
    scale_fill_manual(values=col2, guide="none")+
    ## scale_color_manual(values=c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
    ##    "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3"), guide="none")+
    ## scale_alpha_manual(values=c("Down"=0.5, "Up"=1),
    ##      guide=guide_legend(override.aes=list(fill="red", size=1.5)))+
    ylab("DEGs")+                                                     
    theme_bw()+
    theme(legend.position="none", ##c(0.8, 0.8),
          legend.title=element_blank(),
          legend.key=element_blank(),
          legend.box.background=element_blank(),
          legend.background=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_text(size=12),
          axis.text.x=element_blank(),
          axis.text.y=element_text(size=12),
          axis.ticks.x=element_blank())

figfn <- paste(outdir, "Figure1.1_DEGs_treat.bar.png", sep="")
ggsave(figfn, plot=p2, width=480, height=220, units="px", dpi=120) ## ggsave 


## figfn <- paste(outdir, "Figure1.2_bar_DEG.png", sep="")
## png(figfn, width=550, height=380, res=120)
## p2
## dev.off()


###
### bar plots for row annotation, group_by, cell type
## x <- res2%>%group_by(MCls, direction)%>%summarize(ny=length(unique(gene)), .groups="drop")%>%
##     group_by(MCls)%>%summarize(ny=sum(ny),.groups="drop")

summDF2 <- res2%>%group_by(MCls)%>%summarize(ny=length(unique(gene)),.groups="drop")

p3 <- ggplot(summDF2, aes(x=MCls, y=ny, fill=MCls))+ ##, color=MCls, alpha=direction))+
    geom_bar(stat="identity")+
    coord_flip()+
    ## scale_y_reverse()+ y axis in the other side
    scale_fill_manual(values=col1, guide="none")+
    ## scale_color_manual(values=c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
    ##     "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
    ##     "6_Monocyte"="#984ea3", "7_dnT"="black"), guide="none")+
    ## scale_alpha_manual(values=c("Down"=0.5, "Up"=1),
    ##      guide=guide_legend(override.aes=list(fill="#ffaa00", size=1.5)))+
    ylab("DEGs")+                                                     
    theme_bw()+
    theme(legend.position="none",
          legend.title=element_blank(),
          legend.key=element_blank(),
          legend.box.background=element_blank(),
          legend.background=element_blank(),
          axis.title.x=element_text(size=12),
          axis.title.y=element_blank(),
          axis.text.x=element_text(size=12),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
 
figfn <- paste(outdir, "Figure1.1_DEGs_MCls.bar.png", sep="")
ggsave(figfn, plot=p3, width=250, height=480, units="px", dpi=120) ## ggsave 

## figfn <- paste(outdir, "Figure1.3_bar_DEG.png", sep="")
## png(figfn, width=320, height=450, res=120)
## p3
## dev.off()


####
#### correlation between conditions 

###############
### heatmap ###
###############

rm(list=ls())

outdir <- "./Plots_pub/2_diff_plots/"
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)




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

## opfn <- paste(outdir, "1.x_DEG_corr.rds", sep="")
## write_rds(corr_mat, file=opfn)

corr_mat <- read_rds("./Plots_pub/2_diff_plots/1.x_DEG_corr.rds")
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
p0 <- Heatmap(corr_mat, col=mycol,
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
figfn <- paste(outdir, "FigS2_3_LFC.RNA_corr.heatmap.png", sep="")
## ggsave(figfn, plot=p2, device="png", width=680, height=700, units="px", dpi=300) ## ggsave 
png(figfn, width=700, height=780,res=120)
p0 <- draw(p0)
dev.off()



#################################
### boxplots show correlation ###
#################################


rm(list=ls())

outdir <- "./Plots_pub/2_diff_plots/"
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)



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
corr_mat <- cor(mat2, method="spearman") 
opfn <- paste(outdir, "1.x_DEG_corr.rds", sep="")
write_rds(corr_mat, file=opfn)


### read data
fn <- paste(outdir, "1.x_DEG_corr.rds", sep="")
corr_mat <- read_rds(fn)

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
   rnSel <- cvt2%>%filter(MCls==ii)%>%pull(rn)
   ntr <- length(rnSel)
   ### 
   DF2 <- NULL     
   for ( i in 1:(ntr-1)){
      for ( j in (i+1):ntr){
          ##
          rn_i <- rnSel[i]
          rn_j <- rnSel[j]
          rr0 <- round(corr_mat[rn_i, rn_j], 3)
          ###
          df2 <- data.frame(comb=rn_i, comb2=rn_j, rr=rr0)
          DF2 <- rbind(DF2, df2)
      }    
   }
   DF2$MCls <- ii
   DF2 
})
 


###
### boxplot group by cell types
DFcorr <- DFcorr%>%mutate(rr=round(rr,3))%>%filter(rr<1)
p1 <- ggplot(DFcorr, aes(x=MCls, y=rr, fill=MCls))+
   ##geom_violin()+ 
   geom_boxplot(outlier.shape=NA, lwd=0.25)+
   stat_summary(fun=median, color="grey", geom="point", shape=23, size=2.5, stroke=0.9)+
   geom_jitter(width=0.2, size=0.2)+     
   scale_fill_manual(values=col1)+
   ylab("SCC (DEGs)")+ylim(0, 0.8)+
   theme_bw()+
   theme(legend.position="none",
         axis.title.x=element_blank(),
         axis.title.y=element_text(size=12),
         axis.text.x=element_blank(),   ##element_text(angle=45, hjust=1, size=9),
         axis.text.y=element_text(size=12),
         axis.ticks.x=element_blank())                         

figfn <- paste(outdir, "Figure1.2_corr_MCls.box.png", sep="")
ggsave(figfn, plot=p1, width=480, height=250, units="px", dpi=120) ## ggsave 




###
### pair treatments within the same cell types, each box has 8 points
fn <- paste(outdir, "1.x_DEG_corr.rds", sep="")
mat <- read_rds(fn)
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
 
p3 <- ggplot(DFcorr3, aes(x=comb, y=as.numeric(rr)))+ ##, fill=factor(MCls), group=comb))+
   ## geom_violin()+
   geom_boxplot(outlier.shape=NA, lwd=0.25)+
   stat_summary(fun=median, color="grey", geom="point", shape=23, size=1.8, stroke=0.9)+
   geom_jitter(width=0.2, size=0.5, aes(color=MCls))+ 
   scale_color_manual(values=col1)+
   ylab("SCC (DEGs)")+ylim(0,0.8)+
   theme_bw()+
   theme(legend.position="none",
         axis.title.x=element_blank(),
         axis.title.y=element_text(size=12),
         axis.text.x=element_text(angle=45, hjust=1, size=7.5),
         axis.text.y=element_text(size=12))                         
  
figfn <- paste(outdir, "Figure1.3_pairtreats.box.png", sep="")
ggsave(figfn, plot=p3, width=480, height=350, units="px", dpi=120)



###
### boxplots group by treats
fn <- paste(outdir, "1.x_DEG_corr.rds", sep="")
corr_mat <- read_rds(fn)

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


treats <- unique(cvt2$treats)
DFcorr2 <- map_dfr(treats, function(ii){
   ###
   rnSel <- cvt2%>%filter(treats==ii)%>%pull(rn)
   ntr <- length(rnSel)
   ### 
   DF2 <- NULL     
   for ( i in 1:(ntr-1)){
      for ( j in (i+1):ntr){
          ##
          rn_i <- rnSel[i]
          rn_j <- rnSel[j]
          rr0 <- round(corr_mat[rn_i, rn_j], 3)
          ###
          df2 <- data.frame(comb=rn_i, comb2=rn_j, rr=rr0)
          DF2 <- rbind(DF2, df2)
      }    
   }    
   DF2$treats <- ii
   DF2 
})


DFcorr2 <- DFcorr2%>%mutate(rr=round(rr,3))%>%filter(rr<1)
p4 <- ggplot(DFcorr2, aes(x=treats, y=rr, fill=treats))+
   ##geom_violin()+ 
   geom_boxplot(outlier.shape=NA, lwd=0.25)+
   stat_summary(fun=median, color="grey", geom="point", shape=23, size=2.5, stroke=0.9)+
   geom_jitter(width=0.2, size=0.2)+ 
   scale_fill_manual(values=col2)+
   ylab("SCC (DEGs)")+ylim(0,0.8)+
   theme_bw()+
   theme(legend.position="none",
         axis.title.x=element_blank(),
         axis.title.y=element_text(size=12),
         axis.text.x=element_blank(),  ##element_text(angle=45, hjust=1, size=9),
         axis.text.y=element_text(size=12),
         axis.ticks.x=element_blank())                         
###
figfn <- paste(outdir, "Figure1.4_corr_treat.box.png", sep="")
ggsave(figfn, plot=p4, width=480, height=250, units="px", dpi=120) ## ggsave 





#################################################
### summary Differentially accessible regions ###
#################################################

rm(list=ls())

outdir <- "./Plots_pub/2_diff_plots/"



###
### setting colors
col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
  "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
   "6_Monocyte"="#984ea3", "7_dnT"="black")
col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
       "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3")



###
### DE results

fn <- "./sc_multiome_data/2_Differential/1.2_DiffPeak.outs/2.0_DESeq.results_0.1_cn_sampleID+new_treat.rds"
res <- read_rds(fn)%>%as.data.frame()%>%mutate(comb=paste(MCls, contrast, sep="_"))

res2 <- res%>%filter(p.adjusted<0.1, abs(estimate)>0.5)%>%
   mutate(direction=ifelse(estimate>0, "Up", "Down")) # sign(estimate)))


### report numbers of DARs used for paper 
res2%>%group_by(contrast)%>%summarize(ngene=length(unique(gene)), .groups="drop")%>%arrange(desc(ngene))
res2%>%group_by(MCls)%>%summarize(ngene=length(unique(gene)), .groups="drop")%>%arrange(desc(ngene))
sigs <- unique(res2$gene)
x <- res2%>%filter(comb=="6_Monocyte_vitD")
table(x$direction)


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
    row_names_gp=gpar(fontsize=12), row_names_side="left",      
    column_names_gp=gpar(fontsize=12), column_names_rot=45,
    heatmap_legend_param=list(at=round(seq(olap_min, olap_max, length.out=5),0),
       grid_width=grid::unit(0.4, "cm"), legend_height=grid::unit(6, "cm"), title="",
       title_gp=gpar(fontsize=12), labels_gp=gpar(fontsize=12), direction="horizontal"),
    cell_fun=function(j, i, x, y, width, height, fill){
        grid.text(mat2[i,j], x, y, gp=gpar(fontsize=12))
    })    


figfn <- paste(outdir, "Figure2.1_DARs.heatmap.png", sep="")
png(figfn, width=500, height=520, res=120)
draw(p, heatmap_legend_side="bottom")
dev.off()


## ggsave(figfn, plot=p, device="png", width=580, height=620, units="px", dpi=300) ## ggsave 
## unlink(figfn)


###
### bar plots for top annotation 
summDF <- res2%>%group_by(contrast)%>%summarize(ny=length(unique(gene)),.groups="drop")

## col=list(celltype=c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56", "3_TEM"="blue",
##                    "4_Bcell"="#4daf4a", "5_CD8Naive"="green", "6_Monocyte"="#984ea3", "7_dnT"="black"),
##             treats=)
 
p2 <- ggplot(summDF, aes(x=contrast, y=ny, fill=contrast))+ ##, color=contrast, alpha=direction))+
    geom_bar(stat="identity")+
    scale_fill_manual(values=col2, guide="none")+
    ## scale_color_manual(values=c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
    ##    "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3"), guide="none")+
    ## scale_alpha_manual(values=c("Down"=0.5, "Up"=1),
    ##      guide=guide_legend(override.aes=list(fill="red", size=1.5)))+
    ylab("DARs")+                                                     
    theme_bw()+
    theme(legend.position="none", ##legend.position=c(0.8, 0.8),
          legend.title=element_blank(),
          legend.key=element_blank(),
          legend.box.background=element_blank(),
          legend.background=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_text(size=12),
          axis.text.x=element_blank(),
          axis.text.y=element_text(size=12),
          axis.ticks.x=element_blank())
##
##
figfn <- paste(outdir, "Figure2.1_DARs_treat.bar.png", sep="")
ggsave(figfn, plot=p2, width=480, height=220, units="px", dpi=120) ## ggsave 
##unlink(figfn)

## figfn <- paste(outdir, "Figure2.2_bar_DEG.png", sep="")
## png(figfn, width=550, height=380, res=120)
## p2
## dev.off()


###
### bar plots for row annotation, group_by, cell type

## summDF2 <- res2%>%group_by(MCls, direction)%>%summarize(ny=n(),.groups="drop")%>%
##     group_by(MCls)%>%summarize(ny=sum(ny),.groups="drop")

summDF2 <- res2%>%group_by(MCls)%>%summarize(ny=n(),.groups="drop")

p3 <- ggplot(summDF2, aes(x=MCls, y=ny, fill=MCls))+ ##, color=MCls, alpha=direction))+
    geom_bar(stat="identity")+
    coord_flip()+
    ##scale_y_reverse()+
    scale_fill_manual(values=col1, guide="none")+
    ## scale_color_manual(values=c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
    ##     "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
    ##     "6_Monocyte"="#984ea3", "7_dnT"="black"), guide="none")+
    ## scale_alpha_manual(values=c("Down"=0.5, "Up"=1),
    ##      guide=guide_legend(override.aes=list(fill="#ffaa00", size=1.5)))+
    ylab("DARs")+                                                     
    theme_bw()+
    theme(legend.position="none",
          legend.title=element_blank(),
          legend.key=element_blank(),
          legend.box.background=element_blank(),
          legend.background=element_blank(),
          axis.title.x=element_text(size=12),
          axis.title.y=element_blank(),
          axis.text.x=element_text(size=12),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())

figfn <- paste(outdir, "Figure2.1_DARs_MCls.bar.png", sep="")
ggsave(figfn, plot=p3, width=250, height=480, units="px", dpi=120) ## ggsave 

## figfn <- paste(outdir, "Figure1.3_bar_DEG.png", sep="")
## png(figfn, width=320, height=450, res=120)
## p3
## dev.off()


 

####################################################
#### Use union of peak to calculate correlation ####
####################################################
 
rm(list=ls())

outdir <- "./Plots_pub/2_diff_plots/"


fn <- "./sc_multiome_data/2_Differential/1.2_DiffPeak.outs/2.0_DESeq.results_0.1_cn_sampleID+new_treat.rds"
res <- read_rds(fn)%>%as.data.frame()%>%mutate(comb=paste(MCls, contrast, sep="_"))

sigs <- res%>%filter(p.adjusted<0.1, abs(estimate)>0.5)%>%pull(gene)%>%unique() ## 17774

res2 <- res%>%filter(gene%in%sigs)%>%dplyr::select(comb, gene, estimate, statistic)
 
###
### data format #DEG*combination
### Orginal data is a matrix of 17,774*48 and after filtering missing value, it is 17,445*48
mat <- res2%>%pivot_wider(id_cols=gene, names_from=comb, values_from=estimate)%>%
    column_to_rownames(var="gene")%>%as.matrix()
     
rnz <- rowSums(is.na(mat))
mat2 <- mat[rnz==0,] ### 17,445*48

### correlation data for plot
corr_mat <- cor(mat2, method="spearman")

opfn <- paste(outdir, "2.x_DP_corr.rds", sep="")
write_rds(corr_mat, file=opfn)

###
### peaks
## MCls <- unique(res$MCls)
## DFpeak <- map_dfr(MCls, function(ii){
##    ###
##    peakSel <- res%>%filter(p.adjusted<0.1, abs(estimate)>0.5, MCls==ii)%>%pull(gene)%>%unique() 
##    DF <- data.frame(peak=peakSel)
##    DF$MCls <- ii
##    DF
## })    

## ## summ <- DFpeak%>%group_by(MCls)%>%summarize(npeak=n(), .groups="drop")
## ## write.xlsx(summ, "./2_diff_plots.outs/2_summ.xlsx")


## ###
## ### correlation 
## comb <- sort(unique(res$comb))
## ncomb <- length(comb)
## cvt <- str_split(comb, "_", simplify=T)
## DFcomb <- data.frame(rn=comb, MCls=paste(cvt[,1], cvt[,2], sep="_"), treats=cvt[,3])


## corr_mat <- matrix(0, ncomb, ncomb)
## colnames(corr_mat) <- comb
## rownames(corr_mat) <- comb
## MCls <- sort(unique(DFcomb$MCls))
## for (oneMCl in MCls){
##    ###
##    peakSel <- DFpeak%>%filter(MCls==oneMCl)%>%pull(peak)%>%unique()
##    rn2 <- DFcomb%>%filter(MCls==oneMCl)%>%pull(rn)%>%unique()%>%sort()

##    res2 <- res%>%filter(comb%in%rn2, gene%in%peakSel)%>%select(comb, gene, LFC=estimate)
##    dfmat <- res2%>%pivot_wider(id_cols=gene, names_from=comb, values_from=LFC)
##    rr_mat  <- cor(dfmat[,-1], method="spearman")
##    corr_mat[rn2, rn2] <- rr_mat[rn2, rn2]     
##    cat(oneMCl, "\n")
## }    
     
## opfn <- paste(outdir, "2.x_DP_corr.rds", sep="")
## write_rds(corr_mat, file=opfn)




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
figfn <- paste(outdir, "FigS2_8_LFC.peak_corr.heatmap.png", sep="")
## ggsave(figfn, plot=p2, device="png", width=680, height=700, units="px", dpi=300) ## ggsave 
png(figfn, width=580, height=700,res=120)
p2 <- draw(p2)
dev.off()




###################################################################
#### boxplots show spearman correlation of LFC of union of DARs ###
###################################################################

rm(list=ls())

outdir <- "./Plots_pub/2_diff_plots/"



###
### setting colors
col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
  "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
   "6_Monocyte"="#984ea3", "7_dnT"="black")
col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
       "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3")



###
### plot data group by cell types

fn <- paste(outdir, "2.x_DP_corr.rds", sep="")
corr_mat <- read_rds(fn)

rn<- rownames(corr_mat)
cvt <- str_split(rownames(corr_mat), "_", simplify=T)
cvt2 <- data.frame(rn=rownames(corr_mat), MCls=paste(cvt[,1], cvt[,2], sep="_"), treats=cvt[,3])

###
MCls <- unique(cvt2$MCls)
DFcorr <- map_dfr(MCls, function(ii){
   ###
   rnSel <- cvt2%>%filter(MCls==ii)%>%pull(rn)
   ntr <- length(rnSel)

   ### 
   DF2 <- NULL     
   for ( i in 1:(ntr-1)){
      for ( j in (i+1):ntr){
          ##
          rn_i <- rnSel[i]
          rn_j <- rnSel[j]
          rr0 <- round(corr_mat[rn_i, rn_j], 3)
          ###
          df2 <- data.frame(comb=rn_i, comb2=rn_j, rr=rr0)
          DF2 <- rbind(DF2, df2)
      }    
   }

   ### 
   DF2$MCls <- ii
   DF2 
})
 

 

###
### boxplot group by cell types
DFcorr <- DFcorr%>%mutate(rr=round(rr,3))%>%filter(rr<1)
p1 <- ggplot(DFcorr, aes(x=MCls, y=rr, fill=MCls))+
   ##geom_violin()+ 
   geom_boxplot(outlier.shape=NA, lwd=0.25)+
   stat_summary(fun=median, color="grey", geom="point", shape=23, size=2.5, stroke=0.9)+
   geom_jitter(width=0.2, size=0.2)+     
   scale_fill_manual(values=col1)+
   ylab("SCC (DARs)")+ylim(-0.1,0.8)+
   theme_bw()+
   theme(legend.position="none",
         axis.title.x=element_blank(),
         axis.title.y=element_text(size=12),
         axis.text.x=element_blank(),  ##element_text(angle=45, hjust=1, size=9),
         axis.text.y=element_text(size=12),
         axis.ticks.x=element_blank())                         

###
figfn <- paste(outdir, "Figure2.2_corr_MCls.box.png", sep="")
ggsave(figfn, plot=p1, width=480, height=250, units="px", dpi=120) ## ggsave 




###
### pair treatments within the same cell types, each box has 8 points
fn <- paste(outdir, "2.x_DP_corr.rds", sep="")
mat <- read_rds(fn)
##
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
 
p3 <- ggplot(DFcorr3, aes(x=comb, y=as.numeric(rr)))+ ##, fill=factor(MCls), group=comb))+
   ## geom_violin()+
   geom_boxplot(outlier.shape=NA, lwd=0.25)+
   stat_summary(fun=median, color="grey", geom="point", shape=23, size=1.8, stroke=0.9)+
   geom_jitter(width=0.2, size=0.5, aes(color=MCls))+ 
   scale_color_manual(values=col1)+
   ylab("SCC (DARs)")+ylim(-0.1,0.8)+
   theme_bw()+
   theme(legend.position="none",
         axis.title.x=element_blank(),
         axis.title.y=element_text(size=12),
         axis.text.x=element_text(angle=45, hjust=1, size=7.5),
         axis.text.y=element_text(size=12))                         
  
figfn <- paste(outdir, "Figure2.3_pairtreats.box.png", sep="")
ggsave(figfn, p3,  width=480, height=350, units="px", dpi=120)



###
### boxplot, spearman's correlation  between cell types for each treat

fn <- paste(outdir, "2.x_DP_corr.rds", sep="")
corr_mat <- read_rds(fn)

rn<- rownames(corr_mat)
cvt <- str_split(rownames(corr_mat), "_", simplify=T)
cvt2 <- data.frame(rn=rownames(corr_mat), MCls=paste(cvt[,1], cvt[,2], sep="_"), treats=cvt[,3])


###
treats <- unique(cvt2$treats)
DFcorr2 <- map_dfr(treats, function(ii){
   ###    
   rnSel <- cvt2%>%filter(treats==ii)%>%pull(rn)
   ntr <- length(rnSel)
    
   ### 
   DF2 <- NULL     
   for ( i in 1:(ntr-1)){
      for ( j in (i+1):ntr){
          ##
          rn_i <- rnSel[i]
          rn_j <- rnSel[j]
          rr0 <- round(corr_mat[rn_i, rn_j], 3)
          ###
          df2 <- data.frame(comb=rn_i, comb2=rn_j, rr=rr0)
          DF2 <- rbind(DF2, df2)
      }    
   }
   DF2$treats <- ii
   DF2 
})



## DFcorr2 <- DFcorr2%>%mutate(rr=round(rr,3))%>%filter(rr<1)
p4 <- ggplot(DFcorr2, aes(x=treats, y=rr, fill=treats))+
   ##geom_violin()+ 
   geom_boxplot(outlier.shape=NA, lwd=0.25)+
   stat_summary(fun=median, color="grey", geom="point", shape=23, size=2.5, stroke=0.9)+
   geom_jitter(width=0.2, size=0.2)+ 
   scale_fill_manual(values=col2)+
   ylab("SCC (DARs)")+ylim(0,0.8)+
   theme_bw()+
   theme(legend.position="none",
         axis.title.x=element_blank(),
         axis.title.y=element_text(size=12),
         axis.text.x=element_blank(),  ##element_text(angle=45, hjust=1, size=9),
         axis.text.y=element_text(size=12),
         axis.ticks.x=element_blank())                         
###
figfn <- paste(outdir, "Figure2.4_corr_treat.box.png", sep="")
ggsave(figfn, plot=p4, width=480, height=250, units="px", dpi=120) ## ggsave 





########################################################################### 
### compare treat-treat pair of correlation between LFC on RNA and ATAC ###
###########################################################################


rm(list=ls())

outdir <- "./Plots_pub/2_diff_plots/"


fn <- paste(outdir, "1.x_DEG_corr.rds", sep="") 
mat <- read_rds(fn)
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


fn2 <- paste(outdir, "2.x_DP_corr.rds", sep="")
mat2 <- read_rds(fn2)

        
plotDF <- map_dfr(comb, function(ii){
   ##
   tr1 <- gsub("_.*", "", ii)
   tr2 <- gsub(".*_", "", ii) 
   df <- map_dfr(MCls, function(oneMCl){
      ##
      comb1 <- paste(oneMCl, tr1, sep="_")
      comb2 <- paste(oneMCl, tr2, sep="_")
      df2 <- data.frame(comb=ii, MCls=oneMCl, rr_RNA=mat[comb1, comb2], rr_ATAC=mat2[comb1, comb2])
      df2
   })
   df
})    


col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
    "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
     "6_Monocyte"="#984ea3", "7_dnT"="black")

corr <- cor.test(plotDF$rr_RNA, plotDF$rr_ATAC, method="spearman")
rr <- round(as.numeric(corr$estimate), digits=3)
eq <- deparse(bquote(italic(rho)==.(rr)~" ***"))

  



###
### add identical line and fitting lines
  
p <- ggplot(plotDF, aes(x=rr_ATAC, y=rr_RNA, color=MCls))+
   geom_point(size=1.5)+
   annotate("text", label=eq, x=0.3, y=0.75, parse=T)+ 
   geom_smooth(method="lm", se=F, linewidth=0.6)+
   geom_abline(color="grey")+ 
   scale_color_manual(values=col1)+      
   xlab("pairs' SCC (DARs)")+xlim(0.1, 0.8)+ylab("pairs' SCC (DEGs)")+ylim(0.1, 0.8)+
   theme_bw()+
   theme(legend.position="none",
         axis.title=element_text(size=12),
         axis.text=element_text(size=12))
##
figfn <- paste(outdir, "Figure2.5_compare.scatter.png", sep="")
ggsave(figfn, p, width=320, height=270, units="px", dpi=120)
          




######################
### summary Tables ###
######################


rm(list=ls())

outdir <- "./Plots_pub/2_diff_plots/"
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)


###
### summary DEGs 
fn <- "./sc_multiome_data/2_Differential/1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn_sampleID+new_treat.rds"
res <- read_rds(fn)%>%as.data.frame()%>%
    mutate(comb=paste(MCls, contrast, sep="_"), is_sig=ifelse(p.adjusted<0.1&abs(estimate)>0.5, 1, 0))%>%
    drop_na(estimate)

summ <- res%>%group_by(comb)%>%summarize(test_genes=n(), DEGs=sum(is_sig, na.rm=T), .groups="drop")



###
### summary DARs
fn <- "./sc_multiome_data/2_Differential/1.2_DiffPeak.outs/2.0_DESeq.results_0.1_cn_sampleID+new_treat.rds"
res2 <- read_rds(fn)%>%as.data.frame()%>%
    mutate(comb=paste(MCls, contrast, sep="_"), is_sig=ifelse(p.adjusted<0.1&abs(estimate)>0.5, 1, 0))%>%
    drop_na(estimate)

summ2 <- res2%>%group_by(comb)%>%summarize(test_peaks=n(), DARs=sum(is_sig, na.rm=T), .groups="drop")

df_comb <- summ%>%full_join(summ2, by="comb")

opfn <- paste(outdir, "TableS2_1_summary_differential.xlsx", sep="")
write.xlsx(df_comb, file=opfn)


###
### output summary data

res <- res%>%
    dplyr::select(comb, MCls, contrast, gene, baseMean, estimate, stderror, statistic, p.value, p.adjusted)
##
opfn <- paste(outdir, "TableS2_2_differential_gene.txt.gz", sep="")
fwrite(res, file=opfn, quote=F, sep=" ", na=NA)


###
res2 <- res2%>%
    dplyr::select(comb, MCls, contrast, peak=gene, baseMean, estimate, stderror, statistic, p.value, p.adjusted)
##
opfn <- paste(outdir, "TableS2_3_differential_chromatin.txt.gz", sep="")
fwrite(res2, file=opfn, quote=F, sep=" ", na=NA)
