##
library(tidyverse)
library(qvalue)
library(ComplexHeatmap)
library(cowplot)
library(RColorBrewer)
library(circlize)
library(viridis)
library(openxlsx)
library(glue)
library(ggtext)



#######################################################
### examples of cell-type specific genes
########################################################

rm(list=ls())

outdir2 <- "./3_RNA_summary.outs/Example_genes/"
if ( !file.exists(outdir2) ) dir.create(outdir2, showWarnings=F, recursive=T)


### 1-lrt results
fn <- "./1_inter_RNA.outs/1.2_lrt_results.rds"
res_lrt <- read_rds(fn)
res_lrt2 <- res_lrt%>%filter(p.adjusted<0.1)%>%arrange(p.value)
sigs <- res_lrt2$gene

## gene1 <- unique(res2$gene)
## opfn <- paste(outdir2, "inter1_lrt_genes.txt", sep="") 
## write.table(unique(res2$gene), file=opfn, row.names=FALSE, col.names=FALSE, quote=FALSE)


### 2-compared effects to 0_CD4Naive 
fn <- "./3_RNA_summary.outs/1_MCls.diff.rds"
res_diff <- read_rds(fn)%>%mutate(treats=gsub(".*_", "", conditionX))

res_diff2 <- res_diff%>%filter(p.value<0.01, gene%in%sigs)%>%
    group_by(gene)%>%slice_max(order_by=abs(statistic), n=1, with_ties=FALSE)%>%ungroup()%>%
    arrange(p.value)%>%as.data.frame()



### 3-results in contrast to control, for plot data 
fn <- "./1_inter_RNA.outs/2.2_each_pair.results.rds"
DF <- read_rds(fn)



### colors           
col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
  "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
   "6_Monocyte"="#984ea3", "7_dnT"="black")


###
### Examples plots

for ( i in 1:10){
    
gene0 <- res_diff2$gene[i]
treat2 <- res_diff2$treats[i]
cat(gene0, treat2, "\n")

    
### obtain padj from interaction II 
df2 <- res_diff%>%filter(gene==gene0, treats==treat2)%>%
    dplyr::select(condition=conditionX, p_diff=p.value, padj_diff=padj2)

    
### plot data
plotDF <- DF%>%filter(gene==gene0, grepl(treat2, condition))%>%left_join(df2, by="condition")
x <- str_split(plotDF$condition, "_", simplify=TRUE)
plotDF2 <- plotDF%>%
    mutate(MCls=paste(x[,1],x[,2], sep="_"), treats=x[,3], MCls_val=as.numeric(x[,1]))%>%
    mutate(bhat=estimate, b_upper=bhat+1.96*stderror, b_lower=bhat-1.96*stderror,
           MCl2=ifelse(p_diff<0.01&!is.na(p_diff), glue("<b>{MCls}"), glue("<i>{MCls}")),
           MCl2=fct_reorder(MCl2, MCls_val))


###
### plots
p <- ggplot(plotDF2, aes(x=bhat, y=MCl2))+
   geom_errorbarh(aes(xmax=b_upper, xmin=b_lower, colour=MCls),
       linewidth=0.5, height=0.2)+
   geom_point(aes(colour=MCls), shape=19, size=1.5)+
   scale_colour_manual(values=col1)+
   geom_vline(aes(xintercept=0), size=0.25, linetype="dashed")+
   xlab("Effects of treatment")+ 
   ## xlab("log odds ratio")+xlim(-3, 6)+
   ggtitle(bquote(italic(.(gene0))~"in"~.(treat2)))+ 
   theme_bw()+
   theme(plot.title=element_text(hjust=0.5, size=12),
         axis.title.x=element_text(size=10),
         axis.title.y=element_blank(),
         axis.text.x=element_text(size=10),
         axis.text.y=element_markdown(size=10),         
         legend.position="none")
         ## legend.title=element_blank(),
         ## legend.text=element_text(size=8),
         ## legend.key.size=unit(0.4, "cm"))
         ##legend.position="none")
 
figfn <- paste(outdir2, "Figure", i, "_", gene0, "_", treat2, ".forest.png", sep="")
ggsave(figfn, p,  width=450, height=550, units="px", dpi=120)

}

 
    

######################################
### Examples show tables
#####################################


for ( i in 1:10){
    
gene0 <- res_diff2$gene[i]
treat2 <- res_diff2$treats[i]
cat(gene0, treat2, "\n")

    
### obtain padj from interaction II 
tmp <- res_diff%>%filter(gene==gene0, treats==treat2)%>%as.data.frame()
 


### plot data
DF2 <- DF%>%filter(gene==gene0, grepl(treat2, condition))%>%as.data.frame()
b <- DF2$estimate
names(b) <- DF2$condition
se <- DF2$stderror
names(se) <- DF2$condition
 
tmp2 <- tmp%>%
    mutate(baseMean=round(baseMean, 3), estimate=round(estimate, 3), stderror=round(stderror, 3),
           statistic=round(statistic, 3),
           bx=round(b[conditionX],3), sx=round(se[conditionX], 3),
           by=round(b[conditionY],3), sy=round(se[conditionY],3))
tmp2 <- tmp2%>%dplyr::select(-p.adjusted, -conditionY, -padj2, -treats)
     

opfn <- paste(outdir2, "Table", i, "_", gene0, "_", treat2, ".xlsx", sep="")
write.xlsx(tmp2, file=opfn)

}    

















##############################################################
### Heatmap of different effects size across cell-types 
##############################################################

## rm(list=ls())
## outdir <- "./3_RNA_summary.outs/"
 
 
## fn <- "./3_RNA_summary.outs/1_MCls.diff.rds"
## res <- read_rds(fn)
## sigs <- res%>%filter(p.adjusted<0.1)%>%pull(gene)%>%unique()

## dfmat <- res%>%
##     pivot_wider(id_cols=gene, names_from=conditionX, values_from=statistic, values_fill=NA)%>%
##     column_to_rownames(var="gene")%>%as.matrix()

## dfmat2 <- dfmat[sigs,]



## ### color for heatmap value
## y <- as.vector(dfmat2)
## ##y0 <- y[abs(y)<2] ### quantile(abs(y), probs=0.99)
## mybreak <- c(min(y), seq(-6, 6, length.out=98), max(y)) 
## mycol <- colorRamp2(mybreak, colorRampPalette(rev(brewer.pal(n=7,name="RdBu")))(100))

## quantile(abs(y), probs=c(0.9, 0.95, 0.99))


## ###
## ### annotation columns
## ## col0 <- c("1"="#8c510a", "2"="#d8b365")
## col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
##   "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
##    "6_Monocyte"="#984ea3", "7_dnT"="black")
## col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
##        "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3")


## x <- str_split(colnames(dfmat2), "_", simplify=T)
## x2 <- data.frame(comb=colnames(dfmat2),
##    MCls=paste(x[,1], x[,2], sep="_"), treat=x[,3])
## x2 <- x2%>%arrange(treat)
## df_col <- x2%>%dplyr::select(MCls, treat) ##%>%mutate(Feature=as.character(Feature))

## col_ha <- HeatmapAnnotation(df=df_col, col=list(MCls=col1, treat=col2),
##     annotation_name_gp=gpar(fontsize=10),
##     annotation_legend_param=list(       
##   MCls=list(title_gp=gpar(fontsize=9), labels_gp=gpar(fontsize=9), title="celltype",
##                 grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm")),
##   treat=list(title_gp=gpar(fontsize=9), labels_gp=gpar(fontsize=9), title="treatment",
##                 grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm"))))



## ## show_legend=c(F,F))
## ### main plots
## p2 <- Heatmap(dfmat2, col=mycol,
##    cluster_rows=T, cluster_columns=F, row_km=7,
##    show_row_names=F, row_names_gp=gpar(fontsize=4),
##    show_column_names=T, column_names_gp=gpar(fontsize=6),
##    show_row_dend=T, row_dend_width=grid::unit(2, "cm"),
##    show_column_dend=F,
##    column_names_rot=-45,   
##    top_annotation=col_ha,
##    heatmap_legend_param=list(title=bquote(~italic(Z)~"-score"),
##       title_gp=gpar(fontsize=9),
##       at=seq(-6, 6, by=3), 
##       labels_gp=gpar(fontsize=9),
##       grid_width=grid::unit(0.45, "cm"),
##       legend_height=grid::unit(7.8, "cm")))
 
## ###
## figfn <- paste(outdir, "Figure1_zscore.heatmap.png", sep="")
## ## ggsave(figfn, plot=p2, device="png", width=680, height=700, units="px", dpi=300) ## ggsave 
## png(figfn, width=1000, height=1200,res=120)
## set.seed(0)
## p2 <- draw(p2)
## dev.off()





## ## nna <- rowSums(is.na(dfmat))



## #################
## ### lrt
## ##################

## ###
## fn <- "./3_RNA_summary.outs/1_MCls.diff.rds"
## res <- read_rds(fn)
## sigs <- res%>%filter(p.adjusted<0.05)%>%pull(gene)%>%unique()


## fn <- "./1_inter_RNA.outs/1.2_lrt_results.rds"
## res2 <- read_rds(fn)
## sig2 <- res2%>%filter(p.adjusted<0.05)%>%pull(gene)%>%unique()

## shared <- intersect(sigs, sig2)

## res2 <- res2%>%filter(gene%in%shared)%>%arrange(p.value)



## ###
## ### forest plots

## outdir2 <- paste(outdir, "Example_genes", sep="")
