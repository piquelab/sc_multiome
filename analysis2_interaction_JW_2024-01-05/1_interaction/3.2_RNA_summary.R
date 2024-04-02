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
library(annotables)



outdir <- "./3_RNA_summary.outs/"
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)


##################################################
### Heatmap of number of interaction genes 
##################################################

### first passing LRT 10% FD and then p<0.01 for each contrastXR

rm(list=ls())
outdir <- "./3_RNA_summary.outs/"

### lrt analysis
fn <- "./1_inter_RNA.outs/1.2_lrt_results.rds"
res_lrt <- read_rds(fn)

 
gene0 <- unique(res_lrt$gene)
autosome <- as.character(1:22)
gene1 <- grch38%>%filter(chr%in%autosome)%>%pull(symbol)%>%unique()
gene2 <- intersect(gene0, gene1)


sigs_lrt <- res_lrt%>%filter(p.adjusted<0.05, gene%in%gene2)%>%pull(gene)%>%unique()


### compare contrast
fn <- "./3_RNA_summary.outs/1_MCls.diff.rds"
res <- read_rds(fn)

sigs <-res%>%filter(gene%in%sigs_lrt, p.value<0.01)%>%pull(gene)%>%unique()

    
## summary #genes with interaction effects
summ <- res%>%filter(gene%in%sigs_lrt)%>%group_by(conditionX)%>%
    summarise(nsig=sum(p.value<0.01, na.rm=T), .groups="drop")

x <- str_split(summ$conditionX, "_", simplify=T)
summ <- summ%>%mutate(MCls=paste(x[,1], x[,2], sep="_"), treats=x[,3])
summ2 <- summ%>%arrange(treats)

opfn <- paste(outdir, "2_inter_gene_lrt0.1_p0.01.xlsx", sep="")
write.xlsx(summ2, file=opfn)



###
### Heatmap
fn <- "./3_RNA_summary.outs/2_inter_gene_lrt0.1_p0.01.xlsx"
x <- read.xlsx(fn)

mat2 <- x%>%pivot_wider(id_cols=MCls, names_from=treats, values_from=nsig)%>%
    column_to_rownames(var="MCls")%>%as.matrix()

MCl_sort <- sort(rownames(mat2), decreasing=TRUE)
mat2 <- mat2[MCl_sort,]

### setting colors
olap_min <- min(mat2)
olap_max <- max(mat2)
mycol <- colorRamp2(seq(olap_min, olap_max, length.out=20), colorRampPalette(brewer.pal(n=7,name="YlGnBu"))(20))

 
p <- Heatmap(mat2, name="#genes", col=mycol, cluster_rows=F, cluster_columns=F,             
    row_names_gp=gpar(fontsize=12), row_names_side="left",      
    column_names_gp=gpar(fontsize=12), column_names_rot=45,
    heatmap_legend_param=list(at=round(seq(olap_min, olap_max, length.out=5),0),
       grid_width=grid::unit(0.4, "cm"), legend_height=grid::unit(6, "cm"), title="",
       title_gp=gpar(fontsize=9), labels_gp=gpar(fontsize=12)), ### direction="horizontal"),
    cell_fun=function(j, i, x, y, width, height, fill){
        grid.text(mat2[i,j], x, y, gp=gpar(fontsize=12))
    })    
 
figfn <- paste(outdir, "Figure1.1_genes_lrt0.1_p0.01.heatmap.png", sep="")
png(figfn, width=550, height=500, res=120)
draw(p) ## heatmap_legend_side="bottom")
dev.off()



###
### Barplots
summDF<- res%>%filter(gene%in%sigs_lrt, p.value<0.01)%>%mutate(treats=gsub(".*_", "", conditionX))%>%
    group_by(treats)%>%summarize(ny=length(unique(gene)), .groups="drop")

col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
       "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3")

p2 <- ggplot(summDF, aes(x=treats, y=ny, fill=treats))+ ## , color=contrast))+ ##, alpha=direction))+
    geom_bar(stat="identity")+
    scale_fill_manual(values=col2, guide="none")+
    geom_text(aes(x=treats, y=ny+50, label=ny), size=3)+
    ylab("#genes")+                                                     
    theme_bw()+
    theme(legend.position="none", ##c(0.8, 0.8),
          legend.title=element_blank(),
          legend.key=element_blank(),
          legend.box.background=element_blank(),
          legend.background=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_text(size=10),
          axis.text.x=element_text(size=10, angle=45, hjust=1),
          axis.text.y=element_text(size=10))

figfn <- paste(outdir, "Figure1.2_treat.bar.png", sep="")
ggsave(figfn, plot=p2, width=380, height=320, units="px", dpi=120) ## ggsave                                 



######################################################################
### Heatmap of #intersection genes overlapping with lrt 
######################################################################


## rm(list=ls())
## outdir <- "./3_RNA_summary.outs/"

## ### lrt results 
## fn <- "./1_inter_RNA.outs/1.2_lrt_results.rds"
## res_lrt <- read_rds(fn)
## sigs_lrt <- res_lrt%>%filter(p.adjusted<0.1)%>%pull(gene)%>%unique()
 
## ###
## fn <- "./3_RNA_summary.outs/1_MCls.diff.rds"
## res <- read_rds(fn)

## sigs <- res%>%filter(p.adjusted<0.1)%>%pull(gene)%>%unique()

## olap <- intersect(sigs_lrt, sigs)


## ##
## summ <- res%>%group_by(conditionX)%>%
##     summarise(nsig_0.1=sum(p.adjusted<0.1&gene%in%sigs_lrt, na.rm=T),
##               nsig_0.05=sum(p.adjusted<0.05&gene%in%sigs_lrt, na.rm=T),  .groups="drop")
## x <- str_split(summ$conditionX, "_", simplify=T)
## summ <- summ%>%mutate(MCls=paste(x[,1], x[,2], sep="_"), treats=x[,3])
## summ2 <- summ%>%arrange(treats)

## opfn <- paste(outdir, "2.1_inter_gene_olap.xlsx", sep="")
## write.xlsx(summ2, file=opfn)



## fn <- "./3_RNA_summary.outs/2.1_inter_gene_olap.xlsx"
## x <- read.xlsx(fn)

## mat2 <- x%>%pivot_wider(id_cols=MCls, names_from=treats, values_from=nsig_0.1)%>%
##     column_to_rownames(var="MCls")%>%as.matrix() 

## ### setting colors
## olap_min <- min(mat2)
## olap_max <- max(mat2)
## mycol <- colorRamp2(seq(olap_min, olap_max, length.out=20), colorRampPalette(brewer.pal(n=7,name="YlGnBu"))(20))

 
## p <- Heatmap(mat2, name="#genes", col=mycol, cluster_rows=F, cluster_columns=F,             
##     row_names_gp=gpar(fontsize=12), row_names_side="left",      
##     column_names_gp=gpar(fontsize=12), column_names_rot=45,
##     heatmap_legend_param=list(at=round(seq(olap_min, olap_max, length.out=5),0),
##        grid_width=grid::unit(0.4, "cm"), legend_height=grid::unit(6, "cm"), title="",
##        title_gp=gpar(fontsize=9), labels_gp=gpar(fontsize=12)), ### direction="horizontal"),
##     cell_fun=function(j, i, x, y, width, height, fill){
##         grid.text(mat2[i,j], x, y, gp=gpar(fontsize=12))
##     })    

## figfn <- paste(outdir, "Figure2.1_genes_q0.1_olap.heatmap.png", sep="")
## png(figfn, width=550, height=500, res=120)
## draw(p) ## heatmap_legend_side="bottom")
## dev.off()


## ###
## ### Barplots
## summDF<- res%>%filter(p.adjusted<0.1, gene%in%sigs_lrt)%>%mutate(treats=gsub(".*_", "", conditionX))%>%
##     group_by(treats)%>%summarize(ny=length(unique(gene)), .groups="drop")

## col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
##        "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3")

## p2 <- ggplot(summDF, aes(x=treats, y=ny, fill=treats))+ ## , color=contrast))+ ##, alpha=direction))+
##     geom_bar(stat="identity")+
##     scale_fill_manual(values=col2, guide="none")+
##     geom_text(aes(x=treats, y=ny+50, label=ny), size=3)+
##     ylab("#genes")+                                                     
##     theme_bw()+
##     theme(legend.position="none", ##c(0.8, 0.8),
##           legend.title=element_blank(),
##           legend.key=element_blank(),
##           legend.box.background=element_blank(),
##           legend.background=element_blank(),
##           axis.title.x=element_blank(),
##           axis.title.y=element_text(size=10),
##           axis.text.x=element_text(size=10, angle=45, hjust=1),
##           axis.text.y=element_text(size=10))

## figfn <- paste(outdir, "Figure2.2_genes_q0.1_olap_treat.bar.png", sep="")
## ggsave(figfn, plot=p2, width=380, height=320, units="px", dpi=120) ## ggsave                                 




#######################################################
### examples of cell-type specific genes
########################################################

rm(list=ls())

outdir <- "./3_RNA_summary.outs/"
outdir2 <- paste(outdir, "Example_genes/", sep="")
if ( !file.exists(outdir2) ) dir.create(outdir2, showWarnings=F, recursive=T)


### 1-lrt results
fn <- "./1_inter_RNA.outs/1.2_lrt_results.rds"
res_lrt <- read_rds(fn)
res_lrt2 <- res_lrt%>%filter(p.adjusted<0.1)%>%arrange(p.value)

## gene1 <- unique(res2$gene)
## opfn <- paste(outdir2, "inter1_lrt_genes.txt", sep="") 
## write.table(unique(res2$gene), file=opfn, row.names=FALSE, col.names=FALSE, quote=FALSE)


### 2-compared effects to 0_CD4Naive 
fn <- "./3_RNA_summary.outs/1_MCls.diff.rds"
res_diff <- read_rds(fn)%>%mutate(padj2=qvalue(p.value)$qvalues)
res_diff2 <- res_diff%>%filter(padj2<0.1)%>%arrange(p.value)


## gene2 <- res_diff%>%filter(p.adjusted<0.1)%>%pull(gene)%>%unique()
## opfn <- paste(outdir2, "inter2_genes.txt", sep="") 
## write.table(gene2, file=opfn, row.names=FALSE, col.names=FALSE, quote=FALSE)

## gene2 <- res_diff%>%filter(padj2<0.1)%>%pull(gene)%>%unique()
## opfn <- paste(outdir2, "inter2_genes_padj2.txt", sep="") 
## write.table(gene2, file=opfn, row.names=FALSE, col.names=FALSE, quote=FALSE)


### 3-results in contrast to control 
fn <- "./1_inter_RNA.outs/2.2_each_pair.results.rds"
DF <- read_rds(fn)


################################
### lrt unique
##################################


outdir2 <- paste(outdir, "Example_genes/lrt_unq/", sep="")
if ( !file.exists(outdir2) ) dir.create(outdir2, showWarnings=F, recursive=T)

geneSel <- setdiff(res_lrt2$gene, res_diff2$gene)


###
### Examples plots
  
for ( i in 1:6){
    
gene0 <- geneSel[i]

### obtain padj from interaction II 
df2 <- res_diff%>%filter(gene==gene0)%>%
    mutate(treats=gsub(".*_", "", conditionX))%>%arrange(p.value)
treat2 <- df2$treats[1]

cat(i, gene0, "\n")

## tmp <- df2%>%filter(grepl(treat2, conditionX))%>%as.data.frame()
    
###
df3 <- df2%>%filter(grepl(treat2, conditionX))%>%
     dplyr::select(condition=conditionX, padj_diff=padj2)
### plot data
plotDF <- DF%>%filter(gene==gene0, grepl(treat2, condition))%>%left_join(df3, by="condition")
x <- str_split(plotDF$condition, "_", simplify=T)

plotDF2 <- plotDF%>%mutate(MCls=paste(x[,1],x[,2], sep="_"), treats=x[,3], MCls_val=as.numeric(x[,1]))%>%
    mutate(bhat=estimate, b_upper=bhat+1.96*stderror, b_lower=bhat-1.96*stderror,
           MCl2=ifelse(padj_diff<0.1&!is.na(padj_diff), glue("<b>{MCls}"), glue("<i>{MCls}")),
           MCl2=fct_reorder(MCl2, MCls_val))

           
col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
  "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
   "6_Monocyte"="#984ea3", "7_dnT"="black")

###
### plots
p <- ggplot(plotDF2, aes(x=bhat, y=MCl2))+
   geom_errorbarh(aes(xmax=b_upper, xmin=b_lower, colour=MCls),
       linewidth=0.5, height=0.2)+
   geom_point(aes(colour=MCls), shape=19, size=1.5)+
   scale_colour_manual(values=col1)+
   geom_vline(aes(xintercept=0), size=0.25, linetype="dashed")+ 
   ## xlab("log odds ratio")+xlim(-3, 6)+
   ggtitle(bquote(italic(.(gene0))~"in"~.(treat2)))+ 
   theme_bw()+
   theme(plot.title=element_text(hjust=0.5, size=9),
         axis.title=element_blank(),
         axis.text.x=element_text(size=10),
         axis.text.y=element_markdown(size=10),
         legend.position="none")
         ## legend.title=element_blank(),
         ## legend.text=element_text(size=8),
         ## legend.key.size=unit(0.4, "cm"))
         ##legend.position="none")
 
figfn <- paste(outdir2, "Figure_", gene0, "_", treat2, ".forest.png", sep="")
ggsave(figfn, p, device="png", width=450, height=500, units="px", dpi=120)

}


    

####
### Examples
i <- 2
gene0 <- geneSel[i]
###
df2 <- res_diff%>%filter(gene==gene0)%>%
    mutate(treats=gsub(".*_", "", conditionX))%>%arrange(p.value)
treat2 <- df2$treats[1]


###
tmp <- res_diff%>%filter(gene==gene0, grepl(treat2, conditionX))%>%as.data.frame()


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

opfn <- paste(outdir2, "Table_", gene0, "_", treat2, ".xlsx", sep="")
write.xlsx(tmp2, file=opfn)


    

##########################################################
### only identified in interaction-II 
#########################################################    


outdir2 <- paste(outdir, "Example_genes/analysis2_unq/", sep="")
if ( !file.exists(outdir2) ) dir.create(outdir2, showWarnings=F, recursive=T)

geneSel <- setdiff(res_diff2$gene, res_lrt2$gene)


###
### Examples plots
  
for ( i in 1:6){
    
gene0 <- geneSel[i]

###    
### obtain padj from interaction II, res_diff 
df2 <- res_diff%>%filter(gene==gene0)%>%
    mutate(treats=gsub(".*_", "", conditionX))%>%arrange(p.value)
treat2 <- df2$treats[1]
cat(i, gene0, treat2, "\n")
    
df3 <- res_diff%>%filter(gene==gene0, grepl(treat2, conditionX))%>%
     dplyr::select(condition=conditionX, padj_diff=padj2)

###    
### plot data from DF, before comparison
plotDF <- DF%>%filter(gene==gene0, grepl(treat2, condition))%>%left_join(df3, by="condition")
x <- str_split(plotDF$condition, "_", simplify=T)

plotDF2 <- plotDF%>%mutate(MCls=paste(x[,1],x[,2], sep="_"),  MCls_val=as.numeric(x[,1]))%>%
    mutate(bhat=estimate, b_upper=bhat+1.96*stderror, b_lower=bhat-1.96*stderror,
           MCl2=ifelse(padj_diff<0.1&!is.na(padj_diff), glue("<b>{MCls}"), glue("<i>{MCls}")),           
           MCl2=fct_reorder(MCl2, MCls_val))

           
col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
  "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
   "6_Monocyte"="#984ea3", "7_dnT"="black")

###
### plots
p <- ggplot(plotDF2, aes(x=bhat, y=MCl2))+
   geom_errorbarh(aes(xmax=b_upper, xmin=b_lower, colour=MCls),
       linewidth=0.5, height=0.2)+
   geom_point(aes(colour=MCls), shape=19, size=1.5)+
   scale_colour_manual(values=col1)+
   geom_vline(aes(xintercept=0), size=0.25, linetype="dashed")+ 
   ## xlab("log odds ratio")+xlim(-3, 6)+
   ggtitle(bquote(italic(.(gene0))~"in"~.(treat2)))+ 
   theme_bw()+
   theme(plot.title=element_text(hjust=0.5, size=9),
         axis.title=element_blank(),
         axis.text.x=element_text(size=10),
         axis.text.y=element_markdown(size=10),
         legend.position="none")
         ## legend.title=element_blank(),
         ## legend.text=element_text(size=8),
         ## legend.key.size=unit(0.4, "cm"))
         ##legend.position="none")
 
figfn <- paste(outdir2, "Figure_", gene0, "_", treat2, ".forest.png", sep="")
ggsave(figfn, p, device="png", width=450, height=500, units="px", dpi=120)

}


###
### tables of examples
i <- 2
gene0 <- geneSel[i]
###
df2 <- res_diff%>%filter(gene==gene0)%>%
    mutate(treats=gsub(".*_", "", conditionX))%>%arrange(p.value)
treat2 <- "caffeine"


###
tmp <- res_diff%>%filter(gene==gene0, grepl(treat2, conditionX))%>%as.data.frame()


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

opfn <- paste(outdir2, "Table_", gene0, "_", treat2, ".xlsx", sep="")
write.xlsx(tmp2, file=opfn)




#################################################
### Examples plots of shared 
################################################


###
### Examples of shared

outdir2 <- paste(outdir, "Example_genes/olap/", sep="")
if ( !file.exists(outdir2) ) dir.create(outdir2, showWarnings=F, recursive=T)

geneSel <- intersect(res_lrt2$gene, res_diff2$gene)

  
for ( i in 1:6){
    
gene0 <- geneSel[i]

    
###    
### obtain padj from interaction II, res_diff 
df2 <- res_diff%>%filter(gene==gene0)%>%
    mutate(treats=gsub(".*_", "", conditionX))%>%arrange(p.value)
treat2 <- df2$treats[1]
cat(i, gene0, treat2, "\n")    
    
df3 <- res_diff%>%filter(gene==gene0, grepl(treat2, conditionX))%>%
     dplyr::select(condition=conditionX, padj_diff=p.adjusted)


###    
### plot data
plotDF <- DF%>%filter(gene==gene0, grepl(treat2, condition))%>%left_join(df3, by="condition")
x <- str_split(plotDF$condition, "_", simplify=T)

plotDF2 <- plotDF%>%mutate(MCls=paste(x[,1],x[,2], sep="_"),  MCls_val=as.numeric(x[,1]))%>%
    mutate(bhat=estimate, b_upper=bhat+1.96*stderror, b_lower=bhat-1.96*stderror,
           MCl2=ifelse(padj_diff<0.1&!is.na(padj_diff), glue("<b>{MCls}"), glue("<i>{MCls}")),           
           MCl2=fct_reorder(MCl2, MCls_val))

           
col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
  "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
   "6_Monocyte"="#984ea3", "7_dnT"="black")

###
### plots
p <- ggplot(plotDF2, aes(x=bhat, y=MCl2))+
   geom_errorbarh(aes(xmax=b_upper, xmin=b_lower, colour=MCls),
       linewidth=0.5, height=0.2)+
   geom_point(aes(colour=MCls), shape=19, size=1.5)+
   scale_colour_manual(values=col1)+
   geom_vline(aes(xintercept=0), size=0.25, linetype="dashed")+ 
   ## xlab("log odds ratio")+xlim(-3, 6)+
   ggtitle(bquote(italic(.(gene0))~"in"~.(treat2)))+ 
   theme_bw()+
   theme(plot.title=element_text(hjust=0.5, size=9),
         axis.title=element_blank(),
         axis.text.x=element_text(size=10),
         axis.text.y=element_markdown(size=10),
         legend.position="none")
         ## legend.title=element_blank(),
         ## legend.text=element_text(size=8),
         ## legend.key.size=unit(0.4, "cm"))
         ##legend.position="none")
 
figfn <- paste(outdir2, "Figure_", gene0, "_", treat2, ".forest.png", sep="")
ggsave(figfn, p, device="png", width=450, height=500, units="px", dpi=120)

}


###

###
### Tables of examples 
i <- 3
gene0 <- geneSel[i]
###
df2 <- res_diff%>%filter(gene==gene0)%>%
    mutate(treats=gsub(".*_", "", conditionX))%>%arrange(p.value)
treat2 <- df2$treats[1]


###
tmp <- res_diff%>%filter(gene==gene0, grepl(treat2, conditionX))%>%as.data.frame()


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

opfn <- paste(outdir2, "Table_", gene0, "_", treat2, ".xlsx", sep="")
write.xlsx(tmp2, file=opfn)




###
### boxplots of stderro


rm(list=ls())

outdir <- "./3_RNA_summary.outs/"


### 1-lrt results
fn <- "./1_inter_RNA.outs/1.2_lrt_results.rds"
res_lrt <- read_rds(fn)
res_lrt2 <- res_lrt%>%filter(p.adjusted<0.1)%>%arrange(p.value)
sigs <- unique(res_lrt2$gene)

## gene1 <- unique(res2$gene)
## opfn <- paste(outdir2, "inter1_lrt_genes.txt", sep="") 
## write.table(unique(res2$gene), file=opfn, row.names=FALSE, col.names=FALSE, quote=FALSE)


### 2-compared effects to 0_CD4Naive 
fn <- "./3_RNA_summary.outs/1_MCls.diff.rds"
res_diff <- read_rds(fn)%>%mutate(padj2=qvalue(p.value)$qvalues)
res_diff2 <- res_diff%>%filter(padj2<0.1)%>%arrange(p.value)
sig2 <- unique(res_diff2$gene)

sig0 <- setdiff(sigs, sig2)


## gene2 <- res_diff%>%filter(p.adjusted<0.1)%>%pull(gene)%>%unique()
## opfn <- paste(outdir2, "inter2_genes.txt", sep="") 
## write.table(gene2, file=opfn, row.names=FALSE, col.names=FALSE, quote=FALSE)

## gene2 <- res_diff%>%filter(padj2<0.1)%>%pull(gene)%>%unique()
## opfn <- paste(outdir2, "inter2_genes_padj2.txt", sep="") 
## write.table(gene2, file=opfn, row.names=FALSE, col.names=FALSE, quote=FALSE)


### 3-results in contrast to control 
fn <- "./1_inter_RNA.outs/2.2_each_pair.results.rds"
DF <- read_rds(fn)


plotDF <- DF%>%mutate(gr=case_when(gene%in%sig0~"inter1_sig", gene%in%sig2~"inter2_sig", TRUE~"Not_sigs"))%>%
   dplyr::select(stderror, gr)%>%filter(!grepl("Not_sigs", gr))


###
p <- ggplot(plotDF, aes(x=gr, y=stderror, color=gr))+
   geom_boxplot(outlier.size=0.8)+
   stat_summary(fun=median, geom="point", shape=23, size=4, stroke=0.8)+
   scale_x_discrete(labels=c("inter1_sig"="inter1_sig(only)", "inter2_sig"="inter2_sig"))+
   ylab("Standard error")+
   theme_bw()+
   theme(legend.position="none",
         axis.title.x=element_blank())

figfn <- paste(outdir, "Figure0_stderr.box.png", sep="")
ggsave(figfn, p, width=380, height=320, units="px", dpi=120)





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
