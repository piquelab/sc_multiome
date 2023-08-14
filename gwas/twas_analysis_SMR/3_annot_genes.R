##
library(Matrix)
library(tidyverse)
##library(clusterProfiler)
##library(org.Hs.eg.db)
library(data.table)
library(qvalue)
library(annotables)
library(limma)

##
## library(GenomicRanges)
## library(Seurat)
## library(SeuratDisk)
## library(SeuratData)
## library(Signac)  ##, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
## library(SeuratWrappers)
## library(cicero, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
## library(monocle3)


library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(viridis)
library(openxlsx)
library(cowplot)
library(ggrepel)
##
## options(scipen=18)

rm(list=ls())




####
#### annotation files 
fn <- "../../gtex_v8/torus/torus_input/zzz_torus.annot.gz"
annot <- fread(fn, header=T, data.table=F)

 
traits <- read.table("traits_of_interest.txt")$V1
traits <- sort(traits)

for (trait in traits){


###    
outdir <- paste("./3_examples_output/", trait, "/", sep="")
if ( !file.exists(outdir)) dir.create(outdir, recursive=T)

    
### twas results
fn <- paste("./2_summary_twas/eQTL_dap_conditions/", trait, "_twas.txt", sep="")
res <- read.table(fn, header=T, fill=T)
res2 <- res%>%dplyr::select(Gene=gene2, symbol, id_b38, id_b38_2, 
   pval_eqtl, PIP_conditions, PIP_base,  pval_gwas, FDR_twas=FDR)%>%
   arrange(FDR_twas) ### rank by FDR_twas
    
### colocalization FDR_intact<0.1   
fn <- paste("../enloc_analysis/INTACT_output/", trait, "_intact.txt", sep="")
intact <- read.table(fn, header=T)
intact2 <- intact%>%dplyr::select(Gene, PCG, GLCP, FDR_intact=FDR)

### combine    
res_comb <- res2%>%full_join(intact2, by="Gene")    


### annotation     
annot2 <- annot%>%filter(SNP%in%res_comb$id_b38)
x <- annot2[,-1]
rnz <- rowSums(x)
DF <- data.frame(id_b38=annot2$SNP, nannot=rnz)


res_comb <- res_comb%>%left_join(DF, by="id_b38")
opfn <- paste(outdir, trait, "_comb_infor.txt", sep="")
write.table(res_comb, opfn, quote=F, row.names=F, col.names=T)


### output the genes with FDR_intact<0.1 and also annotated
## res_comb2 <- res_comb%>%filter(nannot>1)
## opfn2 <- paste(outdir, trait, "_intact_sig_annoted.xlsx", sep="")
## write.xlsx(res_comb2, opfn2)

cat(trait, "total:", nrow(res_comb), "annotated:", sum(rnz>0), "\n")
    
}



######################
### summary traits ###
######################

###
### summary traits
summ <- map_dfr(traits, function(trait){
   ##
   fn <- paste("./3_examples_output/", trait, "/", trait, "_comb_infor.txt", sep="")
   res <- read.table(fn, header=T, fill=T)

   ### 
   gene1 <- res%>%dplyr::filter(FDR_twas<0.1)%>%pull(Gene)%>%unique()
   gene2 <- res%>%dplyr::filter(FDR_intact<0.1)%>%pull(Gene)%>%unique() 
   gene3 <- res%>%dplyr::filter(FDR_intact<0.1, nannot>1)%>%pull(Gene)%>%unique()

   df <- data.frame(traits=trait, ntwas=length(gene1), nINTACT=length(gene2), n_annoted=length(gene3))
   df
})

## x <- read.xlsx("../gwas_prepare/traits_of_interest.xlsx")
## names(x)[2] <- "n_ptwas"
## summ2 <- summ%>%left_join(x, by="traits")

opfn <- "./3_examples_output/1_intact_summ.xlsx"
write.xlsx(summ, file=opfn, overwrite=T)

fn <- "./3_examples_output/1_intact_summ.xlsx"
x <- read.xlsx(fn)


##################################################
### Summary-2, which condition response motifs ###
##################################################

fn <- "../../gtex_v8/torus/torus_input/zzz_torus.annot.gz"
annot <- fread(fn, header=T, data.table=F)

rn <- colnames(annot)[-(1:2)]
x <- str_split(rn, "_", simplify=T)
DFcond <- data.frame(comb=rn, MCls=paste(x[,1], x[,2], sep="_"), treats=x[,3])

###
traits <- read.table("traits_of_interest.txt")$V1
traits <- sort(traits)

###
for (trait in traits){    
   cat(trait, "\n") 
   fn <- paste("./3_examples_output/", trait, "/", trait, "_comb_infor.txt", sep="")
   res <- read.table(fn, header=T, fill=T)
   res2 <- res%>%filter(FDR_intact<0.1, nannot>1)%>%dplyr::select(Gene, symbol, SNP=id_b38)
   annot2 <- annot%>%filter(SNP%in%res2$SNP) 
   mat2 <- res2%>%left_join(annot2, by="SNP")
   ##
   outdir <- paste("./3_examples_output/", trait, "/", sep="")
   ##if ( !file.exists(outdir)) dir.create(outdir, recursive=T)
   opfn <- paste(outdir, trait, "_summ_annot.xlsx", sep="") 
   write.xlsx(mat2, file=opfn, overwrite=T)
   ###   
}


###
### summary annotation results


###
traits <- read.table("traits_of_interest.txt")$V1
traits <- sort(traits)
###
trait <- traits[1]
fn <- paste("./3_examples_output/", trait, "/", trait, "_summ_annot.xlsx", sep="")
summ <- read.xlsx(fn)
rn <- colnames(summ)[5:43]
x <- str_split(rn, "_", simplify=T)
DFcond <- data.frame(comb=rn, MCls=paste(x[,1], x[,2], sep="_"), treats=x[,3])


###
### summ for each condition 
summ <- map_dfr(traits, function(trait){
   ###
   cat(trait, "\n") 
   outdir <- paste("./3_examples_output/", trait, "/", sep="")
   fn <- paste(outdir, trait, "_summ_annot.xlsx", sep="")  
   x <- read.xlsx(fn)
   x2 <- x[,c(1,5:43)]
    
   ##
   nx <- ncol(x2) 
   rn <- colnames(x2)[2:nx]
   ngene <- sapply(rn, function(ii){
      ###
      x3 <- x2[,c("Gene",ii)]  
      names(x3) <- c("gene", "Y") 
      gene <- x3%>%filter(Y==1)%>%pull(gene)%>%unique()
      length(gene)
   })
   names(ngene) <- rn 
   df <- data.frame(traits=trait, nannot=length(unique(x2$Gene)))
   df <- cbind(df, t(ngene))
   df
 })   
    
###
opfn <- "./3_examples_output/2.0_summ_conditions.xlsx"
write.xlsx(summ, opfn, overwrite=T)
 

###
### summ for each treatment
treats <- sort(unique(DFcond$treats)) 
summ <- map_dfr(traits, function(trait){
   ###
   cat(trait, "\n") 
   outdir <- paste("./3_examples_output/", trait, "/", sep="")
   fn <- paste(outdir, trait, "_summ_annot.xlsx", sep="")  
   x <- read.xlsx(fn)
   x2 <- x[,c(1,5:43)]
   ##
   ngene <- sapply(treats, function(ii){
      ###
      rn <- DFcond%>%filter(treats==ii)%>%pull(comb) 
      nanno <- rowSums(x2[,rn])
      x3 <- data.frame(gene=x2$Gene, Y=nanno) 
      gene <- x3%>%filter(Y>0)%>%pull(gene)%>%unique()
      length(gene)
   })
   names(ngene) <- treats 
   df <- data.frame(traits=trait, nannot=length(unique(x2$Gene)))
   df <- cbind(df, t(ngene))
   df
 })   
    
###
opfn <- "./3_examples_output/2.1_summ_treats.xlsx"
write.xlsx(summ, opfn, overwrite=T)


###
### summ for each cell types

MCls <- sort(unique(DFcond$MCls)) 
summ <- map_dfr(traits, function(trait){
   ###
   cat(trait, "\n") 
   outdir <- paste("./3_examples_output/", trait, "/", sep="")
   fn <- paste(outdir, trait, "_summ_annot.xlsx", sep="")  
   x <- read.xlsx(fn)
   x2 <- x[,c(1,5:43)]
   ##
   ngene <- sapply(MCls, function(ii){
      ###
      rn <- DFcond%>%filter(MCls==ii)%>%pull(comb) 
      nanno <- rowSums(x2[,rn])
      x3 <- data.frame(gene=x2$Gene, Y=nanno) 
      gene <- x3%>%filter(Y>0)%>%pull(gene)%>%unique()
      length(gene)
   })
   names(ngene) <- MCls 
   df <- data.frame(traits=trait, nannot=length(unique(x2$Gene)))
   df <- cbind(df, t(ngene))
   df
 })   
    
###
opfn <- "./3_examples_output/2.2_summ_MCls.xlsx"
write.xlsx(summ, opfn, overwrite=T)




###
### Proportion of genes annotated in response motifs and heatmap 

outdir <- "./3_examples_output/"

###
###
traitDF <- read.xlsx("traits_name_df.xlsx")

####
fn <- "./3_examples_output/2.0_summ_conditions.xlsx"
x <- read.xlsx(fn)
trait_df <- data.frame(traits=x$traits)%>%
   left_join(traitDF, by=c("traits"="traits_long"))
x2 <- x[,-c(1,2)]
nannot <- x$nannot
x2 <- sweep(x2, 1,  nannot, "/")


### conditions data frame
rn <-  gsub("_d$", "", colnames(x2))
cvt <- str_split(rn, "_", simplify=T)
cvt <- data.frame(comb=rn, MCls=paste(cvt[,1], cvt[,2], sep="_"), treats=cvt[,3])


### plot data
xt <- t(x2)
rownames(xt) <- cvt$comb
colnames(xt)<- trait_df$traits_short


### single colors
olap_min <- min(xt)
olap_max <- max(xt)
mycol <- colorRamp2(seq(0, 1, length.out=100), colorRampPalette(brewer.pal(n=7,name="YlGnBu"))(100))


### annotation row
col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
  "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
   "6_Monocyte"="#984ea3", "7_dnT"="black")
col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
       "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3")

df_col <- data.frame(celltype=cvt$MCls, contrast=cvt$treats)
row_ha <- rowAnnotation(df=df_col, col=list(celltype=col1, contrast=col2),
    annotation_name_gp=gpar(fontsize=10),
    annotation_legend_param=list(
  celltype=list(title_gp=gpar(fontsize=10), labels_gp=gpar(fontsize=8), title="celltype",
                grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm")),
  contrast=list(title_gp=gpar(fontsize=10), labels_gp=gpar(fontsize=8), title="contrast",
                grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm")))) 


###
### main heatmap  
p1 <- Heatmap(xt, col=mycol,
   cluster_rows=F, cluster_columns=F,
   show_row_names=T, row_names_gp=gpar(fontsize=6),
   show_column_names=T, column_names_gp=gpar(fontsize=9),
   ##column_names_rot=-45,   
   left_annotation=row_ha,
   heatmap_legend_param=list(title="Percents",
      title_gp=gpar(fontsize=10),
      at=seq(0, 1, by=0.2), 
      labels_gp=gpar(fontsize=8),
      grid_width=grid::unit(0.45, "cm"),
      legend_height=grid::unit(7.8, "cm")),
   cell_fun=function(j,i,x,y, width, height, fill){
       grid.text(round(xt[i,j], 2), x, y, gp=gpar(fontsize=6))}, 
   use_raster=T, raster_device="png")


###
figfn <- paste(outdir, "Figure3.0_condition_heatmap.png", sep="")
## ggsave(figfn, plot=p2, device="png", width=680, height=700, units="px", dpi=300) ## ggsave 
png(figfn, width=800, height=700,res=120)
p1 <- draw(p1)
dev.off()

res <- as.data.frame(round(xt,3))
res$conditions <- rownames(xt)
opfn <- paste(outdir, "3.0_summ_conditions.xlsx", sep="")
write.xlsx(res, file=opfn, overwrite=T)
### Heatmap



###
### annotation by treats

fn <- "./3_examples_output/2.1_summ_treats.xlsx"
x <- read.xlsx(fn)
trait_df <- data.frame(traits=x$traits)%>%
   left_join(traitDF, by=c("traits"="traits_long"))
x2 <- x[,-c(1,2)]
nannot <- x$nannot
x2 <- sweep(x2, 1,  nannot, "/")
rownames(x2) <- trait_df$traits_short
xt <- as.matrix(x2)

### conditions data frame
cvt <- data.frame(treats=colnames(x2))



### single colors
mycol <- colorRamp2(seq(0, 1, length.out=100), colorRampPalette(brewer.pal(n=7,name="YlGnBu"))(100))


### annotation row
col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
  "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
   "6_Monocyte"="#984ea3", "7_dnT"="black")
col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
       "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3")

df_col <- data.frame(contrast=cvt$treats)
col_ha <- HeatmapAnnotation(df=df_col, col=list(contrast=col2),
    annotation_name_gp=gpar(fontsize=10),
    annotation_legend_param=list(
  ## celltype=list(title_gp=gpar(fontsize=10), labels_gp=gpar(fontsize=8), title="celltype",
  ##               grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm")),
  contrast=list(title_gp=gpar(fontsize=10), labels_gp=gpar(fontsize=8), title="contrast",
                grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm")))) 

 
###
### main heatmap  
p2 <- Heatmap(xt, col=mycol,
   cluster_rows=F, cluster_columns=F,
   show_row_names=T, row_names_gp=gpar(fontsize=9),
   show_column_names=T, column_names_gp=gpar(fontsize=9),
   ##column_names_rot=-45,   
   top_annotation=col_ha,
   heatmap_legend_param=list(title="Percents",
      title_gp=gpar(fontsize=10),
      at=seq(0, 1, by=0.2), 
      labels_gp=gpar(fontsize=8),
      grid_width=grid::unit(0.38, "cm"),
      legend_height=grid::unit(6, "cm")),
   cell_fun=function(j, i, x, y, width, height, fill){
       grid.text(round(xt[i,j], 3), x, y, gp=gpar(fontsize=8))
       }) 
   ## use_raster=T, raster_device="png")


###
figfn <- paste(outdir, "Figure3.1_treats_heatmap.png", sep="")
## ggsave(figfn, plot=p2, device="png", width=680, height=700, units="px", dpi=300) ## ggsave 
png(figfn, width=580, height=620,res=120)
p2 <- draw(p2)
dev.off()

res <- as.data.frame(round(xt,3))
res$traits <- rownames(xt)
opfn <- paste(outdir, "3.1_summ_treats.xlsx", sep="")
write.xlsx(res, file=opfn, overwrite=T)
### Heatmap






###
### annotation by cell types
fn <- "./3_examples_output/2.2_summ_MCls.xlsx"
x <- read.xlsx(fn)
trait_df <- data.frame(traits=x$traits)%>%
   left_join(traitDF, by=c("traits"="traits_long"))
x2 <- x[,-c(1,2)]
nannot <- x$nannot
x2 <- sweep(x2, 1,  nannot, "/")
rownames(x2) <- trait_df$traits_short
xt <- as.matrix(x2)

### conditions data frame
cvt <- data.frame(MCls=colnames(x2))



### single colors

mycol <- colorRamp2(seq(0, 1, length.out=100), colorRampPalette(brewer.pal(n=7,name="YlGnBu"))(100))


### annotation row
col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
  "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
   "6_Monocyte"="#984ea3", "7_dnT"="black")
col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
       "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3")

df_col <- data.frame(celltype=cvt$MCls)
col_ha <- HeatmapAnnotation(df=df_col, col=list(celltype=col1),
    annotation_name_gp=gpar(fontsize=10),
    annotation_legend_param=list(
  celltype=list(title_gp=gpar(fontsize=10), labels_gp=gpar(fontsize=8), title="celltype",
                grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm"))))
  ## contrast=list(title_gp=gpar(fontsize=10), labels_gp=gpar(fontsize=8), title="contrast",
  ##               grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm")))) 


###
### main heatmap  
p3 <- Heatmap(xt, col=mycol,
   cluster_rows=F, cluster_columns=F,
   show_row_names=T, row_names_gp=gpar(fontsize=9),
   show_column_names=T, column_names_gp=gpar(fontsize=9),
   ##column_names_rot=-45,   
   top_annotation=col_ha,
   heatmap_legend_param=list(title="Percents",
      title_gp=gpar(fontsize=10),
      at=seq(0, 1, by=0.2), 
      labels_gp=gpar(fontsize=8),
      grid_width=grid::unit(0.38, "cm"),
      legend_height=grid::unit(6, "cm")),
   cell_fun=function(j,i,x,y, width, height, fill){
       grid.text(round(xt[i,j], 3), x, y, gp=gpar(fontsize=8))}) 


###
figfn <- paste(outdir, "Figure3.2_MCls_heatmap.png", sep="")
## ggsave(figfn, plot=p2, device="png", width=680, height=700, units="px", dpi=300) ## ggsave 
png(figfn, width=600, height=620,res=120)
p3 <- draw(p3)
dev.off()

res <- as.data.frame(round(xt,3))
res$traits <- rownames(xt)
opfn <- paste(outdir, "3.2_summ_MCls.xlsx", sep="")
write.xlsx(res, file=opfn, overwrite=T)
### Heatmap













###
### summary bar plots
traits <- read.table("./traits_ls.txt")$V1

plotDF <- map_dfr(traits, function(ii){
   ###
   fn <- paste("./2_summary_twas/eQTL_conditions/",  ii, "_twas2_sig.xlsx", sep="")
   x1 <- read.xlsx(fn) 
   twas_gene <- unique(x1$gene)
   ##
   df <- data.frame(ntwas=length(twas_gene), traits=ii)
   df 
})    
 
###
my <- max(plotDF$ntwas)+200 
p1 <- ggplot(plotDF, aes(x=traits, y=ntwas))+
    geom_bar(stat="identity", fill="#43a2ca")+
    geom_text(aes(label=ntwas), vjust=0, size=3)+
    ylim(0, my)+
    scale_x_discrete(labels=c("Asthma"="Asthma", "CAD"="CAD", "Height"="Height",
                             "Hypertension"="Hypertension", "UKB_BMI"="BMI"))+
    ylab("#Genes by TWAS")+
    theme_bw()+
    theme(axis.title.x=element_blank(),
          axis.title.y=element_text(size=10),
          axis.text.x=element_text(angle=45, hjust=1, size=10),
          axis.text.y=element_text(size=10))


###
## figfn <- "./3_examples_output/Figure1.0_twas_barplot.png"
## png(figfn, width=480, height=380, res=120)
## p1
## dev.off()


###
### INTACT analysis

plotDF2 <- map_dfr(traits, function(ii){
   ##
    fn <- paste("3_examples_output/", ii, "/", ii, "_intact_sig.xlsx", sep="")
    x <- read.xlsx(fn)
    gene1 <- x%>%filter(nannot>1)%>%pull(Gene)%>%unique()
    n1 <- length(gene1)
    gene2 <- x%>%filter(nannot<=1)%>%pull(Gene)%>%unique()
    n2 <- length(gene2)
    df <- data.frame(ny=c(n1, n2), ntotal=rep(n1+n2, 2), Category=c("gr1", "gr2"), traits=rep(ii,2))
    df
})    

plotDF2 <- plotDF2%>%arrange(traits, desc(Category))%>%
    group_by(traits)%>%
    mutate(label_y=cumsum(ny)-0.5*ny)

###
my <- max(plotDF2$ntotal)+100 
p2 <- ggplot(plotDF2, aes(x=traits, y=ny, fill=Category))+
    geom_bar(stat="identity")+
    geom_text(aes(y=label_y, label=ny), size=3)+
    scale_fill_manual(values=c("gr1"="#f4a582", "gr2"="#bababa"),
         labels=c("gr1"="In Resp motifs", "gr2"="Not Resp motifs"))+
    ylab("#PCGs by INTACT")+    
    ylim(0, my)+
    scale_x_discrete(labels=c("Asthma"="Asthma", "CAD"="CAD", "Height"="Height",
                             "Hypertension"="Hypertension", "UKB_BMI"="BMI"))+
    theme_bw()+
    theme(legend.title=element_blank(),
          legend.position=c(0.3, 0.85),
          ##legend.key=element_blank(),
          legend.key.size=grid::unit(0.4, "cm"),
          legend.background=element_blank(),
          axis.title.x=element_blank(),         
          axis.title.y=element_text(size=10),
          axis.text.x=element_text(angle=45, hjust=1, size=10),
          axis.text.y=element_text(size=10))

figfn <- "./3_examples_output/Figure1.3_comb_barplot.png"
png(figfn, width=680, height=420, res=120)
plot_grid(p1, p2)
dev.off()


    







######################################################
### compare eQTL_base, eQTL_conditions, eQTL_union ###
######################################################

### focus on chromosome 17

## traits <- sort(read.table("traits_ls.txt")$V1)

## trait <- "Asthma"

## ### eQTL_base
## fn <- paste("./2_summary_twas/eQTL_base/", trait, "_twas.txt", sep="")
## x1 <- read.table(fn, header=T)
## x1 <- x1%>%dplyr::filter(FDR<0.1)

## df <- str_split(x1$id_b38, "_", simplify=T)
## df <- df[,1:2]
## colnames(df) <- c("chr", "pos")
## x1 <- cbind(x1, df)

## x1 <- x1%>%filter(chr=="chr17", as.numeric(pos)>39e+6, as.numeric(pos)<41e+6)

## x1 <- x1%>%arrange(gene)%>%
##     dplyr::select(gene, symbol, base_id_b38=id_b38,
##                   base_pval_eqtl=pval_eqtl, base_PIP=PIP, base_pval_gwas=pval_gwas)


## ### eQTL_conditions
## fn <- paste("./2_summary_twas/eQTL_conditions/", trait, "_twas.txt", sep="")
## x2 <- read.table(fn, header=T)
## x2 <- x2%>%dplyr::filter(FDR<0.1)

## df <- str_split(x2$id_b38, "_", simplify=T)
## df <- df[,1:2]
## colnames(df) <- c("chr", "pos")
## x2 <- cbind(x2, df)

## x2 <- x2%>%filter(chr=="chr17", as.numeric(pos)>39e+6, as.numeric(pos)<41e+6)

## x2 <- x2%>%arrange(gene)%>%
##     dplyr::select(gene,  conditions_id_b38=id_b38,
##                   conditions_pval_eqtl=pval_eqtl, conditions_PIP=PIP, conditions_pval_gwas=pval_gwas)


## ### eQTL_union
## fn <- paste("./2_summary_twas/eQTL_union/", trait, "_twas.txt", sep="")
## x3 <- read.table(fn, header=T)
## x3 <- x3%>%dplyr::filter(FDR<0.1)

## df <- str_split(x3$id_b38, "_", simplify=T)
## df <- df[,1:2]
## colnames(df) <- c("chr", "pos")
## x3 <- cbind(x3, df)

## x3 <- x3%>%filter(chr=="chr17", as.numeric(pos)>39e+6, as.numeric(pos)<41e+6)


## x3 <- x3%>%arrange(gene)%>%
##     dplyr::select(gene,  union_id_b38=id_b38,
##                   union_pval_eqtl=pval_eqtl, union_PIP=PIP, union_pval_gwas=pval_gwas)


## ###
## ### combine

## comb <- x1%>%full_join(x2, by="gene")%>%full_join(x3, by="gene")

## ### INTACT results
## fn <- paste("../enloc_analysis/INTACT_output/", trait, "_intact.txt", sep="")
## intact <- read.table(fn, header=T)
## intact <- intact%>%dplyr::select(gene=Gene, conditions_GLCP=GLCP, conditions_PCG=PCG, conditions_FDR_intact=FDR)

## comb <- comb%>%left_join(intact, by="gene")


## opfn <- paste("./3_examples_output/", trait, "_chr17_twas.xlsx", sep="")
## write.xlsx(comb, file=opfn, overwrite=T)


################
### examples ###
################


### annotation

fn <- "../../gtex_v8/torus/torus_input/zzz_torus.annot.gz"
annot <- fread(fn, header=T, data.table=F)


trait <- "Asthma"

fn <- paste("./2_summary_twas/eQTL_conditions/", trait, "_twas2_sig.xlsx", sep="")
x <- read.xlsx(fn)

##coordinates
df <- str_split(x$id_b38, "_", simplify=T)
df <- df[,1:2]
colnames(df) <- c("chr", "pos")
x <- cbind(x, df)

###
### annotation
annot2 <- annot%>%filter(SNP%in%x$id_b38)
rnz <- rowSums(annot2[,-1])
df <- data.frame(id_b38=annot2$SNP, nannot=rnz)

x <- x%>%left_join(df, by="id_b38")

###
### intact results
fn <- paste("../enloc_analysis/INTACT_output/", trait, "_intact.txt", sep="")
intact <- read.table(fn, header=T)
intact <- intact%>%dplyr::select(gene=Gene, GLCP, PCG, FDR_intact=FDR)


###
### combine
comb <- x%>%left_join(intact, by="gene")

opfn <- paste("./2_summary_twas/eQTL_conditions/", trait, "_twas3_sig_comb.xlsx", sep="")
write.xlsx(comb, file=opfn)

### output

comb1 <- comb%>%filter(chr=="chr6")
comb1 <- comb1[1:6,]

comb2 <- comb%>%filter(chr=="chr17", as.numeric(pos)>39e+6, as.numeric(pos)<41e+6)

comb3 <- rbind(comb1, comb2)

comb3 <- comb3%>%dplyr::select(gene, symbol, id_b38, pval_eqtl, PIP, pval_gwas, GLCP, PCG, FDR_intact, nannot)

##
opfn <- paste("./3_examples_output/", trait, "_twas_example.xlsx", sep="")
write.xlsx(comb3, file=opfn)


###
### annotation 
snp <- comb3%>%filter(symbol=="ORMDL3")%>%pull(id_b38)
cvt <- str_split(snp, "_", simplify=T)
chr <- cvt[1,1]
chr_pos2 <- paste(cvt[1,1], cvt[1,2], sep="_")

ann <- annot2%>%filter(SNP==snp)
ann <- unlist(ann[1,3:ncol(ann)])
conditions <- gsub("_d", "", names(ann)[ann==1])


### response motif
fn <- "../../genetic_variant_annot/Response_motif/2.2_response_motif_th0.1.txt"
resp <- read.table(fn, header=T) 
resp2 <- resp%>%filter(comb%in%conditions)

### annotation
fn <- "../../genetic_variant_annot/annot_jaspar2022/zzz_allmotif.bed.gz"
annot <- fread(fn, header=F, data.table=F)

annot2 <- annot%>%filter(V5%in%resp2$motif_ID, V1==chr)%>%mutate(chr_pos=paste(V1, V2, sep="_")) 

annot3 <- annot2%>%filter(chr_pos==chr_pos2)

resp2%>%filter(motif_ID%in%annot3$V5)
