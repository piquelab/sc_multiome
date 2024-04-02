###
###
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
library(GenomicRanges)
library(annotables)

##
library(cowplot)
library(RColorBrewer)
library(ComplexHeatmap)
library(viridis)
library(circlize)
library(ggrepel)
library(ggrastr)
library(ggExtra)
library(openxlsx)


##
rm(list=ls())


###
### sliding windows and visulize dynamic changes
### show example genes 


####################################################
### dynamical changes pseudotime-related genes
###################################################


rm(list=ls())
outdir <- "./3_summary.outs/1_rna_results/"


outdir2 <- "./3_summary.outs/1_rna_results/3.2_dynamic.outs/"
if ( !file.exists(outdir2) ) dir.create(outdir2, showWarnings=F, recursive=T)



###
### sliding window
slideFun <- function(X_sort, lda, win=0.1, step=0.001){
###
  win <- trunc(ncol(X_sort)*win)
  step <- trunc(ncol(X_sort)*step)
    
  X <- X_sort  
  nlen <- ncol(X) 

  Xnew <- NULL
  ldaNew <- NULL  
  s0 <- 1
  while(TRUE){
    ##
    s1 <- s0+win-1
    if (s1>nlen) break
    xi <- apply(X[,s0:s1], 1, mean)
    Xnew <- cbind(Xnew, xi)
    ldai <- mean(lda[s0:s1])
    ldaNew <- c(ldaNew, ldai)  
    s0 <- s0+step
  }
  res <- list(Xnew=Xnew, ldaNew=ldaNew)
  res  
}


###
#####
getDynamicalMat <- function(X, meta, treat, win=0.1, step=0.01){

### re-order cell
  meta_sort <- meta%>%mutate(z2=(z-min(z))/(max(z)-min(z)))%>%arrange(z2)
###     
  x.min <- apply(X, 1, min)
  x.max <- apply(X, 1, max)
  x.r <- x.max-x.min  
  X <- sweep(X, 1, x.min, "-")
  X <- sweep(X, 1, x.r, "/")  
  
## treatment
  treat1 <- treat  
  meta1 <- meta_sort%>%filter(treats==treat1)  
  mat1 <- as.matrix(X[, meta1$barcode])    
  #rownames(mat2) <- gsub("S-|\\..*", "", rownames(mat2))
  tmp1 <- slideFun(mat1, lda=meta1$z2, win=win, step=step)
  mat1 <- tmp1$Xnew
  lda1 <- tmp1$ldaNew
    
  ### contrast treatment
  treat0 <- "control" 
  meta0 <- meta_sort%>%filter(treats==treat0)  
  mat0 <- as.matrix(X[, meta0$barcode])    
  #rownames(mat2) <- gsub("S-|\\..*", "", rownames(mat2))
  tmp0 <- slideFun(mat0, lda=meta0$z2, win=win, step=step)
  mat0 <- tmp0$Xnew
  lda0 <- tmp0$ldaNew  
###
  ## mMax <- apply(mat1, 1, max)
  ## mMin <- apply(mat1, 1, min)
  ## range <- mMax-mMin
  

  ## mat1 <- sweep(mat1, 1, mMin, "-")   
  ## mat1 <- sweep(mat1, 1, range, "/")
  ## mat1 <- as.matrix(mat1)
  ##
  ## mat0 <- sweep(mat0, 1, mMin, "-")
  ## mat0 <- sweep(mat0, 1, mMax, "/")  
  ## mat0 <- sweep(mat0, 1, mMax, "/")
  ## mat0 <- as.matrix(mat0)
###
  mat_ls <- list(mat0=mat0, mat1=mat1, "lda0"=lda0, "lda1"=lda1)  
  mat_ls  
###    
}


### colors
coltreat <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
       "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3", "control"="grey")



####
### Differential results
#### DEGs
fn <- "../sc_multiome_data/2_Differential/1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn_sampleID+new_treat.rds"
res <- read_rds(fn)%>%as.data.frame()%>%mutate(comb=paste(MCls, contrast, sep="_"))
## resDG2 <- resDG%>%filter(abs(estimate)>0.5, p.adjusted<0.1)
###
### select gene
gene0 <- unique(res$gene)
gene1 <- grch38%>%dplyr::filter(chr%in%as.character(1:22))%>%pull(symbol)%>%unique()
gene2 <- intersect(gene0, gene1)

###
geneSel <- res%>%
  filter(abs(estimate)>0.5, p.adjusted<0.1, gene%in%gene2)%>%pull(gene)%>%unique()


MCls <- sort(unique(res$MCls))
for (oneMCl in MCls){       
###
### sc object
fn <- paste("./1_input_data/", oneMCl, ".multiome.rds", sep="")
sc <- read_rds(fn)


###
#### DLDA scores    
fn <- paste("./2_DLDA.outs/1_RNA.", oneMCl, ".DLDA.rds", sep="")
df0 <- read_rds(fn)
rownames(df0) <- NULL
###
df2 <- df0%>%
   pivot_wider(names_from=LDA,  values_from=zscore)%>%
   as.data.frame()




comb <- sort(names(df2)[4:9])
    
###
### 0_CD4Naive_zinc

for (ii in comb){ 

    
   cat(ii,  "\n") 


    
   sc2 <- subset(sc, features=geneSel)
   data2 <- as.matrix(sc2[["SCT"]]$data)


   ### meta data
   treat0 <- gsub(".*_", "", ii)
    
   treatSel <- c(treat0, "control")
   meta <- df2[,c("barcode", "MCls", "treat2",  ii)]%>%
      filter(treat2%in%treatSel)
   names(meta)[3] <- "treats"
   names(meta)[4] <- "z"

   
   ### expression data
   data2 <- data2[, meta$barcode]


   ### sliding windows    
   mat_ls <- getDynamicalMat(data2, meta, treat=treat0, win=0.1, step=0.01)

   opfn <- paste(outdir2, "1_", ii, "_sliding.win.rds", sep="")
   write_rds(mat_ls, file=opfn)


   }
}




############################
### Example genes
#############################


rm(list=ls())

outdir <- "./3_summary.outs/1_rna_results/3.2_dynamic.outs/"
outdir2 <- "./3_summary.outs/1_rna_results/3.2_dynamic.outs/Example2/"
if ( !file.exists(outdir2) ) dir.create(outdir2, showWarnings=F, recursive=T)


###
### setting colors
col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
  "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
   "6_Monocyte"="#984ea3", "7_dnT"="black")
col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
       "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3", "control"="grey")




###
### differential results
fn <- "../sc_multiome_data/2_Differential/1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn_sampleID+new_treat.rds"
res <- read_rds(fn)%>%as.data.frame()%>%
    mutate(comb=paste(MCls, contrast, sep="_"),
           is_sig=ifelse(abs(estimate)>0.5&p.adjusted<0.1, 1, 0),
           is_sig=ifelse(is.na(is_sig), 0, is_sig))


###
### select gene
gene0 <- unique(res$gene)
gene1 <- grch38%>%dplyr::filter(chr%in%as.character(1:22))%>%pull(symbol)%>%unique()
gene2 <- intersect(gene0, gene1)

## sigs <- res%>%filter(is_sig==1, gene%in%gene2)%>%pull(gene)%>%unique()

treats <- sort(unique(res$contrast))


###
### DEGs in treatment across cell types 
treat0 <- treats[6]
sigs <- res%>%filter(is_sig==1, gene%in%gene2, contrast==treat0)%>%pull(gene)%>%unique()

res2 <- res%>%filter(gene%in%sigs, contrast==treat0)%>%dplyr::select(gene, MCls, is_sig, estimate)


df_gene <- res2%>%pivot_wider(id_cols=gene, names_from=MCls, values_from=is_sig)%>%
    as.data.frame()

df_gene$ny <- rowSums(df_gene[,-1])

df2 <- df_gene%>%filter(ny>=4) ##4 for zinc
  

###
### 
fn <- "./3_summary.outs/1_rna_results/3_0_CD4Naive.rr.xlsx"
dfcorr <- read.xlsx(fn)

ii <- which(grepl(treat0, colnames(dfcorr)))
 
dfcorr2 <- dfcorr%>%filter(abs(dfcorr[,ii])>0.2, gene%in%df2$gene)
dfcorr2 <- dfcorr2[,c(ii,7)]
names(dfcorr2)[1] <- "rr"
dfcorr2 <- dfcorr2%>%arrange(desc(rr))
geneSel <- dfcorr2$gene

### caffeine
## geneSel <- c(geneSel[c(1:4, 8, 10:14, 20:23)], "CD55")  
  
MCls <- sort(unique(res$MCls))
for (i in 1:length(geneSel)){
   ##
   gene0 <- geneSel[i]
   cat(gene0, "\n")
   ### plot DF 
   plotDF <- map_dfr(MCls, function(oneMCl){
       ##
       ii <- paste(oneMCl, treat0, sep="_")
       fn <- paste(outdir, "1_", ii,  "_sliding.win.rds", sep="")
       mat_res <- read_rds(fn)

       ###
       mat <- mat_res$mat1
       lda <- mat_res$lda1
       plotDF0 <- data.frame(x=lda, y=as.vector(mat[gene0,]))%>%
            mutate(x2=(x-min(x))/(max(x)-min(x)), treats=treat0)

       ###
       mat <- mat_res$mat0
       lda <- mat_res$lda0
       plotDF2 <- data.frame(x=lda, y=as.vector(mat[gene0,]))%>%
            mutate(x2=(x-min(x))/(max(x)-min(x)), treats="control")
              
       plotDF0 <- rbind(plotDF0, plotDF2)

       plotDF0 <- plotDF0%>%
           mutate(y2=(y-min(y))/(max(y)-min(y)), x2=(x-min(x))/(max(x)-min(x)))
       plotDF0$MCls <- oneMCl
       ###
       plotDF0
   })    
       
   p <- ggplot(plotDF, aes(x=x2, y=y2, color=treats))+
       geom_point(size=0.8)+ 
       geom_smooth(method="loess", se=F, linewidth=0.8)+
       scale_color_manual(values=col2, guide=guide_legend(override.aes=list(size=2)) ) +
       facet_wrap(~MCls, nrow=2)+
       xlab(paste("pseudotime-", treat0, sep=""))+ylab("Relative changes")+
       ggtitle(bquote(~italic(.(gene0))))+
       ylim(0, 1)+
       theme_bw()+
       theme(plot.title=element_text(hjust=0.5, size=12),
           axis.title=element_text(size=10),
           axis.text=element_text(size=9),
           legend.title=element_blank(),
           legend.text=element_text(size=10),
           legend.key=element_blank(),
           legend.background=element_blank(),           
           strip.text=element_text(size=12))
 
figfn <- paste(outdir2, "Figure_", treat0, "_",  i,  "_", gene0, ".fitting.png", sep="")
ggsave(figfn, p, width=1200, height=620, units="px", dpi=120)

}

                     









######
###

## for (oneMCl in MCls){

## ###
## ### read correlation data    
## fn <- paste("./3_summary.outs/3_", oneMCl, ".rr.xlsx", sep="")
## df0 <- read.xlsx(fn)

## ntreat <- ncol(df0)
## geneSel <- lapply(2:ntreat, function(i){
##     ##
##     df2 <- df0[,c(1, i)]
##     names(df2)[2] <- "rr"
##     ###
##     gene0 <- df2%>%arrange(desc(abs(rr)))%>%pull(gene)
##     gene0[1:10]
## })
## geneSel <- unique(unlist(geneSel))



## ###
## ### plot data    
## mat <- df0%>%column_to_rownames(var="gene")
## mat <- as.matrix(mat)
## colnames(mat) <- gsub(".*_", "", colnames(mat))

## mat2 <- mat[geneSel,]


## mycol <- colorRamp2(seq(-1, 1, length.out=50), colorRampPalette(rev(brewer.pal(n=7,name="RdBu")))(50))
## ## mycol <- colorRamp2(seq(0, 1, length.out=12), colorRampPalette(c("white", "red"))(12))


## ###
## ### Heatmap plots
## p <- Heatmap(mat2, name="PCC", na_col="grey90", 
##     col=mycol, cluster_rows=T, cluster_columns=F,
##     show_row_dend=T, row_dend_width=grid::unit(2, "cm"),
##     row_names_gp=gpar(fontsize=8), column_names_gp=gpar(fontsize=10), ##column_names_rot=-45,
##     heatmap_legend_param=list(at=seq(-1, 1, by=0.5),
##         grid_width=grid::unit(0.38, "cm"), legend_height=grid::unit(4, "cm"),
##         title_gp=gpar(fontsize=8), labels_gp=gpar(fontsize=8))
##     )
 
## figfn <- paste(outdir, "Figure3_", oneMCl, "_corr.heatmap.png", sep="")
## ###ggsave(figfn, p, width=520, height=520, units="px",dpi=120)
## png(figfn, height=900, width=450, res=120)
## print(p)
## dev.off()    

## cat(oneMCl, "genes", length(geneSel), "\n")
    
## }    
