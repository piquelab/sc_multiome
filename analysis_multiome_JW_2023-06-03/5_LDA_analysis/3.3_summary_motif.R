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


##
library(cowplot)
library(RColorBrewer)
library(ComplexHeatmap)
library(viridis)
library(circlize)
library(ggrepel)
library(ggrastr)
library(openxlsx)


##
rm(list=ls())

outdir <- "./3_summary.outs/3_motif_results/"
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)


#####################################
### pseudotime related genes
#####################################

######################
### prepare data 
########################

##
rm(list=ls())
outdir <- "./3_summary.outs/3_motif_results/"


###
### extract motif activity data 
fn <- "../sc_multiome_data/3_motif/2_motif.activities.outs/1_scATAC.motifActivities.rds"
sc <- read_rds(fn)
meta2 <- sc@meta.data%>%
    mutate(MCls=paste(wsnn_res.0.1, MCls, sep="_"),
           treat2=ifelse(treats%in%c("etOH", "water"), "control", treats))

sc <- AddMetaData(sc, meta2)
DefaultAssay(sc) <- "ATAC"

data0 <- as.matrix(sc[["chromvar"]]$data)


motif <- unlist(Motifs(sc)@motif.names)
motif_df <- data.frame(motif_ID=names(motif), motif_name=motif)
motifSel <- unique(motif_df$motif_name)

###    
data2 <- lapply(motifSel, function(ii){
    ##
    id2 <- motif_df%>%filter(motif_name==ii)%>%pull(motif_ID)%>%unique()
    d0 <- data0[id2,,drop=F]
    d0 <- apply(d0, 2, mean)
    d0 <- matrix(d0, nrow=1, ncol=length(d0))
    d0
})
data2 <- do.call(rbind, data2)    
colnames(data2) <- colnames(data0)    
rownames(data2) <- motifSel

###
opfn <- paste(outdir, "0_motif_activity.data.rds", sep="")
write_rds(data2, file=opfn)

##
opfn <- paste(outdir, "0_meta.data.rds", sep="")
write_rds(sc@meta.data, file=opfn)




#############################################################
### correlation between motif activity and pseudotime  
#############################################################


##
rm(list=ls())
outdir <- "./3_summary.outs/3_motif_results/"


fn <- "./3_summary.outs/3_motif_results/0_meta.data.rds"
meta <- read_rds(fn)

fn <- "./3_summary.outs/3_motif_results/0_motif_activity.data.rds"
data <- read_rds(fn)
data2 <- t(data)



###
MCls <- sort(unique(meta$MCls))[c(1:2, 4:9)]
for (oneMCl in MCls){
    
#### DLDA scores    
fn <- paste("./2_DLDA.outs/1_RNA.", oneMCl, ".DLDA.rds", sep="")
df0 <- read_rds(fn)
rownames(df0) <- NULL

###
df2 <- df0%>%
   pivot_wider(names_from=LDA,  values_from=zscore)%>%
   as.data.frame()

###
### correlation between pseudotime and gene expression

rnSel <- sort(names(df2)[4:9])
rr_mat <- map_dfc(1:length(rnSel), function(i){
   ##
   rn <- rnSel[i] 
   treatSel <- c(gsub(".*_", "", rn), "control") 
   df3 <- df2%>%dplyr::filter(treat2%in%treatSel)

   cellSel <- df3$barcode 
   data0 <- data2[cellSel,]

   cat(rn, "cell,", identical(rownames(data0), cellSel), " ")
    
   pseu <- df3[,rn]
   rr0 <- cor(data0, pseu)
   colnames(rr0) <- paste("rr_", rn, sep="")
    
   cat("gene,", identical(rownames(rr0), colnames(data2)), "\n")
    
   rr0 
})

##    
rr_df <- rr_mat%>%as.data.frame()
rr_df$motif_name <- colnames(data2)
                     
###
### output 
opfn <- paste(outdir, "3_", oneMCl, "_motif.rr.xlsx", sep="")
write.xlsx(rr_df, file=opfn)

cat(oneMCl, "\n")
###    
}    
    


####################################################
### dynamical changes pseudotime-related genes
###################################################


rm(list=ls())
## outdir <- "./3_summary.outs/3_motif_results/"


outdir2 <- "./3_summary.outs/3_motif_results/3_dynamic.outs/"
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
    xi <- apply(X[,s0:s1], 1, sum)
    Xnew <- cbind(Xnew, xi)
    ldai <- sum(lda[s0:s1])
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
  mMax <- apply(mat1, 1, max)
  mMin <- apply(mat1, 1, min)
  range <- mMax-mMin
  

  mat1 <- sweep(mat1, 1, mMin, "-")   
  mat1 <- sweep(mat1, 1, range, "/")
  mat1 <- as.matrix(mat1)
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



###
### Input data
fn <- "../4_motif_plots.outs/3.2_motif_diff_fitInd_JW.rds"
res <- read_rds(fn)%>%as.data.frame()%>%mutate(comb=paste(MCls, contrast, sep="_"))
res <- res%>%mutate(zscore=beta/stderr)%>%
    dplyr::select(comb, MCls, contrast,  motif_name, beta, stderr, zscore, pval, qval, comb)
res <- res%>%group_by(comb, motif_name)%>%slice_max(order_by=abs(zscore), n=1)%>%as.data.frame()


###
fn <- "../4_motif_plots.outs/response_motif_correct/2.2_response_motif_th0.1.txt"
resp_motif <- read.table(fn, header=T)%>%pull(motif_name)%>%unique()



## fn <- "../sc_multiome_data/3_motif/2_motif.activities.outs/1_scATAC.motifActivities.rds"
## sc <- read_rds(fn)
## meta2 <- sc@meta.data%>%
##     mutate(MCls=paste(wsnn_res.0.1, MCls, sep="_"),
##            treat2=ifelse(treats%in%c("etOH", "water"), "control", treats))

## sc <- AddMetaData(sc, meta2)
## DefaultAssay(sc) <- "ATAC"


## ###
## motif <- unlist(Motifs(sc)@motif.names)
## motif_df <- data.frame(motif_ID=names(motif), motif_name=motif)

## motifSel <- unique(motif_df$motif_name)


fn <- "./3_summary.outs/3_motif_results/0_motif_activity.data.rds"
data <- read_rds(fn)


 

MCls <- sort(unique(res$MCls))
for (oneMCl in MCls[c(1,6,7)]){

###
#### DLDA scores    
fn <- paste("./2_DLDA.outs/1_RNA.", oneMCl, ".DLDA.rds", sep="")
df0 <- read_rds(fn)
rownames(df0) <- NULL
###
df2 <- df0%>%
   pivot_wider(names_from=LDA,  values_from=zscore)%>%
   as.data.frame()


###
### correlation 
fn <- paste("./3_summary.outs/3_motif_results/3_", oneMCl,  "_motif.rr.xlsx", sep="")
rr_df <- read.xlsx(fn)


comb <- sort(names(df2)[4:9])
    
###
### 0_CD4Naive_zinc

for (ii in comb[6]){
    
 
   rr2 <- rr_df%>%dplyr::select(motif_name, paste("rr_", ii, sep=""))
   names(rr2)[2] <- "rr"
    
   
   geneSel2 <- rr2%>%filter(motif_name%in%resp_motif, abs(rr)>0.1)%>%pull(motif_name)
   nsel <- length(geneSel2)
    
   cat(ii, nsel, "\n") 

   if ( nsel<2) next
    

    
   ###
   ### plot data

   ### meta data 
   treatSel <- c(gsub(".*_", "", ii), "control")
   meta <- df2[,c("barcode", "MCls", "treat2",  ii)]%>%
      filter(treat2%in%treatSel)
   names(meta)[3] <- "treats"
   names(meta)[4] <- "z"

   
   ### expression data 
   data2 <- data[geneSel2, meta$barcode]
   treat0 <- gsub(".*_", "", ii)

   ### sliding windows    
   mat_ls <- getDynamicalMat(data2, meta, treat=treat0, win=0.1, step=0.01)

   mat1 <- mat_ls$mat1

   ## if ( nrow(mat1)<5) next 



   #####################
   ### Heatmap 
   #####################

   ### 
   val <- as.vector(mat1)
   q0 <- quantile(val, probs=c(0.025, 0.975))
   break0 <- c(0, seq(q0[1], q0[2], length.out=98), 1) 
   ###setting colors
   colset0 <- colorRampPalette(rev(brewer.pal(11, "Spectral")))(100) 
   col_fun <-  colorRamp2(break0, colset0)

   ###
   ### column annotation
   z <- mat_ls$lda1
   z2 <- (z-min(z))/(max(z)-min(z))
   q0 <- quantile(z2, probs=0.5)
   break2 <- c(seq(0, q0, length.out=50), seq(q0+0.01, 1, length.out=50))
   ### setting colors
   colset2 <- colorRampPalette(c("white", coltreat[[treat0]]))(100)  
   col2 <- colorRamp2(break2, colset2)
      

   col_ha <- HeatmapAnnotation(pseudotime=z2, col=list(pseudotime=col2),
      show_legend=FALSE, simple_anno_size=unit(0.4, "cm"), show_annotation_name=F)
      
      ## annotation_legend_param=list(LFC=list(grid_width=unit(0.3,"cm"),
      ## labels_gp=gpar(fontsize=8),
      ## title_gp=gpar(fontsize=12),
      ## title=paste("pseudotime-",treat0, sep=""))),
      ## annotation_name_gp=gpar(fontsize=10),                         
      ## simple_anno_size=unit(0.3,"cm") )


   ###   
   ### row annotation 
   # <- colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100)
   res0 <- res%>%filter(comb==ii, motif_name%in%geneSel2)%>%select(motif_name, beta) 
   b <- res0$beta
   names(b) <- res0$motif_name
   b <- b[rownames(mat1)]     
   ## th0 <- as.numeric(quantile(abs(b),probs=0.99)) 
   ## b2 <- b[abs(b)<th0] ### quantile(abs(b),probs=0.99)
   ## breaks <- c(min(b),quantile(b2,probs=seq(0,1,length.out=98)),max(b))
   bmax <- max(abs(b)) 
   breaks <- seq(-bmax, bmax, length.out=10) 
   col3 <- colorRamp2(breaks,
      colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(10) ) 
   row_ha <- rowAnnotation(LFC=b, col=list(LFC=col3),
      annotation_legend_param=list(LFC=list(grid_width=unit(0.3,"cm"), legend_height=unit(3, "cm"),                        labels_gp=gpar(fontsize=8),
      title_gp=gpar(fontsize=10), title="motif activity", title_position="leftcenter-rot")),    
      width=unit(0.4,"cm"),
      annotation_name_gp=gpar(fontsize=9), annotation_name_rot=30) 

    ###
    fsize <- case_when(nsel>=60~4,
                       nsel>=45&nsel<60~5,
                       nsel>=20&nsel<45~6,
                       TRUE~7)
    
    p1 <- Heatmap(mat1, col=col_fun,
       cluster_rows=TRUE, cluster_columns=FALSE,
       show_row_dend=T, row_dend_width=grid::unit(1.5, "cm"),
       show_row_names=nsel<100, row_names_gp=gpar(fontsize=fsize),
       show_column_names=FALSE,
       ###
       column_title=paste("pseudotime-", treat0, sep=""),
       column_title_gp=gpar(fontsize=9),
       top_annotation=col_ha,
       right_annotation=row_ha,
       heatmap_legend_param=list(title="Expression",
          title_gp=gpar(fontsize=10),
          title_position="leftcenter-rot",
          labels_gp=gpar(fontsize=8), at=seq(0, 1, by=0.25), 
          grid_width=unit(0.3, "cm"),
          legend_height=unit(4, "cm")),
       use_raster=TRUE, raster_device="png")


    ### output directory
    height0 <- 580  ### ifelse(nsel>50, 620, 580)
    figfn <- paste(outdir2, "Figure_", ii, ".heatmap.png", sep="")
    png(figfn, height=height0, width=560, res=120)
    set.seed(0)
    p1 <- draw(p1)
    ## r.list <- row_order(p1)
    ## r.dend <- row_dend(p2) 
    dev.off()

    ###
    ### genes 
    ## cl <- row_order(p1)
    ## summ <- data.frame(motif_name=rownames(mat1)[cl])%>%left_join(res0, by="motif_name")
    ## opfn <- paste(outdir2, "1_", ii, "_motifname.heatmap.xlsx", sep="")
    ## write.xlsx(summ, opfn)
    }
}

 


###########################################################
### Motif activity dynamic changes along pseudotime     
###########################################################


###
### motif activity correlated with pseudotime 

outdir <- "./3_summary.outs/1_rna_results/"




MCls <- sort(unique(resDG$MCls))
for (oneMCl in MCls){


#### DLDA scores    
fn <- paste("./2_DLDA.outs/1_RNA.", oneMCl, ".DLDA.rds", sep="")
df0 <- read_rds(fn)
rownames(df0) <- NULL

###
df2 <- df0%>%
   pivot_wider(names_from=LDA,  values_from=zscore)%>%
   as.data.frame()


###
### select DEGs for each cell type
geneSel <- resDG%>%dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5, MCls==oneMCl)%>%pull(gene)%>%unique()

###
### seurat object     
fn <- paste("./1_input_data/", oneMCl, ".multiome.rds", sep="")
sc <- read_rds(fn)
sc2 <- subset(sc, features=geneSel)
data2 <- as.matrix(sc2[["SCT"]]$data)
data2 <- t(data2)


###
### correlation between pseudotime and gene expression

rnSel <- sort(names(df2)[4:9])
rr_mat <- map_dfc(1:length(rnSel), function(i){
   ##
   rn <- rnSel[i] 
   treatSel <- c(gsub(".*_", "", rn), "control") 
   df3 <- df2%>%dplyr::filter(treat2%in%treatSel)

   cellSel <- df3$barcode 
   data0 <- data2[cellSel,]

   cat(rn, "cell,", identical(rownames(data0), cellSel), " ")
    
   pseu <- df3[,rn]
   rr0 <- cor(data0, pseu)
   colnames(rr0) <- paste("rr_", rn, sep="")
    
   cat("gene,", identical(rownames(rr0), colnames(data2)), "\n")
    
   rr0 
})

##    
rr_df <- rr_mat%>%as.data.frame()
rr_df$gene <- colnames(data2)

###
### output 
opfn <- paste(outdir, "4_", oneMCl, "_motif.rr.xlsx", sep="")
write.xlsx(rr_df, file=opfn)

cat(oneMCl, "genes:", length(geneSel), "\n")
###    
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
