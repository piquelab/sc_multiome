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
###


outdir <- "./3_summary.outs/1_rna_results/"
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)



######################
### scatter plots
######################


###
### setting colors
col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
  "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
   "6_Monocyte"="#984ea3", "7_dnT"="black")
col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
       "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3", "control"="grey")


###
###
contrast.list <- list("caffeine"=c("caffeine", "control"),
                      "nicotine"=c("nicotine", "control"),
                      "vitA"=c("vitA", "control"),
                      "vitD"=c("vitD", "control"),
                      "vitE"=c("vitE", "control"),
                      "zinc"=c("zinc", "control"))

MCls <- names(col1)
treats <- names(contrast.list)

###
for (oneMCl in MCls){
    ###
    fn <- paste("./2_DLDA.outs/1_RNA.", oneMCl, ".DLDA.rds", sep="")
    df0 <- read_rds(fn)
    rownames(df0) <- NULL

    df2 <- df0%>%
        pivot_wider(names_from=LDA, names_prefix="zscore_", values_from=zscore)%>%
        as.data.frame()

    ###
    comb <- names(df2)[4:9]
    dfcomb <- NULL
    for (i in 1:(length(comb)-1)){
       for ( j in (i+1):length(comb)){
          ##
          x <- data.frame(x1=comb[i], x2=comb[j])
          dfcomb <- rbind(dfcomb,x)
       }    
    }
    

    ###
    ### Figures
    figs <- lapply(1:nrow(dfcomb), function(i){
        ##
        sel1 <- dfcomb$x1[i]
        sel2 <- dfcomb$x2[i]

        treatSel <- c(gsub(".*_", "", sel1), gsub(".*_", "", sel2), "control")        
        df3 <- df2%>%filter(treat2%in%treatSel)
        
        ###
        plotDF <- df3[, c("barcode", "MCls", "treat2", sel1, sel2)]
        names(plotDF)[4] <- "x"
        names(plotDF)[5] <- "y"
        plotDF <- plotDF%>%mutate(gr=ifelse(grepl("control", treat2), "gr0", "gr1")) 
        ##
        p0 <- ggplot(plotDF, aes(x=x, y=y, color=treat2, alpha=gr))+
            geom_point(size=0.1)+
            scale_colour_manual(values=col2)+
            scale_alpha_manual(values=c("gr0"=0.5, "gr1"=1))+
            xlab(gsub(".*_", "", sel1))+ylab(gsub(".*_", "", sel2))+
            theme_bw()+
            theme(legend.position="none",
                  axis.title=element_text(size=9),
                  axis.text=element_text(size=9))
        p0 <- ggMarginal(p0, groupColour=T, groupFill=F, size=4)
    })

    ###
    ###
    figfn <- paste(outdir, "Figure1_LDA_", oneMCl, ".scatters.png", sep="")
    p2 <- plot_grid(plotlist=figs, nrow=3, ncol=5) 
    ggsave(figfn, plot=p2, width=1200, height=720, units="px", dpi=120)

    cat(oneMCl, "\n")
}



######################
### Heatmap
######################


for (oneMCl in MCls){
    

fn <- paste("./2_DLDA.outs/1_RNA.", oneMCl, ".DLDA.rds", sep="")
df0 <- read_rds(fn)
rownames(df0) <- NULL


df2 <- df0%>%
   pivot_wider(names_from=LDA, values_from=zscore)%>%
   as.data.frame()


###
### correlation data
 
rnSel <- sort(names(df2)[4:9])    
mat <- matrix(0, length(rnSel), length(rnSel))
for (i in 1:length(rnSel)){
   ##
   for (j in 1:length(rnSel)){

       rn1 <- rnSel[i]
       rn2 <- rnSel[j]
       treatSel <- c(gsub(".*_", "", rn1), gsub(".*_", "", rn2), "control")
       df3 <- df2%>%filter(treat2%in%treatSel)
       df3 <- df3[,c(rn1, rn2)]
       ### corr
       mat[i,j] <- cor(df3[,1], df3[,2])
   }
}    
     
colnames(mat) <- gsub(".*_", "", rnSel)
rownames(mat) <- gsub(".*_", "", rnSel)


mycol <- colorRamp2(seq(-1, 1, length.out=50), colorRampPalette(rev(brewer.pal(n=7,name="RdBu")))(50))
## mycol <- colorRamp2(seq(0, 1, length.out=12), colorRampPalette(c("white", "red"))(12))



p <- Heatmap(mat, name="PCC", na_col="grey90", 
    col=mycol, cluster_rows=F, cluster_columns=F,
    row_names_gp=gpar(fontsize=10), column_names_gp=gpar(fontsize=10), ##column_names_rot=-45,
    heatmap_legend_param=list(at=seq(-1, 1, by=0.5),
        grid_width=grid::unit(0.38, "cm"), legend_height=grid::unit(4, "cm"),
        title_gp=gpar(fontsize=8), labels_gp=gpar(fontsize=8)),
    cell_fun=function(j, i, x, y, width, height, fill){
          ##
          ## grid.rect(x,y,width, height, gp=gpar(fill=fill, col=fill))         
          grid.text(round(mat[i,j],digits=3), x, y, gp=gpar(fontsize=10))
    }    
    )

figfn <- paste(outdir, "Figure2_LDA_", oneMCl, "_corr_heatmap.png", sep="")
###ggsave(figfn, p, width=520, height=520, units="px",dpi=120)
png(figfn, height=400, width=450, res=120)
print(p)
dev.off()    

cat(oneMCl, "\n")    
}




#####################################
### pseudotime related genes
#####################################

##
rm(list=ls())


outdir <- "./3_summary.outs/1_rna_results/"





###############################################################################################################
### calculate correlation between pseudotime and gene expression across cells in treatment and control  
###############################################################################################################


####
### Differential results
#### DEGs
fn <- "../sc_multiome_data/2_Differential/1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn_sampleID+new_treat.rds"
res <- read_rds(fn)%>%as.data.frame()%>%mutate(comb=paste(MCls, contrast, sep="_"))

###
### select gene
gene0 <- unique(res$gene)
gene1 <- grch38%>%dplyr::filter(chr%in%as.character(1:22))%>%pull(symbol)%>%unique()
gene2 <- intersect(gene0, gene1)




MCls <- sort(unique(res$MCls))
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
geneSel <- res%>%
    dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5, MCls==oneMCl, gene%in%gene2)%>%pull(gene)%>%unique()

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
opfn <- paste(outdir, "3_", oneMCl, ".rr.xlsx", sep="")
write.xlsx(rr_df, file=opfn)

cat(oneMCl, "genes:", length(geneSel), "\n")
###    
}    
    


####################################################
### dynamical changes pseudotime-related genes
###################################################


rm(list=ls())
outdir <- "./3_summary.outs/1_rna_results/"


outdir2 <- "./3_summary.outs/1_rna_results/3_dynamic.outs/"
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


###
### correlation 
fn <- paste("./3_summary.outs/1_rna_results/3_", oneMCl,  ".rr.xlsx", sep="")
rr_df <- read.xlsx(fn)


comb <- sort(names(df2)[4:9])
    
###
### 0_CD4Naive_zinc

for (ii in comb){


    
   ### 
   geneSel <- res%>%
       filter(comb==ii, abs(estimate)>0.5, p.adjusted<0.1, gene%in%gene2)%>%pull(gene)
   rr2 <- rr_df%>%filter(gene%in%geneSel)%>%dplyr::select(gene, paste("rr_", ii, sep=""))
   names(rr2)[2] <- "rr"
     
   
   geneSel2 <- rr2%>%filter(abs(rr)>0)%>%pull(gene)
   nsel <- length(geneSel2)
    
   cat(ii, nsel, "\n") 

   if ( nsel<2) next

    
    
   sc2 <- subset(sc, features=geneSel2)
   data2 <- as.matrix(sc2[["SCT"]]$data)


    
   ###
   ### plot data

   ### meta data 
   treatSel <- c(gsub(".*_", "", ii), "control")
   meta <- df2[,c("barcode", "MCls", "treat2",  ii)]%>%
      filter(treat2%in%treatSel)
   names(meta)[3] <- "treats"
   names(meta)[4] <- "z"

   
   ### expression data
   data2 <- data2[, meta$barcode]
   treat0 <- gsub(".*_", "", ii)

   ### sliding windows    
   mat_ls <- getDynamicalMat(data2, meta, treat=treat0, win=0.1, step=0.01)

   ## opfn <- paste(outdir2, "1_", ii, "_sliding.win.rds", sep="")

   ## write_rds(mat_ls, file=opfn)
    
   mat1 <- mat_ls$mat1 



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
   res0 <- res%>%filter(comb==ii, gene%in%geneSel2)%>%select(gene, estimate) 
   b <- res0$estimate
   names(b) <- res0$gene
   b <- b[rownames(mat1)]     
   ## th0 <- as.numeric(quantile(abs(b),probs=0.99)) 
   ## b2 <- b[abs(b)<th0] ### quantile(abs(b),probs=0.99)
   ## breaks <- c(min(b),quantile(b2,probs=seq(0,1,length.out=98)),max(b))
   bmax <- max(abs(b)) 
   breaks <- seq(-bmax, bmax, length.out=10) 
   col3 <- colorRamp2(breaks,
      colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(10) ) 
   row_ha <- rowAnnotation(LFC=b, col=list(LFC=col3),
      annotation_legend_param=list(LFC=list(grid_width=unit(0.3,"cm"), legend_height=unit(3, "cm"),
      labels_gp=gpar(fontsize=8),
      title_gp=gpar(fontsize=10), title="LFC", title_position="leftcenter-rot")),    
      width=unit(0.4,"cm"),
      annotation_name_gp=gpar(fontsize=9), annotation_name_rot=30) 

    fsize <- case_when(nsel>=60~4,
                       nsel>=45&nsel<60~5,
                       nsel>=20&nsel<45~6,
                       TRUE~7)
    ###
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
    height0 <- 580 ##ifelse(nsel>50, 620, 580)
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
    ## summ <- data.frame(gene=rownames(mat1)[cl])%>%left_join(res0, by="gene")
    ## opfn <- paste(outdir2, "1_", ii, "_genes.heatmap.xlsx", sep="")
    ## write.xlsx(summ, opfn)
    }
}





                     


