###
###
library(Matrix)
library(tidyverse)
library(data.table)
## library(irlba)
## library(Seurat)
## library(SeuratDisk)
## library(SeuratData)
## library(Signac)
## library(SeuratWrappers)
## library(SeuratObject)
## library(GenomicRanges)


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


###
###


outdir <- "./3_summary.outs/2_atac_results/"
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
    oneMCl <- MCls[2]
    fn <- paste("./2_DLDA.outs/2_ATAC.", oneMCl, ".DLDA.rds", sep="")
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
        plotDF <- df3[, c("NEW_BARCODE", "MCls", "treat2", sel1, sel2)]
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
                  axis.title=element_text(size=8),
                  axis.text=element_text(size=8))
        p0
    })

    ###
    ###
    figfn <- paste(outdir, "Figure1_LDA_", oneMCl, ".scatters.png", sep="")
    p2 <- plot_grid(plotlist=figs, nrow=3, ncol=5) 
    ggsave(figfn, plot=p2, width=950, height=520, units="px", dpi=120)

    cat(oneMCl, "\n")
}



######################
### Heatmap
######################


for (oneMCl in MCls){
    

fn <- paste("./2_DLDA.outs/2_ATAC.", oneMCl, ".DLDA.rds", sep="")
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
       grid.text(round(mat[i,j],digits=3), x, y, gp=gpar(fontsize=10))
    })

figfn <- paste(outdir, "Figure2_LDA_", oneMCl, "_corr_heatmap.png", sep="")
###ggsave(figfn, p, width=520, height=520, units="px",dpi=120)
png(figfn, height=400, width=450, res=120)
print(p)
dev.off()    

cat(oneMCl, "\n")    
}


