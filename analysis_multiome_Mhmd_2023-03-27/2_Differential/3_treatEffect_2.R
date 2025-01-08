##
library(Matrix)
library(tidyverse)
library(DESeq2)
library(annotables)
##
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(ggExtra)
library(gtable)
library(ggsignif)
library(pheatmap)
library(ComplexHeatmap)
library(corrplot)
library(gtable)
library(RColorBrewer)
library(viridis)
library(ggrastr)
library(ggsci)
library(circlize)

rm(list=ls())

##
outdir <- "./3_treatEffect_2.outs/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=F)


#####################################################################
###            heatmap for resDP                                  ###
#####################################################################
## res <- read_rds("./1.2_DiffPeak.outs/2.0_DESeq.results.rds")%>%
##    as.data.frame()%>%
##    mutate(comb=paste(MCls, contrast, sep="_"))%>%
##    dplyr::filter(MCls!="DC") 

## DP <- res%>%drop_na(p.adjusted)%>%
##    dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)%>%
##    dplyr::pull(gene)
## DP <- as.character(unique(DP))


#----resDP----     
resDP <- read_rds("../2_Differential/1.2_DiffPeak.outs/2.0_DESeq.results_0.1_cn.rds") %>% as.data.frame()

nrow(resDP)

nrow(resDP %>% dplyr::filter(is.na(p.adjusted)))

#resDP <- resDP %>% mutate(MCls=gsub("_", "-", MCls)) %>% mutate(comb=paste(MCls, contrast, sep="_"))%>% dplyr::rename("peak"="gene")
#%>% drop_na(p.value)
resDP <- resDP %>% mutate(MCls=gsub("_", "-", MCls)) %>% mutate(comb=paste(contrast, MCls, sep="_"))%>% dplyr::rename("peak"="gene")
#%>% drop_na(p.value)
head(resDP)
nrow(resDP)
max(resDP$estimate)
min(resDP$estimate)


resDP[!is.na(resDP$p.adjusted) & resDP$p.adjusted<0.1 & resDP$estimate > 0.5, "Significance"] <- "Up"
resDP[!is.na(resDP$p.adjusted) & resDP$p.adjusted<0.1 & resDP$estimate < -0.5, "Significance"] <- "Down"
resDP[!is.na(resDP$p.adjusted) & resDP$p.adjusted>=0.1, "Significance"] <- "Not Significant"
resDP[is.na(resDP$p.adjusted),"Significance"] <- "NA"
resDP <- resDP %>% mutate(comb_2=paste(comb, Significance, sep="_"))
head(resDP)
table(resDP$Significance)

table(resDP$MCls, resDP$contrast)



all.DP  <- resDP %>% dplyr::pull(peak)
all.DP <- as.character(unique(all.DP))
head(all.DP)
length(all.DP)


#----sig_resDP----
sig_resDP <- resDP %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<0.1, abs(estimate)>0)
head(sig_resDP)
nrow(sig_resDP)
max(sig_resDP$estimate)
min(sig_resDP$estimate)

sig <- table(sig_resDP$MCls, sig_resDP$contrast)

sig

#to check
nrow(sig_resDP %>% dplyr::filter(abs(estimate)>5))

#peakAll <- unique(sig_resDP$peak)
#head(peakAll)
#length(peakAll)

DP  <- sig_resDP %>% dplyr::pull(peak)
DP <- as.character(unique(DP))
head(DP)
length(DP)



#----resDP significant genes table----
library(reshape2)

sig.matrix <- as.matrix(sig)
sig.melt <- melt(sig.matrix)


sig.p <- ggplot(sig.melt, aes(x = Var2, y = Var1)) +
          geom_tile(aes(fill=value), color="black")+
          scale_fill_gradient(low = "white", high = "white")+
          geom_text(aes(label = value), color = "black", size = 5) +
          ggtitle("Significant_0.1_0.1-0.5") +
          theme(axis.text.x = element_text(angle = 90),
                          axis.title.x = element_blank(),
                          axis.title.y = element_blank(),
                          legend.position="none")



library("ggpubr")

figure <- ggarrange(sig.p +
                    font("x.text", size = 14))
                   # ncol = 2, nrow = 2)

#png(paste0("./3_treatEffect.outs/plot_sig_treatandind_filt_all_0.1_0.1_0.5_cn.png"), width=500, height=500, pointsize=16, res=125)
png(paste0("./3_treatEffect.outs/plot_sig_treatandind_filt_DP_all_0.1_0.1_0_cn_mitofilt.png"), width=500, height=500, pointsize=16, res=125)
print(figure)
dev.off()




#################
#----mat data----
#################
comb <- unique(resDP$comb)
comb


#----build matrix for heatmap for all genes----
mat <- map_dfc(comb, function(ii){
##
  beta0 <- rep(NA, length(all.DP))
  names(beta0) <- all.DP
  d0 <- resDP%>%dplyr::filter(comb==ii, peak%in%DP)
  peakSel <- as.character(d0$peak)
  beta0[peakSel] <- d0$estimate
  beta0
})
head(mat)

mat <- as.matrix(mat)
rownames(mat) <- all.DP
colnames(mat) <- comb
head(mat)
max(mat)
min(mat)

opfn <- paste0("./3_treatEffect.outs/mat_peaks_macs2_0.1_cn.rds")
write_rds(mat, opfn)


#----build matrix for heatmap for  only significant----
mat <- map_dfc(comb, function(ii){
##
  beta0 <- rep(NA, length(DP))
  names(beta0) <- DP
  d0 <- resDP%>%dplyr::filter(comb==ii, peak%in%DP)
  peakSel <- as.character(d0$peak)
  beta0[peakSel] <- d0$estimate
  beta0
})
head(mat)

mat <- as.matrix(mat)
rownames(mat) <- DP
colnames(mat) <- comb
head(mat)
max(mat)
min(mat)


#----build matrix for heatmap for top 10 significant only----
comb

head(sig_resDP)

mat.10topdps.qval <- map_dfr(comb, function(ii){
      # res2 <- res.topmotif %>% dplyr::filter(condition==ii)
       res2 <- sig_resDP %>% dplyr::filter(comb==ii)
       res2 <- res2 %>% arrange(p.adjusted)
       res2 <- res2 [1:10, ]
       #b <- res2$beta
       #names(b) <- res2$motif
       #b[topmotif]
    })
head(mat.10topdps.qval)
nrow(mat.10topdps.qval)
table(mat.10topdps.qval$MCls, mat.10topdps.qval$contrast)
#length(unique(mat.5topmotifs.qval$motif))

topdps.qval <- unique(mat.10topdps.qval$peak)
head(topdps.qval)
length(topdps.qval)


mat.10topdps.qval.2 <- map_dfc(comb, function(ii){
        res2 <- resDP %>% dplyr::filter(comb==ii)
            #res2 <- res2 %>% arrange(desc(qval))
            #res2 <- res2[1:5, ]
            #res2
            b <- res2$estimate
            names(b) <- res2$peak
            b[topdps.qval]
            })
head(mat.10topdps.qval.2)
nrow(mat.10topdps.qval.2)




mat <- mat.10topdps.qval.2



## mat <- map_dfc(comb, function(ii){
## ##
##   beta0 <- rep(NA, length(DP))
##   names(beta0) <- DP
##   d0 <- resDP%>%dplyr::filter(comb==ii, peak%in%DP) %>% dplyr::arrange(p.adjusted) %>% head(n=10)
##   peakSel <- as.character(d0$peak)
##   beta0[peakSel] <- d0$estimate
##   beta0
## })
## head(mat)

mat <- as.matrix(mat)
rownames(mat) <- DP
colnames(mat) <- comb
head(mat)
max(mat)
min(mat)




#----sum filt----
ii <- rowSums(is.na(mat))
ii > 0
head(ii)
length(ii)

#ii <- ii < 56
#head(ii)
#length(ii)
#table (ii)

#EITHER
mat2 <- mat[ii==0,]

#OR
mat2 <- mat[ii<length(comb),]


#mat2[!is.na(mat2)] = 1
mat2[is.na(mat2)] = 0
#head(mat2)
#nrow(mat2)
#max(mat2)
#min(mat2)
#mean(mat2)
#quantile(mat2,probs=seq(0,1,length.out=101))



#OR
mat[is.na(mat)] = 0

mat2 <- mat
head(mat2)


head(mat2)
nrow(mat2)
max(mat2)
min(mat2)
mean(mat2)


comb.2 <- sort(comb)
comb.2

mat3 <- mat2[, comb.2]
head(mat3)




##############
#----plot----
##############
#----breaks and color scale----
y <- as.numeric(mat2)
head(y)
min(y)
max(y)
mean(y)

quantile(abs(y),probs=0.99) 

#y0 <- y[abs(y)<1.8] #99% percent quantile(abs(y),probs=0.99)
#head(y0)

y0 <- y[abs(y)< quantile(abs(y),probs=0.99)] #99% percent quantile(abs(y),probs=0.99)
head(y0)

quantile(y0,probs=seq(0,1,length.out=98))

mybreaks <- c(min(y),quantile(y0,probs=seq(0,1,length.out=98)),max(y))
names(mybreaks) <- NULL
head(mybreaks)

mycol <- colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100)


#----annotation----
x <- str_split(comb, "_", simplify=T)
x


x <- str_split(comb.2, "_", simplify=T)
x


tmp_column <- data.frame(treatment=x[,1], celltype=x[,2])
rownames(tmp_column) <- comb.2
head(tmp_column)


## col1 <- c("LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
##    "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
## col2 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
##    "NKcell"="#aa4b56", "Tcell"="#ffaa00")
col1 <- c("caffeine"="red", "nicotine"="tan",
            "vitA"="tan4", "vitD"="seagreen4",
            "vitE"="salmon3", "water"="grey", "zinc"="maroon3")
col2 <- c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
           "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00", "1-TCM"="pink",
           "3-TEM"="blue", "5-CD8Naive"="green", "7-dnT"="black")

tmp_colors <- list(treatment=col1, celltype=col2) #brewer.pal(4,"Set1")
#tmp_colors

#mycol <- viridisLite::viridis(100)
#mycol <- viridisLite::cividis(100, direction=1)


#----plot_sig----
## fig1 <- pheatmap(mat2, col=mycol, breaks=mybreaks, border_color="NA",
##    cluster_rows=T, cluster_cols=F,
##    annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend=T,
##    show_colnames=T, show_rownames=F, na_col="white",
##    fontsize_row=12)

#fig1 <- pheatmap(mat2,  col=mycol, breaks=mybreaks, border_color="NA",
fig1 <- pheatmap(mat3,  col=mycol, breaks=mybreaks, border_color="NA",
   cluster_rows=T, cluster_cols=F,
   annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend=T,
   show_colnames=T,
   show_rownames=T,
   na_col="white",
   fontsize_col=10,
   fontsize_row=10)


figfn <- "./3_treatEffect.outs/Figure1.1_heatmap.beta.DP.sig.treat.png"
png(figfn, width=1800, height=2100,res=225)
print(fig1)
dev.off() 


#----plot_top10----
#fig1 <- pheatmap(mat2, border_color="NA",
fig1 <- pheatmap(mat3, border_color="NA",
   cluster_rows=F, cluster_cols=F,
   annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend=T,
   show_colnames=T,
   show_rownames=T,
   na_col="white",
   fontsize_col=5,
   fontsize_row=2)


figfn <- "./3_treatEffect.outs/Figure1.1_heatmap.beta.DP.top10.allvalues.treat.png"
png(figfn, width=1800, height=2100,res=225)
print(fig1)
dev.off() 




#######################################
#----plotting as done in motif file----
#######################################
#----annotation as cell type----
x <- str_split(comb, "_", simplify=T)
x
column_ha_ct <- HeatmapAnnotation(
        treatment=x[,1],
        celltype=x[,2],
        col=list(treatment=c("caffeine"="red", "nicotine"="tan",
                             "vitA"="tan4", "vitD"="seagreen4",
                             "vitE"="salmon3", "water"="grey", "zinc"="maroon3"),
                 celltype=c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
                            "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00",
                            "1-TCM"="pink", "3-TEM"="blue",
                            "5-CD8Naive"="green", "7-dnT"="black"))
)


#----annotation as treatment----
x <- str_split(comb.2, "_", simplify=T)
x
column_ha_treat <- HeatmapAnnotation(
        treatment=x[,1],
        celltype=x[,2],
        col=list(treatment=c("caffeine"="red", "nicotine"="tan",
                             "vitA"="tan4", "vitD"="seagreen4",
                             "vitE"="salmon3", "water"="grey", "zinc"="maroon3"),
                 celltype=c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
                            "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00",
                            "1-TCM"="pink", "3-TEM"="blue",
                            "5-CD8Naive"="green", "7-dnT"="black"))
)


#----breaks and color----
#breaks <- quantile(mat2,probs=seq(0,1,length.out=101))

#col_fun <-

#breaks <- c(seq(0,1,length.out=50), quantile(b2, probs=seq(0, 1, length.out=49)), 38.93)
#breaks
#col_fun <-  colorRamp2(breaks,
#   colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100))


#quantile(mat2, probs=seq(0, 1, length.out=101))

mat2.vector <- as.vector(mat2)
mat2.vector <- sort(mat2.vector)
mat2.vector

breaks <- (seq(-8, 8, length.out=17))
breaks


#col_fun <-  colorRamp2(breaks,
#   colorRampPalette(rev(brewer.pal(n=11, name="RdBu")))(13))

library(ggsci)
library(circlize)

col_fun <- colorRamp2(breaks, colorRampPalette(rev(brewer.pal(n=11, name="RdBu")))(17))


#----figure as celltype----
fig <- Heatmap(mat2,
               col=col_fun,
               cluster_rows=T, cluster_columns=F,
               top_annotation=column_ha_ct,
               heatmap_legend_param=list(title="fold.enrichment",
                                         title_gp=gpar(fontsize=10),
                                         labels_gp=gpar(fontsize=10)),
               show_row_names=T,
               show_column_names=T,
               row_names_gp=gpar(fontsize=2),
               column_names_gp=gpar(fontsize=5),
               raster_device="png")

#figfn <- "./3_treatEffect.outs/Figure1.1.2_heatmap.beta.DP.sig.cell.png"
figfn <- "./3_treatEffect.outs/Figure1.1.2_heatmap.beta.DP.top10.cell_2.png"
png(figfn, width=1800, height=2100,res=225)
print(fig)
dev.off()



#----figure as treatment---- 
fig <- Heatmap(mat3,
#               col=col_fun,
               cluster_rows=T, cluster_columns=F,
               top_annotation=column_ha_treat,
               heatmap_legend_param=list(title="fold.enrichment",
                                         title_gp=gpar(fontsize=10),
                                         labels_gp=gpar(fontsize=10)),
               show_row_names=T,
               show_column_names=T,
               row_names_gp=gpar(fontsize=2),
               column_names_gp=gpar(fontsize=5),
               raster_device="png")

#figfn <- "./3_treatEffect.outs/Figure1.1.2_heatmap.beta.DP.sig.treat.png"
figfn <- "./3_treatEffect.outs/Figure1.1.2_heatmap.beta.DP.top10.allvalues.treat.png"
png(figfn, width=1800, height=2100,res=225)
print(fig)
dev.off()

end
#--------




###########################
### correlation heatmap ###
###########################
#Neworder <- c("Monocyte_LPS", "Monocyte_PHA", "Bcell_LPS", "Bcell_PHA",
   "NKcell_LPS", "NKcell_PHA", "Tcell_LPS", "Tcell_PHA",
   "Monocyte_LPS-DEX", "Monocyte_PHA-DEX", "Bcell_LPS-DEX", "Bcell_PHA-DEX",
   "NKcell_LPS-DEX", "NKcell_PHA-DEX", "Tcell_LPS-DEX", "Tcell_PHA-DEX")

Neworder <- c(comb.2)
Neworder

corr <- cor(mat3)[Neworder, Neworder]
head(corr)

mycol <- colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100)
##

###
x <- str_split(colnames(corr), "_", simplify=T)

tmp_column <- data.frame(celltype=x[,2], treatment=x[,1])
rownames(tmp_column) <- colnames(corr)

tmp_colors <- list(celltype=col2, treatment=col1)

fig2 <- pheatmap(corr, col=mycol, scale="none", border_color="NA",
   cluster_rows=F, cluster_cols=F,
   annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend =T,
   show_colnames=T, show_rownames=T,
   fontsize_row=5,
   fontsize_col=5,
   na_col="white")

###
figfn <- "./3_treatEffect.outs/Figure1.2_heatmap.corr.DP.sig.treat.png"
png(figfn, width=1800, height=1800,res=225)
print(fig2)
dev.off() 























#####################################################################
#### resDE heatmap   ###
#####################################################################
#----resDE----     
resDE <- read_rds("../2_Differential/1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn.rds") %>% as.data.frame()

#resDE <- resDE %>% mutate(MCls=gsub("_", "-", MCls)) %>% mutate(comb=paste(MCls, contrast, sep="_"))
#%>% drop_na(p.value)

resDE <- resDE %>% mutate(MCls=gsub("_", "-", MCls)) %>% mutate(comb=paste(contrast, MCls, sep="_"))
#%>% drop_na(p.value)
head(resDE)
nrow(resDE)
max(resDE$estimate)
min(resDE$estimate)


#resDE$Significance <- NULL
#head(resDE)

resDE[!is.na(resDE$p.adjusted) & resDE$p.adjusted<0.1 & resDE$estimate > 0.5, "Significance"] <- "Up"
resDE[!is.na(resDE$p.adjusted) & resDE$p.adjusted<0.1 & resDE$estimate < -0.5, "Significance"] <- "Down"
resDE[!is.na(resDE$p.adjusted) & resDE$p.adjusted>=0.1, "Significance"] <- "Not Significant"
resDE[is.na(resDE$p.adjusted),"Significance"] <- "NA"
head(resDE)
nrow(resDE)
table(resDE$Significance)


resDE[!is.na(resDE$p.adjusted) & resDE$p.adjusted<0.1 & resDE$estimate > 0.25, "Significance.25"] <- "Up"
resDE[!is.na(resDE$p.adjusted) & resDE$p.adjusted<0.1 & resDE$estimate < -0.25, "Significance.25"] <- "Down"
resDE[!is.na(resDE$p.adjusted) & resDE$p.adjusted>=0.1, "Significance.25"] <- "Not Significant"
resDE[is.na(resDE$p.adjusted),"Significance.25"] <- "NA"
#head(resDE)
table(resDE$Significance.25)


resDE[!is.na(resDE$p.adjusted) & resDE$p.adjusted<0.1 & resDE$estimate > 0, "Significance.0"] <- "Up"
resDE[!is.na(resDE$p.adjusted) & resDE$p.adjusted<0.1 & resDE$estimate < 0, "Significance.0"] <- "Down"
resDE[!is.na(resDE$p.adjusted) & resDE$p.adjusted>=0.1, "Significance.0"] <- "Not Significant"
resDE[is.na(resDE$p.adjusted),"Significance.0"] <- "NA"
#head(resDE)
table(resDE$Significance.0)


resDE <- resDE %>% mutate(comb_2=paste(comb, Significance, sep="_"))
head(resDE)
table(resDE$Significance)
length(unique(resDE$gene))


resDE <- resDE %>% mutate(contrast.gene=paste0(contrast, ".", gene))
head(resDE)
nrow(resDE)


#dplyr::filter does not work for vector
#mit.genes <- resDE$gene %>% dplyr::filter(substr(resDE$gene, 1, 2)=="MT") 
#the foll works but dont need as doing it other way
#mit.genes <- resDE$gene[substr(resDE$gene, 1, 3)=="MT-"] 
#mit.genes

#this select all the rows with "MT" present anywhere in string
#library(data.table)
#mit.resDE <- resDE[resDE$gene %like% "MT", ]
#mit.resDE$gene

resDE <- resDE %>% mutate(mito=ifelse(substr(resDE$gene, 1, 3)=="MT-", 1, 0))
#head(resDE)
nrow(resDE)
table(resDE$mito)


resDE <- resDE %>% mutate(z_score=estimate/stderror)
head(resDE)
nrow(resDE)


#----mitofilt resDE----
resDE <- resDE %>% dplyr::filter(resDE$mito==0)
#head(resDE)
nrow(resDE)
length(unique(resDE$gene))

table(resDE$mito)

table(resDE$Significance)
table(resDE$Significance.25)
table(resDE$Significance.0)




#---------













#----all.DE----
all.DE  <- resDE %>% dplyr::pull(gene)
all.DE <- as.character(unique(all.DE))
head(all.DE)
length(all.DE)







## #----count genes for diff treatments and cell types----
## ## resDE.caff <- resDE %>% dplyr::filter(contrast=="caffeine")
## ## head(resDE.caff)
## ## nrow(resDE.caff)
## ## table(resDE.caff$Significance)
## ## table(resDE.caff$Significance.25)


## ## resDE.caff.sum <- resDE.caff %>% arrange(p.adjusted)
## ## resDE.caff.sum <- resDE.caff.sum %>% distinct(gene, .keep_all = TRUE)
## ## #%>% summarise_each(funs(mean))
## ## head(resDE.caff.sum)
## ## nrow(resDE.caff.sum)
## ## table(resDE.caff.sum$Significance)
## ## table(resDE.caff.sum$Significance.25)

## head(resDE)
## nrow(resDE)


## MCls <- unique(resDE$MCls)
## MCls

## treats <- unique(resDE$contrast)
## treats


## for(i in MCls){
##     for(j in treats){
##         resDE2 <- resDE %>% dplyr::filter(contrast==j, MCls==i)
##         print(paste0(i, "_", j))
##         print(nrow(resDE2))
##         print(table(resDE2$Significance))
##         print(table(resDE2$Significance.25))
## }}        
        






## end







## #----sig_resDE----
sig_resDE.5 <- resDE %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)
sig_resDE.25 <- resDE %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<0.1, abs(estimate)>0.25)
sig_resDE.0 <- resDE %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<0.1, abs(estimate)>0)

sig <- table(sig_resDE.5$MCls, sig_resDE.5$contrast)


sig <- table(sig_resDE.25$MCls, sig_resDE.25$contrast)

sig <- table(sig_resDE.0$MCls, sig_resDE.0$contrast)




end


## sig_resDE <- resDE %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)
## head(sig_resDE)
## nrow(sig_resDE)
## max(sig_resDE$estimate)
## min(sig_resDE$estimate)
## length(unique(sig_resDE$gene))

## sig <- table(sig_resDE$MCls, sig_resDE$contrast)
## sig

## table(sig_resDE$MCls, sig_resDE$Significance)
## table(sig_resDE$contrast, sig_resDE$Significance)

## #to check
## nrow(sig_resDE %>% dplyr::filter(abs(estimate)>5))



#----resDE significant genes table----
library(reshape2)

sig.matrix <- as.matrix(sig)
sig.melt <- melt(sig.matrix)


sig.p <- ggplot(sig.melt, aes(x = Var2, y = Var1)) +
          geom_tile(aes(fill=value), color="black")+
          scale_fill_gradient(low = "white", high = "white")+
          geom_text(aes(label = value), color = "black", size = 5) +
          ggtitle("Significant_0.1_0.1-0.5") +
          theme(axis.text.x = element_text(angle = 90),
                          axis.title.x = element_blank(),
                          axis.title.y = element_blank(),
                          legend.position="none")



library("ggpubr")

figure <- ggarrange(sig.p +
                    font("x.text", size = 14))
                   # ncol = 2, nrow = 2)

#png(paste0("./3_treatEffect.outs/plot_sig_treatandind_filt_all_0.1_0.1_0.5_cn.png"), width=500, height=500, pointsize=16, res=125)
png(paste0("./3_treatEffect.outs/plot_sig_treatandind_filt_all_0.1_0.1_0_cn_mitofilt.png"), width=500, height=500, pointsize=16, res=125)
print(figure)
dev.off()



## #----DE----
## #peakAll <- unique(sig_resDE$peak)
## #head(peakAll)
## #length(peakAll)

DE  <- sig_resDE %>% dplyr::pull(gene)
DE <- as.character(unique(DE))
head(DE)
length(DE)





#################
#----mat data----
#################
comb <- unique(resDE$comb)
comb

## beta0 <- rep(NA, length(DP))
## beta0

## names(beta0) <- DP
## beta0
## length(beta0)

## d0 <- resDE%>%dplyr::filter(comb=="0-CD4Naive_caffeine", peak%in%DP)
## peakSel <- as.character(d0$peak)
## beta0[peakSel] <- d0$estimate
## beta0
## length(beta0)



#----build matrix for heatmap for all genes----
mat <- map_dfc(comb, function(ii){
##
  beta0 <- rep(NA, length(all.DE))
  names(beta0) <- all.DE
#  d0 <- resDE%>%dplyr::filter(comb==ii, gene%in%all.DE)
  d0 <- resDE%>%dplyr::filter(comb==ii)  
  geneSel <- as.character(d0$gene)
  beta0[geneSel] <- d0$estimate
  beta0
})
head(mat)


mat <- as.matrix(mat)
rownames(mat) <- all.DE
colnames(mat) <- comb
head(mat)
max(mat)
min(mat)

#opfn <- paste0("./3_treatEffect.outs/mat_genes_macs2_0.1_cn.rds")
#write_rds(mat, opfn)




#----build matrix for heatmap for  only significant----
sig_resDE <- sig_resDE.5


mat.sigdgs.qval <- map_dfr(comb, function(ii){
      # res2 <- res.topmotif %>% dplyr::filter(condition==ii)
       res2 <- sig_resDE %>% dplyr::filter(comb==ii)
       res2 <- res2 %>% arrange(p.adjusted)
       res2 <- res2 [1:10, ]
       #b <- res2$beta
       #names(b) <- res2$motif
       #b[topmotif]
    })
head(mat.sigdgs.qval)
nrow(mat.sigdgs.qval)
table(mat.sigdgs.qval$MCls, mat.sigdgs.qval$contrast)
#length(unique(mat.5topmotifs.qval$motif))

topdgs.qval <- unique(mat.sigdgs.qval$gene)
head(topdgs.qval)
length(topdgs.qval)


mat.sigdgs.qval.2 <- map_dfc(comb, function(ii){
        res2 <- resDE %>% dplyr::filter(comb==ii)
            #res2 <- res2 %>% arrange(desc(qval))
            #res2 <- res2[1:5, ]
            #res2
            b <- res2$estimate
            names(b) <- res2$gene
            b[topdgs.qval]
            })
head(mat.sigdgs.qval.2)
nrow(mat.sigdgs.qval.2)

mat <- mat.sigdgs.qval.2



## mat <- map_dfc(comb, function(ii){
## ##
##   beta0 <- rep(NA, length(DE))
##   names(beta0) <- DE
##   d0 <- resDE%>%dplyr::filter(comb==ii, gene%in%DE)
##   geneSel <- as.character(d0$gene)
##   beta0[geneSel] <- d0$estimate
##   beta0
## })
## head(mat)

mat <- as.matrix(mat)
rownames(mat) <- DE
colnames(mat) <- comb
head(mat)
max(mat)
min(mat)

nrow(mat)

ncol(mat)

opfn <- paste0("./3_treatEffect.outs/mat_genes_macs2_0.1_cn.rds")
write_rds(mat, opfn)

mat <- read_rds("./3_treatEffect.outs/mat_genes_macs2_0.1_cn.rds")


#----build matrix for heatmap for top 10 significant only----
sig_resDE <- sig_resDE.5

comb
head(sig_resDE)

mat.10topdgs.qval <- map_dfr(comb, function(ii){
      # res2 <- res.topmotif %>% dplyr::filter(condition==ii)
       res2 <- sig_resDE %>% dplyr::filter(comb==ii)
       res2 <- res2 %>% arrange(p.adjusted)
       res2 <- res2 [1:10, ]
       #b <- res2$beta
       #names(b) <- res2$motif
       #b[topmotif]
    })
head(mat.10topdgs.qval)
nrow(mat.10topdgs.qval)
table(mat.10topdgs.qval$MCls, mat.10topdgs.qval$contrast)
#length(unique(mat.5topmotifs.qval$motif))

topdgs.qval <- unique(mat.10topdgs.qval$gene)
head(topdgs.qval)
length(topdgs.qval)


mat.10topdgs.qval.2 <- map_dfc(comb, function(ii){
        res2 <- resDE %>% dplyr::filter(comb==ii)
            #res2 <- res2 %>% arrange(desc(qval))
            #res2 <- res2[1:5, ]
            #res2
            b <- res2$estimate
            names(b) <- res2$gene
            b[topdgs.qval]
            })
head(mat.10topdgs.qval.2)
nrow(mat.10topdgs.qval.2)




mat <- mat.10topdgs.qval.2



## mat <- map_dfc(comb, function(ii){
## ##
##   beta0 <- rep(NA, length(DE))
##   names(beta0) <- DE
##   d0 <- resDE%>%dplyr::filter(comb==ii, gene%in%DE) %>% dplyr::arrange(p.adjusted) %>% head(n=10)
##   geneSel <- as.character(d0$gene)
##   beta0[geneSel] <- d0$estimate
##   beta0
## })
## head(mat)

mat <- as.matrix(mat)
rownames(mat) <- DE
colnames(mat) <- comb
head(mat)
max(mat)
min(mat)

#opfn <- paste0("./3_treatEffect.outs/mat_genes_macs2_0.1_cn.rds")
#write_rds(mat, opfn)



#----sum filt----
ii <- rowSums(is.na(mat))
ii > 0
head(ii)
length(ii)

#ii <- ii < 56
#head(ii)
#length(ii)
#table (ii)

#EITHER
mat2 <- mat[ii==0,]

#OR
mat2 <- mat[ii<length(comb),]

head(mat2)
nrow(mat2)
max(mat2)
min(mat2)

#mat2[!is.na(mat2)] = 1
mat2[is.na(mat2)] = 0
#head(mat2)
#nrow(mat2)
#max(mat2)
#min(mat2)


#OR
mat[is.na(mat)] = 0

mat2 <- mat

comb.2 <- sort(comb)
comb.2

mat3 <- mat2[, comb.2]
head(mat3)



#################
#----plot----
#################
#----breaks and color----
y <- as.numeric(mat2)
head(y)
min(y)
max(y)

quantile(abs(y),probs=0.99) 

#y0 <- y[abs(y)<1.8] #99% percent quantile(abs(y),probs=0.99)
#head(y0)

y0 <- y[abs(y) quantile(abs(y),probs=0.99)] #99% percent quantile(abs(y),probs=0.99)
head(y0)

quantile(y0,probs=seq(0,1,length.out=98))

mybreaks <- c(min(y),quantile(y0,probs=seq(0,1,length.out=98)),max(y))
names(mybreaks) <- NULL
head(mybreaks)


mycol <- colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100)



#----annotation and colors----
x <- str_split(comb, "_", simplify=T)
x
tmp_column_ct <- data.frame(treatment=x[,1], celltype=x[,2])
rownames(tmp_column_ct) <- comb
head(tmp_column_ct)

x <- str_split(comb.2, "_", simplify=T)
x
tmp_column_treat <- data.frame(treatment=x[,1], celltype=x[,2])
rownames(tmp_column_treat) <- comb.2
head(tmp_column_treat)


## col1 <- c("LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
##    "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
## col2 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
##    "NKcell"="#aa4b56", "Tcell"="#ffaa00")
col1 <- c("caffeine"="red", "nicotine"="tan",
            "vitA"="tan4", "vitD"="seagreen4",
            "vitE"="salmon3", "water"="grey", "zinc"="maroon3")
col2 <- c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
           "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00", "1-TCM"="pink",
           "3-TEM"="blue", "5-CD8Naive"="green", "7-dnT"="black")
tmp_colors <- list(treatment=col1, celltype=col2) #brewer.pal(4,"Set1")

#mycol <- viridisLite::viridis(100)
#mycol <- viridisLite::cividis(100, direction=1)


#----fig----
#cell type
fig1 <- pheatmap(mat2,
   col=mycol, breaks=mybreaks, border_color="NA",
   cluster_rows=T,
   cluster_cols=F,
   annotation_col=tmp_column_ct,
   annotation_colors=tmp_colors, annotation_legend=T,
   show_colnames=T,
   show_rownames=F,
   na_col="white",
   fontsize_col=5,
#   fontsize_row=5
   )


figfn <- "./3_treatEffect.outs/Figure1.1_heatmap.beta.DG.sig.treat_mitofilt.png"
#figfn <- "./3_treatEffect.outs/Figure1.1_heatmap.beta.DG.top10.allvalues.treat.png"
png(figfn, width=1800, height=2100,res=225)
print(fig1)
dev.off() 




#treat
fig1 <- pheatmap(mat3,
   col=mycol, breaks=mybreaks, border_color="NA",
   cluster_rows=T,
   cluster_cols=F,
   annotation_col=tmp_column_treat,
   annotation_colors=tmp_colors, annotation_legend=T,
   show_colnames=T,
   show_rownames=F,
   na_col="white",
   fontsize_col=5,
#   fontsize_row=2
   )

figfn <- "./3_treatEffect.outs/Figure1.1_heatmap.beta.DG.sig.treat_mitofilt.png"
#figfn <- "./3_treatEffect.outs/Figure1.1_heatmap.beta.DG.top10.all.values.treat.png"
png(figfn, width=1800, height=2100,res=225)
print(fig1)
dev.off() 





#######################################
#----plotting as done in motif file----
#######################################
#----annotation as cell type----
x <- str_split(comb, "_", simplify=T)
x
column_ha_ct <- HeatmapAnnotation(
        treatment=x[,1],
        celltype=x[,2],
        col=list(treatment=c("caffeine"="red", "nicotine"="tan",
                             "vitA"="tan4", "vitD"="seagreen4",
                             "vitE"="salmon3", "water"="grey", "zinc"="maroon3"),
                 celltype=c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
                            "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00",
                            "1-TCM"="pink", "3-TEM"="blue",
                            "5-CD8Naive"="green", "7-dnT"="black"))
)


#----annotation as treatment----
x <- str_split(comb.2, "_", simplify=T)
x
column_ha_treat <- HeatmapAnnotation(
        treatment=x[,1],
        celltype=x[,2],
        col=list(treatment=c("caffeine"="red", "nicotine"="tan",
                             "vitA"="tan4", "vitD"="seagreen4",
                             "vitE"="salmon3", "water"="grey", "zinc"="maroon3"),
                 celltype=c("4-Bcell"="#4daf4a", "6-Monocyte"="#984ea3",
                            "2-NKcell"="#aa4b56", "0-CD4Naive"="#ffaa00",
                            "1-TCM"="pink", "3-TEM"="blue",
                            "5-CD8Naive"="green", "7-dnT"="black"))
)


#----breaks and color----
#breaks <- quantile(mat2,probs=seq(0,1,length.out=101))

#col_fun <-

#breaks <- c(seq(0,1,length.out=50), quantile(b2, probs=seq(0, 1, length.out=49)), 38.93)
#breaks
#col_fun <-  colorRamp2(breaks,
#   colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100))


#quantile(mat2, probs=seq(0, 1, length.out=101))

mat2.vector <- as.vector(mat2)
mat2.vector <- sort(mat2.vector)
mat2.vector

breaks <- (seq(-8, 8, length.out=17))
breaks


#col_fun <-  colorRamp2(breaks,
#   colorRampPalette(rev(brewer.pal(n=11, name="RdBu")))(13))

library(ggsci)
library(circlize)

col_fun <- colorRamp2(breaks, colorRampPalette(rev(brewer.pal(n=11, name="RdBu")))(17))


#----figure as celltype----
fig <- Heatmap(mat2,
#               col=col_fun,
               cluster_rows=T, cluster_columns=F,
               top_annotation=column_ha_ct,
               heatmap_legend_param=list(title="fold.enrichment",
                                         title_gp=gpar(fontsize=10),
                                         labels_gp=gpar(fontsize=10)),
               show_row_names=T,
               show_column_names=T,
               row_names_gp=gpar(fontsize=2),
               column_names_gp=gpar(fontsize=5),
               raster_device="png")

#figfn <- "./3_treatEffect.outs/Figure1.1.2_heatmap.beta.DG.sig.cell.png"
figfn <- "./3_treatEffect.outs/Figure1.1.2_heatmap.beta.DG.top10.cell_2.png"
png(figfn, width=1800, height=2100,res=225)
print(fig)
dev.off()



#----figure as treatment---- 
fig <- Heatmap(mat3,
#               col=col_fun,
               cluster_rows=T, cluster_columns=F,
               top_annotation=column_ha_treat,
               heatmap_legend_param=list(title="fold.enrichment",
                                         title_gp=gpar(fontsize=10),
                                         labels_gp=gpar(fontsize=10)),
               show_row_names=T,
               show_column_names=T,
               row_names_gp=gpar(fontsize=3),
               column_names_gp=gpar(fontsize=8),
               raster_device="png")

#figfn <- "./3_treatEffect.outs/Figure1.1.2_heatmap.beta.DG.sig.treat.png"
figfn <- "./3_treatEffect.outs/Figure1.1.2_heatmap.beta.DG.top10.treat_mitofilt.png"
png(figfn, width=3000, height=3500,res=300)
print(fig)
dev.off()

end
#--------




###########################
### correlation heatmap ###
###########################
#Neworder <- c("Monocyte_LPS", "Monocyte_PHA", "Bcell_LPS", "Bcell_PHA",
   "NKcell_LPS", "NKcell_PHA", "Tcell_LPS", "Tcell_PHA",
   "Monocyte_LPS-DEX", "Monocyte_PHA-DEX", "Bcell_LPS-DEX", "Bcell_PHA-DEX",
   "NKcell_LPS-DEX", "NKcell_PHA-DEX", "Tcell_LPS-DEX", "Tcell_PHA-DEX")

Neworder <- c(comb.2)
Neworder

corr <- cor(mat3)[Neworder, Neworder]
head(corr)

mycol <- colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100)
##

###
x <- str_split(colnames(corr), "_", simplify=T)

tmp_column <- data.frame(celltype=x[,2], treatment=x[,1])
rownames(tmp_column) <- colnames(corr)

tmp_colors <- list(celltype=col2, treatment=col1)

fig2 <- pheatmap(corr, col=mycol, scale="none", border_color="NA",
   cluster_rows=F, cluster_cols=F,
   annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend =T,
   show_colnames=T, show_rownames=F,
   fontsize_row=5,
   fontsize_col=5,
   na_col="white")

###
figfn <- "./3_treatEffect.outs/Figure1.2_heatmap.corr.DG.sig.treat.png"
png(figfn, width=600, height=600,res=120)
print(fig2)
dev.off() 































#####################################################################
#### Sarah's data <- bulk RNA seq ###
#####################################################################
## deg.counts <- read.table("/wsu/home/fq/fq93/fq9336/rprdata/mohammedhusain/DEG/after_correcting_script/batch_effect_ethnicity_only_protein_coding/dge_counts.txt")
## head(deg.counts)
## nrow(deg.counts)
## ncol(deg.counts)

#----resDEG----
resDEG <- read_rds("../2_Differential/1_DEG.outs/2.0_DESeq.results.rds") %>% as.data.frame()

head(resDEG)
nrow(resDEG)
max(resDEG$estimate)
min(resDEG$estimate)


#resDEG$Significance <- NULL
#head(resDEG)

resDEG[!is.na(resDEG$p.adjusted) & resDEG$p.adjusted<0.1 & resDEG$estimate > 0.5, "Significance"] <- "Up"
resDEG[!is.na(resDEG$p.adjusted) & resDEG$p.adjusted<0.1 & resDEG$estimate < -0.5, "Significance"] <- "Down"
resDEG[!is.na(resDEG$p.adjusted) & resDEG$p.adjusted>=0.1, "Significance"] <- "Not Significant"
resDEG[is.na(resDEG$p.adjusted),"Significance"] <- "NA"
table(resDEG$Significance)


resDEG[!is.na(resDEG$p.adjusted) & resDEG$p.adjusted<0.1 & resDEG$estimate > 0.25, "Significance.25"] <- "Up"
resDEG[!is.na(resDEG$p.adjusted) & resDEG$p.adjusted<0.1 & resDEG$estimate < -0.25, "Significance.25"] <- "Down"
resDEG[!is.na(resDEG$p.adjusted) & resDEG$p.adjusted>=0.1, "Significance.25"] <- "Not Significant"
resDEG[is.na(resDEG$p.adjusted),"Significance.25"] <- "NA"
table(resDEG$Significance.25)


resDEG[!is.na(resDEG$p.adjusted) & resDEG$p.adjusted<0.1 & resDEG$estimate > 0, "Significance.0"] <- "Up"
resDEG[!is.na(resDEG$p.adjusted) & resDEG$p.adjusted<0.1 & resDEG$estimate < 0, "Significance.0"] <- "Down"
resDEG[!is.na(resDEG$p.adjusted) & resDEG$p.adjusted>=0.1, "Significance.0"] <- "Not Significant"
resDEG[is.na(resDEG$p.adjusted),"Significance.0"] <- "NA"
table(resDEG$Significance.0)


resDEG <- resDEG %>% mutate(comb_2=paste(contrast, Significance, sep="_"))
head(resDEG)
length(unique(resDEG$gene))


resDEG <- resDEG %>% mutate(contrast=gsub("Caffeine", "caffeine", contrast)) %>% mutate(contrast=gsub("VitaminA", "vitA", contrast)) %>% mutate(contrast=gsub("Water", "water", contrast)) %>% mutate(contrast=gsub("Zinc", "zinc", contrast))



#dplyr::filter does not work for vector
#mit.genes <- resDEG$gene %>% dplyr::filter(substr(resDEG$gene, 1, 2)=="MT") 
#the foll works but dont need as doing it other way
#mit.genes <- resDEG$gene[substr(resDEG$gene, 1, 3)=="MT-"] 
#mit.genes

#this select all the rows with "MT" present anywhere in string
#library(data.table)
#mit.resDEG <- resDEG[resDEG$gene %like% "MT", ]
#mit.resDEG$gene

#----no need of this step removed all MT genes----
#resDEG <- resDEG %>% mutate(mito=ifelse(substr(resDEG$gene, 1, 3)=="MT-", 1, 0))
#head(resDEG)
#nrow(resDEG)
#table(resDEG$mito)

#----mitofilt resDEG----
## resDEG <- resDEG %>% dplyr::filter(resDEG$mito==0)
## #head(resDEG)
## nrow(resDEG)
## table(resDEG$Significance)
## length(unique(resDEG$gene))
## table(resDEG$mito)




#----all.DEG----
all.DEG  <- resDEG %>% dplyr::pull(gene)
length(all.DEG)
all.DEG <- as.character(unique(all.DEG))
head(all.DEG)
length(all.DEG)




#----can get gene names from bioMART----
library(biomaRt)

mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
#mart

genes.table <- try(biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description", "gene_biotype", "chromosome_name", "start_position"), mart = mart, useCache = F))

head(genes.table)

#colnames(genes.table)
#table(genes.table$chromosome_name)
#nrow(genes.table)

genes.table.2 <- genes.table%>%dplyr::filter(!chromosome_name %in% c("X", "Y", "MT"), ensembl_gene_id %in% all.DEG) %>% mutate(gene=ensembl_gene_id) %>% dplyr::select(gene, external_gene_name) 
#table(genes.table.2$chromosome_name)
#colnames(genes.table.2)
head(genes.table.2)
nrow(genes.table.2)

#no_genes <- genes.table.2$external_gene_name
#length(no_genes)

#check <- no_genes %in% count_genes
#check


#----annotate gene names----
#head(resDEG)
#nrow(resDEG)

resDEG.join <- left_join(resDEG, genes.table.2, by = "gene") %>% mutate(ensemble=gene) %>% mutate(gene=external_gene_name)
head(resDEG.join)
nrow(resDEG.join)

#table(resDEG.join$contrast)
#table(resDE$contrast)

#already done above
#resDEG.join <- resDEG.join %>% mutate(contrast=gsub("Caffeine", "caffeine", contrast)) %>% mutate(contrast=gsub("VitaminA", "vitA", contrast)) %>% mutate(contrast=gsub("Water", "water", contrast)) %>% mutate(contrast=gsub("Zinc", "zinc", contrast))

resDEG.join <- resDEG.join %>% mutate(contrast.gene=paste0(contrast, ".", gene))

resDEG.join <- resDEG.join %>% mutate(z_score=estimate/stderror)

head(resDEG.join)
nrow(resDEG.join)
table(resDEG.join$contrast)

table(resDEG.join$Significance)
table(resDEG.join$Significance.25)
table(resDEG.join$Significance.0)




#----counting genes for all treatments----
#~~~~caff~~~~
resDEG.join.caff <- resDEG.join %>% dplyr::filter(contrast=="caffeine")
head(resDEG.join.caff)
nrow(resDEG.join.caff)
table(resDEG.join.caff$Significance)
table(resDEG.join.caff$Significance.25)


resDEG.join.caff.sum <- resDEG.join.caff %>% arrange(p.adjusted)
resDEG.join.caff.sum <- resDEG.join.caff.sum %>% distinct(ensemble, .keep_all = TRUE)
#%>% summarise_each(funs(mean))
head(resDEG.join.caff.sum)
nrow(resDEG.join.caff.sum)
table(resDEG.join.caff.sum$Significance)
table(resDEG.join.caff.sum$Significance.25)


resDEG.join.caff.sum <- resDEG.join.caff %>% arrange(p.adjusted)
resDEG.join.caff.sum <- resDEG.join.caff.sum %>% distinct(gene, .keep_all = TRUE)
#%>% summarise_each(funs(mean))
head(resDEG.join.caff.sum)
nrow(resDEG.join.caff.sum)
table(resDEG.join.caff.sum$Significance)
table(resDEG.join.caff.sum$Significance.25)


#~~~~vitA~~~~
resDEG.join.vitA <- resDEG.join %>% dplyr::filter(contrast=="vitA")
head(resDEG.join.vitA)
nrow(resDEG.join.vitA)
table(resDEG.join.vitA$Significance)
table(resDEG.join.vitA$Significance.25)


resDEG.join.vitA.sum <- resDEG.join.vitA %>% arrange(p.adjusted)
resDEG.join.vitA.sum <- resDEG.join.vitA.sum %>% distinct(ensemble, .keep_all = TRUE)
#%>% summarise_each(funs(mean))
head(resDEG.join.vitA.sum)
nrow(resDEG.join.vitA.sum)
table(resDEG.join.vitA.sum$Significance)
table(resDEG.join.vitA.sum$Significance.25)


resDEG.join.vitA.sum <- resDEG.join.vitA %>% arrange(p.adjusted)
resDEG.join.vitA.sum <- resDEG.join.vitA.sum %>% distinct(gene, .keep_all = TRUE)
#%>% summarise_each(funs(mean))
head(resDEG.join.vitA.sum)
nrow(resDEG.join.vitA.sum)
table(resDEG.join.vitA.sum$Significance)
table(resDEG.join.vitA.sum$Significance.25)


#~~~~zinc~~~~
resDEG.join.zinc <- resDEG.join %>% dplyr::filter(contrast=="zinc")
head(resDEG.join.zinc)
nrow(resDEG.join.zinc)
table(resDEG.join.zinc$Significance)
table(resDEG.join.zinc$Significance.25)


resDEG.join.zinc.sum <- resDEG.join.zinc %>% arrange(p.adjusted)
resDEG.join.zinc.sum <- resDEG.join.zinc.sum %>% distinct(ensemble, .keep_all = TRUE)
#%>% summarise_each(funs(mean))
head(resDEG.join.zinc.sum)
nrow(resDEG.join.zinc.sum)
table(resDEG.join.zinc.sum$Significance)
table(resDEG.join.zinc.sum$Significance.25)


resDEG.join.zinc.sum <- resDEG.join.zinc %>% arrange(p.adjusted)
resDEG.join.zinc.sum <- resDEG.join.zinc.sum %>% distinct(gene, .keep_all = TRUE)
#%>% summarise_each(funs(mean))
head(resDEG.join.zinc.sum)
nrow(resDEG.join.zinc.sum)
table(resDEG.join.zinc.sum$Significance)
table(resDEG.join.zinc.sum$Significance.25)











## #----sig_resDEG----
## sig_resDEG <- resDEG.join %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)
## head(sig_resDEG)

## nrow(sig_resDEG)
## max(sig_resDEG$estimate)
## min(sig_resDEG$estimate)
## length(unique(sig_resDEG$gene))

## sigDEG <- table(sig_resDEG$contrast)
## sigDEG
## table(sig_resDEG$Significance)
## table(sig_resDEG$contrast, sig_resDEG$Significance)

## #to check
## nrow(sig_resDEG %>% dplyr::filter(abs(estimate)>5))



## #----checking for duplicated genes----
## nrow(resDE)
## length(unique(resDE$gene))

## head(resDEG)
## nrow(resDEG)
## length(unique(resDEG$gene))

## resDEG.caff <- resDEG %>% dplyr::filter(contrast=="caffeine")
## head(resDEG.caff)

## nrow(resDEG.caff)
## length(unique(resDEG.caff$gene))

## nrow(resDEG.join.caff)
## length(unique(resDEG.join.caff$gene))




## #----for significant only----
## #----table plots resDE and resDEG----
## #----
## head(sig_resDE)
## sig <- table(sig_resDE$MCls, sig_resDE$contrast)
## sig

## sigDE <- table(sig_resDE$contrast)
## sigDE

## library(reshape2)
## sig.matrix <- as.matrix(sigDE)
## sig.matrix
## sig.melt <- melt(sig.matrix)
## sig.melt

## sig.p <- ggplot(sig.melt, aes(x = Var2, y = Var1)) +
##           geom_tile(aes(fill=value), color="black")+
##           scale_fill_gradient(low = "white", high = "white")+
##           geom_text(aes(label = value), color = "black", size = 5) +
##           ggtitle("Significant_0.1_0.1-0.5") +
##           theme(axis.text.x = element_text(angle = 90),
##                           axis.title.x = element_blank(),
##                           axis.title.y = element_blank(),
##                           legend.position="none")

## library("ggpubr")
## figure <- ggarrange(sig.p +
##                     font("x.text", size = 14))
##                    # ncol = 2, nrow = 2)

## png(paste0("./1_DEG.outs/plot_sigDE.png"), width=500, height=500, pointsize=16, res=125)
## print(figure)
## dev.off()


## #----
## head(sig_resDEG)
## sigDEG <- table(sig_resDEG$contrast)
## sigDEG

## #library(reshape2)
## sig.matrix <- as.matrix(sigDEG)
## sig.matrix
## sig.melt <- melt(sig.matrix)
## sig.melt

## sig.p <- ggplot(sig.melt, aes(x = Var2, y = Var1)) +
##           geom_tile(aes(fill=value), color="black")+
##           scale_fill_gradient(low = "white", high = "white")+
##           geom_text(aes(label = value), color = "black", size = 5) +
##           ggtitle("Significant_0.1_0.1-0.5") +
##           theme(axis.text.x = element_text(angle = 90),
##                           axis.title.x = element_blank(),
##                           axis.title.y = element_blank(),
##                           legend.position="none")



## #library("ggpubr")
## figure <- ggarrange(sig.p +
##                     font("x.text", size = 14))
##                    # ncol = 2, nrow = 2)

## png(paste0("./1_DEG.outs/plot_sigDEG.png"), width=500, height=500, pointsize=16, res=125)
## print(figure)
## dev.off()






## #----scatter plots for resDE and resDEG----
## table(sig_resDE$contrast)

## sig_resDE.caff <- sig_resDE %>% dplyr::filter(contrast=="caffeine")

## head(sig_resDE.caff)

## #nrow(sig_resDE.caff)
## #length(sig_resDE.caff$gene)

## table(sig_resDE.caff$MCls)

## sig_resDE.nic <- sig_resDE %>% dplyr::filter(contrast=="nicotine")
## sig_resDE.vitA <- sig_resDE %>% dplyr::filter(contrast=="vitA")
## sig_resDE.vitD <- sig_resDE %>% dplyr::filter(contrast=="vitD")
## sig_resDE.vitE <- sig_resDE %>% dplyr::filter(contrast=="vitE")
## sig_resDE.water <- sig_resDE %>% dplyr::filter(contrast=="water")
## sig_resDE.zinc <- sig_resDE %>% dplyr::filter(contrast=="zinc")
## sig_resDE.caff.vitA.zinc <- sig_resDE %>% dplyr::filter(contrast %in% c("caffeine", "vitA", "zinc"))
## #nrow(sig_resDE.caff.vitA.zinc)

## sig.resDE.caff.genes <- unique(sig_resDE.caff$gene)
## length(sig.resDE.caff.genes)


## #sig_resDE.caff.cd4 <- sig_resDE %>% dplyr::filter(contrast=="caffeine", MCls=="0-CD4Naive")
## #length(sig_resDE.caff.cd4$gene)
## #length(unique(sig_resDE.caff.cd4$gene))


## #----
## table(sig_resDEG$contrast)

## sig_resDEG.caff <- sig_resDEG %>% dplyr::filter(contrast=="Caffeine")

## head(sig_resDEG.caff)

## #length(sig_resDEG.caff$gene)
## sig_resDEG.vitA <- sig_resDEG %>% dplyr::filter(contrast=="VitaminA")
## sig_resDEG.zinc <- sig_resDEG %>% dplyr::filter(contrast=="Zinc")

## sig.resDEG.caff.genes <- unique(sig_resDEG.caff$gene) 
## length(sig.resDEG.caff.genes)
## sig.resDEG.vitA.genes <- unique(sig_resDEG.vitA$gene) 
## length(sig.resDEG.vitA.genes)
## sig.resDEG.zinc.genes <- unique(sig_resDEG.zinc$gene) 
## length(sig.resDEG.zinc.genes)


## #----caff----
## head(sig_resDE.caff)

## MCls <- unique(sig_resDE.caff$MCls)
## MCls

## comb <- MCls

## sig_comb.caff <- map_dfr(comb, function(ii){
##             resDE2 <- sig_resDE.caff%>%dplyr::filter(MCls==ii)
##             DEG <- sig.resDEG.caff.genes
##             x <- resDE2 %>% dplyr::filter(gene%in%DEG)
##             if(nrow(x)==0){
##                 res1 <- NULL
##             }else{
## #                res1 <- x%>%mutate(gene=gene_name)
##                 res1 <- left_join(x,sig_resDEG.caff, by="gene")
##             }
##             #res <- list(DEG=res1)
##             res1
## })

## head(sig_comb.caff)
## nrow(sig_comb.caff)
## length(unique(sig_comb.caff$gene))
## table(sig_comb.caff$comb)


## comb.2 <- unique(sig_comb.caff$MCls)
## comb.2

## for(i in comb.2){
##             res1 <- sig_comb.caff %>% dplyr::filter(MCls==i)
##             png(paste0("./1_DEG.outs/scatter_caff_", i, ".png"))
##             p <- ggplot(res1, aes(x=estimate.y, y=estimate.x, color=MCls)) +
##                     geom_point(show.legend = FALSE)+
##                     xlab("Sarah_data") + ylab("sc_data")+
## #                    xlim(-4, 10) + ylim(-4, 5)+
##                     ggtitle(paste0("logFC_caff_", i))+
##                     theme_bw()+
##                     theme(plot.title = element_text(color = "black", size = 18, face = "bold"),
##                           axis.text = element_text(size = 18),
##                           axis.title = element_text(size = 18))
##             print(p)
##             dev.off()
##         }


## #----vitA----
## head(sig_resDE.vitA)

## MCls <- unique(sig_resDE.vitA$MCls)
## MCls

## comb <- MCls

## sig_comb.vitA <- map_dfr(comb, function(ii){
##             resDE2 <- sig_resDE.vitA%>%dplyr::filter(MCls==ii)
##             DEG <- sig.resDEG.vitA.genes
##             x <- resDE2 %>% dplyr::filter(gene%in%DEG)
##             if(nrow(x)==0){
##                 res1 <- NULL
##             }else{
## #                res1 <- x%>%mutate(gene=gene_name)
##                 res1 <- left_join(x,sig_resDEG.vitA, by="gene")
##             }
##             #res <- list(DEG=res1)
##             res1
## })

## head(sig_comb.vitA)
## nrow(sig_comb.vitA)
## length(unique(sig_comb.vitA$gene))
## table(sig_comb.vitA$comb)


## comb.2 <- unique(sig_comb.vitA$MCls)
## comb.2

## for(i in comb.2){
##             res1 <- sig_comb.vitA %>% dplyr::filter(MCls==i)
##             png(paste0("./1_DEG.outs/scatter_vitA_", i, ".png"))
##             p <- ggplot(res1, aes(x=estimate.y, y=estimate.x, color=MCls)) +
##                     geom_point(show.legend = FALSE)+
##                     xlab("Sarah_data") + ylab("sc_data")+
## #                    xlim(-4, 10) + ylim(-4, 5)+
##                     ggtitle(paste0("logFC_vitA_", i))+
##                     theme_bw()+
##                     theme(plot.title = element_text(color = "black", size = 18, face = "bold"),
##                           axis.text = element_text(size = 18),
##                           axis.title = element_text(size = 18))
##             print(p)
##             dev.off()
##         }



## #----zinc----
## head(sig_resDE.zinc)

## MCls <- unique(sig_resDE.zinc$MCls)
## MCls

## comb <- MCls

## sig_comb.zinc <- map_dfr(comb, function(ii){
##             resDE2 <- sig_resDE.zinc%>%dplyr::filter(MCls==ii)
##             DEG <- sig.resDEG.zinc.genes
##             x <- resDE2 %>% dplyr::filter(gene%in%DEG)
##             if(nrow(x)==0){
##                 res1 <- NULL
##             }else{
## #                res1 <- x%>%mutate(gene=gene_name)
##                 res1 <- left_join(x,sig_resDEG.zinc, by="gene")
##             }
##             #res <- list(DEG=res1)
##             res1
## })

## head(sig_comb.zinc)
## nrow(sig_comb.zinc)
## length(unique(sig_comb.zinc$gene))
## table(sig_comb.zinc$comb)


## comb.2 <- unique(sig_comb.zinc$MCls)
## comb.2

## for(i in comb.2){
##             res1 <- sig_comb.zinc %>% dplyr::filter(MCls==i)
##             png(paste0("./1_DEG.outs/scatter_zinc_", i, ".png"))
##             p <- ggplot(res1, aes(x=estimate.y, y=estimate.x, color=MCls)) +
##                     geom_point(show.legend = FALSE)+
##                     xlab("Sarahdata") + ylab("scdata")+
## #                    xlim(-4, 10) + ylim(-4, 5)+
##                     ggtitle(paste0("logFC_zinc_", i))+
##                     theme_bw()+
##                     theme(plot.title = element_text(color = "black", size = 18, face = "bold"),
##                           axis.text = element_text(size = 18),
##                           axis.title = element_text(size = 18))
##             print(p)
##             dev.off()
##         }


## ## png("./1_DEG.outs/scatter_caff.png")
## ## #p <- ggplot(comb, aes(x=estimate.y, y=estimate.x, color=comb.x)) +
## ## p <- ggplot(aes(x=estimate.y, y=estimate.x)) +
## ##     geom_point()+
## ##     xlab("genes") + ylab("filtered_peaks")+
## ##     #    xlim(-4, 10) + ylim(-4, 10)+
## ##     ggtitle("fold.enrichment")
## ## print(p)
## ## dev.off()








#----for all genes----
#----joining resDE and resDEG----
head(resDE)
nrow(resDE)

resDE.DEG <- left_join(resDE,resDEG.join, by="contrast.gene")
head(resDE.DEG)
nrow(resDE.DEG)

table(resDE.DEG$contrast.x)
table(resDE.DEG$contrast.y)

table(resDE.DEG$Significance.x)
table(resDE.DEG$Significance.y)


## resDE.DEG.merge <- merge(resDE,resDEG.join, by="contrast.gene")


resDE.DEG <- resDE.DEG %>% filter(!is.na(estimate.y))
#%>% filter(!is.na(estimate.x))
head(resDE.DEG)
nrow(resDE.DEG)

table(resDE.DEG$contrast.x)
table(resDE.DEG$contrast.y)

table(resDE.DEG$Significance.x)
table(resDE.DEG$Significance.y)


resDE.DEG <- resDE.DEG %>% mutate(Significance=
             ifelse(Significance.y %in% c("Up", "Down") &  Significance.x %in% c("Up", "Down"), "both",
             ifelse(Significance.y %in% c("Up", "Down") & (Significance.x %in% c("Not Significant", "NA") | is.na(Significance.x)), "Sarah",
             ifelse((Significance.y %in% c("Not Significant", "NA") | is.na(Significance.y)) & Significance.x %in% c("Up", "Down"), "sc", NA))))
head(resDE.DEG)
nrow(resDE.DEG)
table(resDE.DEG$Significance)



resDE.DEG <- resDE.DEG %>% mutate(Significance.25=
             ifelse(Significance.25.y %in% c("Up", "Down") &  Significance.25.x %in% c("Up", "Down"), "both",
             ifelse(Significance.25.y %in% c("Up", "Down") & (Significance.25.x %in% c("Not Significant", "NA") | is.na(Significance.25.x)), "Sarah",
             ifelse((Significance.25.y %in% c("Not Significant", "NA") | is.na(Significance.25.y)) & Significance.25.x %in% c("Up", "Down"), "sc", NA))))
head(resDE.DEG)
nrow(resDE.DEG)
table(resDE.DEG$Significance.25)



resDE.DEG <- resDE.DEG %>% mutate(Significance.0=
             ifelse(Significance.0.y %in% c("Up", "Down") &  Significance.0.x %in% c("Up", "Down"), "both",
             ifelse(Significance.0.y %in% c("Up", "Down") & (Significance.0.x %in% c("Not Significant", "NA") | is.na(Significance.0.x)), "Sarah",
             ifelse((Significance.0.y %in% c("Not Significant", "NA") | is.na(Significance.0.y)) & Significance.0.x %in% c("Up", "Down"), "sc", NA))))
head(resDE.DEG)
nrow(resDE.DEG)
table(resDE.DEG$Significance.0)







table(resDE.DEG$contrast.x)
table(resDE.DEG$contrast.y)
#colnames(resDE.DEG)

#resDE.DEG <- resDE.DEG %>% dplyr::filter(!is.na(Significance))
#head(resDE.DEG)
#nrow(resDE.DEG)
#table(resDE.DEG$Significance)


resDE.DEG.caff <- resDE.DEG %>% dplyr::filter(contrast.x=="caffeine")
#head(resDE.DEG.caff)
#nrow(resDE.DEG.caff)
resDE.DEG.vitA <- resDE.DEG %>% dplyr::filter(contrast.x=="vitA")
resDE.DEG.zinc <- resDE.DEG %>% dplyr::filter(contrast.x=="zinc")


## table(resDE.DEG.caff$Significance.x)
## table(resDE.DEG.caff$Significance.y)

## table(resDE.DEG.vitA$Significance.x)
## table(resDE.DEG.vitA$Significance.y)

## table(resDE.DEG.zinc$Significance.x)
## table(resDE.DEG.zinc$Significance.y)

## table(resDE.DEG.caff$Significance)
## table(resDE.DEG.vitA$Significance)
## table(resDE.DEG.zinc$Significance)




#----checking the outliers----
resDE.DEG.caff.2 <- resDE.DEG.caff %>% dplyr::filter(!Significance %in% c("sc", "Sarah", "both")) %>% dplyr::filter(!estimate.y  %in% c("NA")) %>% dplyr::filter(!is.na(estimate.y))
                                                     #estimate.y > 10)


head(resDE.DEG.caff.2)

min(resDE.DEG.caff.2$z_score.y)
max(resDE.DEG.caff.2$z_score.y)


#----conditions----
logFC <- 0.5



#----caff----
MCls <- unique(resDE.DEG.caff$MCls)
MCls

comb.2 <- MCls

for(i in comb.2){
            res1 <- resDE.DEG.caff %>% dplyr::filter(MCls==i)
            png(paste0("./1_DEG.outs/logFC_", logFC, "_all_scatter_z_caff_", i, ".png"))
#            p <- ggplot(res1, aes(x=estimate.y, y=estimate.x, color=Significance)) +
            p <- ggplot(res1, aes(x=z_score.y, y=z_score.x, color=Significance.25)) +
                    geom_point(show.legend = TRUE)+
                    xlab("Sarah_data") + ylab("sc_data")+
#                    xlim(-4, 10) + ylim(-4, 5)+
                    ggtitle(paste0("logFC_", logFC, "_all_zscore_caff_", i))+
                    theme_bw()+
                    theme(plot.title = element_text(color = "black", size = 18, face = "bold"),
                          axis.text = element_text(size = 18),
                          axis.title = element_text(size = 18))
            print(p)
            dev.off()
        }


#----vitA----
MCls <- unique(resDE.DEG.vitA$MCls)
MCls

comb.2 <- MCls

for(i in comb.2){
            res1 <- resDE.DEG.vitA %>% dplyr::filter(MCls==i)
            png(paste0("./1_DEG.outs/logFC_", logFC, "_all_scatter_z_vitA_", i, ".png"))
#            p <- ggplot(res1, aes(x=estimate.y, y=estimate.x, color=Significance)) +
            p <- ggplot(res1, aes(x=z_score.y, y=z_score.x, color=Significance.25)) +
                    geom_point(show.legend = TRUE)+
                    xlab("Sarah_data") + ylab("sc_data")+
#                    xlim(-4, 10) + ylim(-4, 5)+
                    ggtitle(paste0("logFC_", logFC, "_all_zscore_vitA_", i))+
                    theme_bw()+
                    theme(plot.title = element_text(color = "black", size = 18, face = "bold"),
                          axis.text = element_text(size = 18),
                          axis.title = element_text(size = 18))
            print(p)
            dev.off()
        }



#----zinc----
MCls <- unique(resDE.DEG.zinc$MCls)
MCls

comb.2 <- MCls

for(i in comb.2){
            res1 <- resDE.DEG.zinc %>% dplyr::filter(MCls==i)
            png(paste0("./1_DEG.outs/logFC_", logFC, "_all_scatter_z_zinc_", i, ".png"))
#            p <- ggplot(res1, aes(x=estimate.y, y=estimate.x, color=Significance)) +
            p <- ggplot(res1, aes(x=z_score.y, y=z_score.x, color=Significance.25)) +
                    geom_point(show.legend = TRUE)+
                    xlab("Sarahdata") + ylab("scdata")+
#                    xlim(-4, 10) + ylim(-4, 5)+
                    ggtitle(paste0("logFC_", logFC, "_all_zscore_zinc_", i))+
                    theme_bw()+
                    theme(plot.title = element_text(color = "black", size = 18, face = "bold"),
                          axis.text = element_text(size = 18),
                          axis.title = element_text(size = 18))
            print(p)
            dev.off()
        }




#----spearman_corr----
scorr_caff <- c()
scorr_vitA <- c()
scorr_zinc <- c()

#scorr_caff
#head(resDE.DEG.caff)

MCls <- unique(resDE.DEG.caff$MCls)
MCls
comb.2 <- MCls


conditions <- c("both", "sc", "Sarah")
include <- "union"

conditions <- c("both")
include <- "both"

logFC <- 0.5
 




#cor.test(res1$z_score.y, res1$z_score.x, method = "spearman")
#as.numeric(cor.test(res1$z_score.y, res1$z_score.x, method = "spearman")$estimate)


for(i in comb.2){
            res1 <- resDE.DEG.caff %>% dplyr::filter(MCls==i)                                 %>% dplyr::filter(Significance %in% conditions)  
            #scorr <- cor(res1$z_score.y, res1$z_score.x, method = "spearman")
            scorr <- as.numeric(cor.test(res1$z_score.y, res1$z_score.x, method = "spearman")$estimate)
            scorr_caff <- c(scorr_caff, scorr)
        }


for(i in comb.2){
            res1 <- resDE.DEG.vitA %>% dplyr::filter(MCls==i)                                 %>% dplyr::filter(Significance %in% conditions)
            #scorr <- cor(res1$z_score.y, res1$z_score.x, method = "spearman")
            scorr <- as.numeric(cor.test(res1$z_score.y, res1$z_score.x, method = "spearman")$estimate)
            scorr_vitA <- c(scorr_vitA, scorr)
        }


for(i in comb.2){
            res1 <- resDE.DEG.zinc %>% dplyr::filter(MCls==i)                                 %>% dplyr::filter(Significance %in% conditions)
            #scorr <- cor(res1$z_score.y, res1$z_score.x, method = "spearman")
            scorr <- as.numeric(cor.test(res1$z_score.y, res1$z_score.x, method = "spearman")$estimate)
            scorr_zinc <- c(scorr_zinc, scorr)
        }




#scorr_caff

scorr_data <- data.frame(MCls,
                         scorr_caff,
                         scorr_vitA)
#                         scorr_zinc)
                            



is.num <- sapply(scorr_data, is.numeric)
scorr_data[is.num] <- lapply(scorr_data[is.num], round, 2)

head(scorr_data)


library(reshape2)
#sig.matrix <- as.matrix(scorr_data)
#sig.matrix

sig.melt <- melt(scorr_data)
sig.melt

sig.p <- ggplot(sig.melt, aes(x = variable, y = MCls)) +
          geom_tile(aes(fill=value), color="black")+
          scale_fill_gradient(low = "white", high = "white")+
          geom_text(aes(label = value), color = "black", size = 5) +
          ggtitle(paste0("logFC_", logFC, "_", include, "_z_score_corr_Sarah_vs_sc")) +
          theme(axis.text.x = element_text(angle = 90),
                          axis.title.x = element_blank(),
                          axis.title.y = element_blank(),
                          legend.position="none")

library("ggpubr")

figure <- ggarrange(sig.p +
                    font("x.text", size = 14))
                   # ncol = 2, nrow = 2)

png(paste0("./1_DEG.outs/logFC_", logFC, "_", include, "_z_score_corr.png"), width=500, height=500, pointsize=16, res=125)
print(figure)
dev.off()




#----heatmap----
mat <- scorr_data
rownames(mat) <- mat$MCls
mat <- mat %>% dplyr::select(!MCls)
mat

mat2 <- as.matrix(mat)
mat2 <- mat2[nrow(mat2):1, ] 
mat2


fig1 <- pheatmap(mat2, border_color="NA",
                    cluster_rows=F, cluster_cols=F,
#                    annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend=T,
                    show_colnames=T,
                    show_rownames=T,
                    na_col="white",
                    fontsize_col=20,
                    fontsize_row=10)


figfn <- "./1_DEG.outs/z_score_corr_heatmap.png"
png(figfn, width=1800, height=2100,res=225)
print(fig1)
dev.off() 
















#####################################################################
#### GxE browser data <- bulk RNA seq ###
#####################################################################
deg.counts <- read.table("/wsu/home/fq/fq93/fq9336/rprdata/mohammedhusain/DEG/after_correcting_script/batch_effect_ethnicity_only_protein_coding/dge_counts.txt")
head(deg.counts)
nrow(deg.counts)
ncol(deg.counts)

#----resGXE----
#resGXE <- read_rds("../2_Differential/1_GXE.outs/2.0_DESeq.results.rds") %>% as.data.frame()


#----caff----
## resGXE.caff.1 <- read.table("../2_Differential/DEG_results/treatments/DP1_DEG_stats_T13C1.txt", header=T)
## resGXE.caff.2 <- read.table("../2_Differential/DEG_results/treatments/DP4_DEG_stats_T13C1.txt", header=T)
## resGXE.caff.3 <- read.table("../2_Differential/DEG_results/treatments/DP6_DEG_stats_T13C1.txt", header=T)
## resGXE.caff.4 <- read.table("../2_Differential/DEG_results/treatments/DP9_DEG_stats_T13C1.txt", header=T)
## resGXE.caff.5 <- read.table("../2_Differential/DEG_results/treatments/DP11_DEG_stats_T13C1.txt", header=T)

## head(resGXE.caff.1)
## head(resGXE.caff.2)
## head(resGXE.caff.3)
## head(resGXE.caff.4)
## head(resGXE.caff.5)



## resGXE.caff <- rbind(resGXE.caff.1, resGXE.caff.2, resGXE.caff.3, resGXE.caff.4, resGXE.caff.5)
## head(resGXE.caff)

## resGXE.caff <- bind_rows(resGXE.caff.1, resGXE.caff.2, resGXE.caff.3, resGXE.caff.4, resGXE.caff.5) %>% group_by(t.id) %>% summarise_each(funs(mean))
## head(resGXE.caff)

## resGXE.caff <- resGXE.caff %>% group_by(t.id) %>% mutate(logFC_avg=mean(logFC))
## #%>% summarise_each(funs(mean))
## head(resGXE.caff)



resGXE.caff <- read.table("../2_Differential/DEG_results/treatments/DP6_DEG_stats_T13C1.txt", header=T)

head(resGXE.caff)

resGXE.caff[!is.na(resGXE.caff$padj) & resGXE.caff$padj<0.1 & resGXE.caff$logFC > 0.5, "Significance"] <- "Up"
resGXE.caff[!is.na(resGXE.caff$padj) & resGXE.caff$padj<0.1 & resGXE.caff$logFC < -0.5, "Significance"] <- "Down"
resGXE.caff[!is.na(resGXE.caff$padj) & resGXE.caff$padj>=0.1, "Significance"] <- "Not Significant"
resGXE.caff[is.na(resGXE.caff$padj),"Significance"] <- "NA"
head(resGXE.caff)
nrow(resGXE.caff)
table(resGXE.caff$Significance)


resGXE.caff[!is.na(resGXE.caff$padj) & resGXE.caff$padj<0.1 & resGXE.caff$logFC > 0.25, "Significance.25"] <- "Up"
resGXE.caff[!is.na(resGXE.caff$padj) & resGXE.caff$padj<0.1 & resGXE.caff$logFC < -0.25, "Significance.25"] <- "Down"
resGXE.caff[!is.na(resGXE.caff$padj) & resGXE.caff$padj>=0.1, "Significance.25"] <- "Not Significant"
resGXE.caff[is.na(resGXE.caff$padj),"Significance.25"] <- "NA"
head(resGXE.caff)
nrow(resGXE.caff)
table(resGXE.caff$Significance.25)


resGXE.caff[!is.na(resGXE.caff$padj) & resGXE.caff$padj<0.1 & resGXE.caff$logFC > 0, "Significance.0"] <- "Up"
resGXE.caff[!is.na(resGXE.caff$padj) & resGXE.caff$padj<0.1 & resGXE.caff$logFC < 0, "Significance.0"] <- "Down"
resGXE.caff[!is.na(resGXE.caff$padj) & resGXE.caff$padj>=0.1, "Significance.0"] <- "Not Significant"
resGXE.caff[is.na(resGXE.caff$padj),"Significance.0"] <- "NA"

head(resGXE.caff)
nrow(resGXE.caff)
table(resGXE.caff$Significance.0)


## length(unique(resGXE.caff$g.id))
## length(unique(resGXE.caff$ensg))

## table(resGXE.caff$ensg)
## resGXE.caff %>% dplyr::filter(ensg=="ENSG00000273492")

write.csv(resGXE.caff,"./1_GXE.outs/resGXE_caff.csv", row.names = FALSE)

#----summarize----
## head(resDE)
## head(all.DE)
## length(all.DE)

## resGXE.caff.filt <- resGXE.caff %>% dplyr::filter(g.id %in% all.DE)
## head(resGXE.caff.filt)
## nrow(resGXE.caff.filt)
## length(unique(resGXE.caff.filt$g.id))

## genes_summary <- as.data.frame(table(resGXE.caff$g.id))
## head(genes_summary)
## nrow(genes_summary)
## #genes_summary <- resGXE.caff %>% group_by(g.id) %>% summarise(sum = sum(g.id))

## genes_summary <- genes_summary %>% dplyr::filter(Freq==1)
## head(genes_summary)
## nrow(genes_summary)





#----taking mean and checking----
## resGXE.caff.2 <- resGXE.caff %>% group_by(ensg) %>% summarise_each(funs(mean)) 
## resGXE.caff.2 <- resGXE.caff.2 %>% as.data.frame()
## head(resGXE.caff.2)

## resGXE.caff.2$Significance.25 <- NULL

## resGXE.caff.2[!is.na(resGXE.caff.2$padj) & resGXE.caff.2$padj<0.1 & resGXE.caff.2$logFC > 0.25, "Significance.25"] <- "Up"
## resGXE.caff.2[!is.na(resGXE.caff.2$padj) & resGXE.caff.2$padj<0.1 & resGXE.caff.2$logFC < -0.25, "Significance.25"] <- "Down"
## resGXE.caff.2[!is.na(resGXE.caff.2$padj) & resGXE.caff.2$padj>=0.1, "Significance.25"] <- "Not Significant"
## resGXE.caff.2[is.na(resGXE.caff.2$padj),"Significance.25"] <- "NA"
## head(resGXE.caff.2)
## nrow(resGXE.caff.2)
## table(resGXE.caff.2$Significance.25)

#----validating with publication----
head(resGXE.caff)
nrow(resGXE.caff)




#----GXE sum----
resGXE.caff.sum <- resGXE.caff %>% arrange(padj)
resGXE.caff.sum <- resGXE.caff.sum %>% distinct(ensg, .keep_all = TRUE)
#%>% summarise_each(funs(mean))
head(resGXE.caff.sum)
nrow(resGXE.caff.sum)
table(resGXE.caff.sum$Significance)
table(resGXE.caff.sum$Significance.25)
table(resGXE.caff.sum$Significance.0)


resGXE.caff.sum <- resGXE.caff %>% arrange(padj)
resGXE.caff.sum <- resGXE.caff.sum %>% distinct(g.id, .keep_all = TRUE)
#%>% summarise_each(funs(mean))
head(resGXE.caff.sum)
nrow(resGXE.caff.sum)
table(resGXE.caff.sum$Significance)
table(resGXE.caff.sum$Significance.25)
table(resGXE.caff.sum$Significance.0)


table(resGXE.caff$Significance.25)
nrow(resGXE.caff %>% dplyr::filter(is.na(Significance.25)))


resGXE.caff.sum <- resGXE.caff %>% mutate(sig=ifelse(Significance.25 %in% c("Up", "Down"), "sig", "nonsig"))
#head(resGXE.caff.sum)
#table(resGXE.caff.sum$sig)
resGXE.caff.sum <- resGXE.caff.sum %>% arrange(desc(logFC))
resGXE.caff.sum <- resGXE.caff.sum %>% arrange(padj)
resGXE.caff.sum <- resGXE.caff.sum %>% arrange(desc(sig))
resGXE.caff.sum <- resGXE.caff.sum %>% distinct(ensg, .keep_all = TRUE)
head(resGXE.caff.sum)
nrow(resGXE.caff.sum)
table(resGXE.caff.sum$sig)
table(resGXE.caff.sum$Significance)
table(resGXE.caff.sum$Significance.25)
table(resGXE.caff.sum$Significance.0)




## newdf <- resGXE.caff %>% group_by(g.id) %>% summarise(mean = mean(logFC))
## head(newdf)

#resGXE.caff.sum$Significance <- NULL

## resGXE.caff.sum[!is.na(resGXE.caff.sum$padj) & resGXE.caff.sum$padj<0.1 & resGXE.caff.sum$logFC > 0.5, "Significance"] <- "Up"
## resGXE.caff.sum[!is.na(resGXE.caff.sum$padj) & resGXE.caff.sum$padj<0.1 & resGXE.caff.sum$logFC < -0.5, "Significance"] <- "Down"
## resGXE.caff.sum[!is.na(resGXE.caff.sum$padj) & resGXE.caff.sum$padj>=0.1, "Significance"] <- "Not Significant"
## resGXE.caff.sum[is.na(resGXE.caff.sum$padj),"Significance"] <- "NA"
## head(resGXE.caff.sum)
## nrow(resGXE.caff.sum)
## table(resGXE.caff.sum$Significance)


resGXE.caff.sum <- resGXE.caff.sum %>% mutate(contrast.gene=paste0("caffeine.", g.id))
head(resGXE.caff.sum)
nrow(resGXE.caff.sum)
#table(resGXE.caff.sum$Significance)



#----joining resDE and resGXE.caff----
#resDE <- resDE %>% mutate(contrast.gene=paste0(contrast, ".", gene))
#table(resDE$contrast)
#head(resDE)
#nrow(resDE)

resDE.caff <- resDE %>% dplyr::filter(contrast=="caffeine")
nrow(resDE.caff)

resDE.GXE.caff.sum <- left_join(resDE.caff,resGXE.caff.sum, by="contrast.gene")
head(resDE.GXE.caff.sum)
nrow(resDE.GXE.caff.sum)
table(resDE.GXE.caff.sum$Significance.x)
table(resDE.GXE.caff.sum$Significance.y)


resDE.GXE.caff.sum <- resDE.GXE.caff.sum %>% filter(!is.na(logFC))
#%>% filter(!is.na(estimate.x))
head(resDE.GXE.caff.sum)
nrow(resDE.GXE.caff.sum)
table(resDE.GXE.caff.sum$Significance.x)
table(resDE.GXE.caff.sum$Significance.y)




resDE.GXE.caff.sum <- resDE.GXE.caff.sum %>% mutate(Significance=
             ifelse(Significance.y %in% c("Up", "Down") &  Significance.x %in% c("Up", "Down"), "both",
             ifelse(Significance.y %in% c("Up", "Down") & (Significance.x %in% c("Not Significant", "NA") | is.na(Significance.x)), "GXE",
             ifelse((Significance.y %in% c("Not Significant", "NA") | is.na(Significance.y)) & Significance.x %in% c("Up", "Down"), "sc", "NA"))))
head(resDE.GXE.caff.sum)
nrow(resDE.GXE.caff.sum)
table(resDE.GXE.caff.sum$Significance)



resDE.GXE.caff.sum <- resDE.GXE.caff.sum %>% mutate(Significance.25=
             ifelse(Significance.25.y %in% c("Up", "Down") &  Significance.25.x %in% c("Up", "Down"), "both",
             ifelse(Significance.25.y %in% c("Up", "Down") & (Significance.25.x %in% c("Not Significant", "NA") | is.na(Significance.25.x)), "GXE",
             ifelse((Significance.25.y %in% c("Not Significant", "NA") | is.na(Significance.25.y)) & Significance.25.x %in% c("Up", "Down"), "sc", "NA"))))
head(resDE.GXE.caff.sum)
nrow(resDE.GXE.caff.sum)
table(resDE.GXE.caff.sum$Significance.25)



resDE.GXE.caff.sum <- resDE.GXE.caff.sum %>% mutate(Significance.0=
             ifelse(Significance.0.y %in% c("Up", "Down") &  Significance.0.x %in% c("Up", "Down"), "both",
             ifelse(Significance.0.y %in% c("Up", "Down") & (Significance.0.x %in% c("Not Significant", "NA") | is.na(Significance.0.x)), "GXE",
             ifelse((Significance.0.y %in% c("Not Significant", "NA") | is.na(Significance.0.y)) & Significance.0.x %in% c("Up", "Down"), "sc", "NA"))))

head(resDE.GXE.caff.sum)
colnames(resDE.GXE.caff.sum)
nrow(resDE.GXE.caff.sum)
table(resDE.GXE.caff.sum$Significance.0)











#resDE.GXE.caff.sum <- resDE.GXE.caff.sum %>% dplyr::filter(!is.na(Significance))
#head(resDE.GXE.caff.sum)
#nrow(resDE.GXE.caff.sum)
#table(resDE.GXE.caff.sum$Significance)




#----caff plot----
logFC <- 0.5

MCls <- unique(resDE.GXE.caff.sum$MCls)
MCls
comb.2 <- MCls

for(i in comb.2){
            res1 <- resDE.GXE.caff.sum %>% dplyr::filter(MCls==i)
            png(paste0("./1_GXE.outs/logFC_", logFC, "_all_zvslogfc_scatter_caff_", i, ".png"))
            p <- ggplot(res1, aes(x=logFC, y=z_score, color=Significance)) +
                    scale_colour_manual(values = c("#F8766D", "#00BA38", "darkgrey", "#619CFF")) +
                    geom_point(show.legend = TRUE)+
                    xlab("GXE_data") + ylab("sc_data")+
#                    xlim(-4, 10) + ylim(-4, 5)+
                    ggtitle(paste0("logFC_", logFC, "_all_zvslogFC_caff_", i))+
                    theme_bw()+
                    theme(plot.title = element_text(color = "black", size = 18, face = "bold"),
                          axis.text = element_text(size = 18),
                          axis.title = element_text(size = 18))
            print(p)
            dev.off()
        }












## #----trying plotting in loop for all----

## #up <- resDE.GXE.
## treats <- c("caff", "nic")

## for(i in treats){
##     res1 <- paste0("resDE.GXE.", i, ".sum")
##     print(res1)
##     }





## logFC <- 0

## MCls <- unique(resDE.GXE.caff.sum$MCls)
## MCls
## comb.2 <- MCls

## for(i in comb.2){
##             res1 <- resDE.GXE.caff.sum %>% dplyr::filter(MCls==i)
##             png(paste0("./1_GXE.outs/logFC_", logFC, "_all_zvslogfc_scatter_caff_", i, ".png"))
##             p <- ggplot(res1, aes(x=logFC, y=z_score, color=Significance.0)) +
##                     scale_colour_manual(values = c("#F8766D", "#00BA38", "darkgrey", "#619CFF")) +
##                     geom_point(show.legend = TRUE)+
##                     xlab("GXE_data") + ylab("sc_data")+
## #                    xlim(-4, 10) + ylim(-4, 5)+
##                     ggtitle(paste0("logFC_", logFC, "_all_zvslogFC_caff_", i))+
##                     theme_bw()+
##                     theme(plot.title = element_text(color = "black", size = 18, face = "bold"),
##                           axis.text = element_text(size = 18),
##                           axis.title = element_text(size = 18))
##             print(p)
##             dev.off()
##         }














#----nic----
resGXE.nic <- read.table("../2_Differential/DEG_results/treatments/DP6_DEG_stats_T14C1.txt", header=T)
head(resGXE.nic)

resGXE.nic[!is.na(resGXE.nic$padj) & resGXE.nic$padj<0.1 & resGXE.nic$logFC > 0.5, "Significance"] <- "Up"
resGXE.nic[!is.na(resGXE.nic$padj) & resGXE.nic$padj<0.1 & resGXE.nic$logFC < -0.5, "Significance"] <- "Down"
resGXE.nic[!is.na(resGXE.nic$padj) & resGXE.nic$padj>=0.1, "Significance"] <- "Not Significant"
resGXE.nic[is.na(resGXE.nic$padj),"Significance"] <- "NA"
head(resGXE.nic)
nrow(resGXE.nic)
table(resGXE.nic$Significance)


resGXE.nic[!is.na(resGXE.nic$padj) & resGXE.nic$padj<0.1 & resGXE.nic$logFC > 0.25, "Significance.25"] <- "Up"
resGXE.nic[!is.na(resGXE.nic$padj) & resGXE.nic$padj<0.1 & resGXE.nic$logFC < -0.25, "Significance.25"] <- "Down"
resGXE.nic[!is.na(resGXE.nic$padj) & resGXE.nic$padj>=0.1, "Significance.25"] <- "Not Significant"
resGXE.nic[is.na(resGXE.nic$padj),"Significance.25"] <- "NA"
head(resGXE.nic)
nrow(resGXE.nic)
table(resGXE.nic$Significance.25)


resGXE.nic[!is.na(resGXE.nic$padj) & resGXE.nic$padj<0.1 & resGXE.nic$logFC > 0, "Significance.0"] <- "Up"
resGXE.nic[!is.na(resGXE.nic$padj) & resGXE.nic$padj<0.1 & resGXE.nic$logFC < 0, "Significance.0"] <- "Down"
resGXE.nic[!is.na(resGXE.nic$padj) & resGXE.nic$padj>=0.1, "Significance.0"] <- "Not Significant"
resGXE.nic[is.na(resGXE.nic$padj),"Significance.0"] <- "NA"
head(resGXE.nic)
nrow(resGXE.nic)
table(resGXE.nic$Significance.0)


## head(resDE)
## head(all.DE)
## length(all.DE)

## resGXE.nic.filt <- resGXE.nic %>% dplyr::filter(g.id %in% all.DE)
## head(resGXE.nic.filt)
## nrow(resGXE.nic.filt)
## length(unique(resGXE.nic.filt$g.id))

## genes_summary <- as.data.frame(table(resGXE.nic$g.id))
## head(genes_summary)
## nrow(genes_summary)
## #genes_summary <- resGXE.nic %>% group_by(g.id) %>% summarise(sum = sum(g.id))

## genes_summary <- genes_summary %>% dplyr::filter(Freq==1)
## head(genes_summary)
## nrow(genes_summary)

#resGXE.nic.sum <- resGXE.nic %>% group_by(g.id) %>% summarise_each(funs(mean))

resGXE.nic.sum <- resGXE.nic %>% arrange(padj)
resGXE.nic.sum <- resGXE.nic.sum %>% distinct(ensg, .keep_all = TRUE)
head(resGXE.nic.sum)
nrow(resGXE.nic.sum)
table(resGXE.nic.sum$Significance)
table(resGXE.nic.sum$Significance.25)
table(resGXE.nic.sum$Significance.0)


resGXE.nic.sum <- resGXE.nic %>% arrange(padj)
resGXE.nic.sum <- resGXE.nic.sum %>% distinct(g.id, .keep_all = TRUE)
head(resGXE.nic.sum)
nrow(resGXE.nic.sum)
table(resGXE.nic.sum$Significance)
table(resGXE.nic.sum$Significance.25)
table(resGXE.nic.sum$Significance.0)


resGXE.nic.sum <- resGXE.nic %>% mutate(sig=ifelse(Significance.25 %in% c("Up", "Down"), "sig", "nonsig"))
#head(resGXE.nic.sum)
#table(resGXE.nic.sum$sig)
resGXE.nic.sum <- resGXE.nic.sum %>% arrange(desc(logFC))
resGXE.nic.sum <- resGXE.nic.sum %>% arrange(padj)
resGXE.nic.sum <- resGXE.nic.sum %>% arrange(desc(sig))
resGXE.nic.sum <- resGXE.nic.sum %>% distinct(ensg, .keep_all = TRUE)
head(resGXE.nic.sum)
nrow(resGXE.nic.sum)
table(resGXE.nic.sum$sig)
table(resGXE.nic.sum$Significance)
table(resGXE.nic.sum$Significance.25)
table(resGXE.nic.sum$Significance.0)


## newdf <- resGXE.nic %>% group_by(g.id) %>% summarise(mean = mean(logFC))
## head(newdf)

## resGXE.nic.sum$Significance <- NULL

## resGXE.nic.sum[!is.na(resGXE.nic.sum$padj) & resGXE.nic.sum$padj<0.1 & resGXE.nic.sum$logFC > 0.5, "Significance"] <- "Up"
## resGXE.nic.sum[!is.na(resGXE.nic.sum$padj) & resGXE.nic.sum$padj<0.1 & resGXE.nic.sum$logFC < -0.5, "Significance"] <- "Down"
## resGXE.nic.sum[!is.na(resGXE.nic.sum$padj) & resGXE.nic.sum$padj>=0.1, "Significance"] <- "Not Significant"
## resGXE.nic.sum[is.na(resGXE.nic.sum$padj),"Significance"] <- "NA"
## head(resGXE.nic.sum)
## nrow(resGXE.nic.sum)
## table(resGXE.nic.sum$Significance)


resGXE.nic.sum <- resGXE.nic.sum %>% mutate(contrast.gene=paste0("nicotine.", g.id))
head(resGXE.nic.sum)
nrow(resGXE.nic.sum)
table(resGXE.nic.sum$Significance)



#----joining resDE and resGXE.nic----
#resDE <- resDE %>% mutate(contrast.gene=paste0(contrast, ".", gene))
table(resDE$contrast)

resDE.nic <- resDE %>% dplyr::filter(contrast=="nicotine")
head(resDE.nic)
nrow(resDE.nic)
#nrow(sig_resDE.nic)
#table(sig_resDE.nic$contrast)


resDE.GXE.nic.sum <- left_join(resDE.nic,resGXE.nic.sum, by="contrast.gene")
head(resDE.GXE.nic.sum)
nrow(resDE.GXE.nic.sum)
table(resDE.GXE.nic.sum$Significance.x)
table(resDE.GXE.nic.sum$Significance.y)


resDE.GXE.nic.sum <- resDE.GXE.nic.sum %>% filter(!is.na(logFC))
#%>% filter(!is.na(estimate.x))
head(resDE.GXE.nic.sum)
nrow(resDE.GXE.nic.sum)
table(resDE.GXE.nic.sum$Significance.x)
table(resDE.GXE.nic.sum$Significance.y)


resDE.GXE.nic.sum <- resDE.GXE.nic.sum %>% mutate(Significance=
             ifelse(Significance.y %in% c("Up", "Down") &  Significance.x %in% c("Up", "Down"), "both",
             ifelse(Significance.y %in% c("Up", "Down") & (Significance.x %in% c("Not Significant", "NA") | is.na(Significance.x)), "GXE",
             ifelse((Significance.y %in% c("Not Significant", "NA") | is.na(Significance.y)) & Significance.x %in% c("Up", "Down"), "sc", "NA"))))
head(resDE.GXE.nic.sum)
nrow(resDE.GXE.nic.sum)
table(resDE.GXE.nic.sum$Significance)


resDE.GXE.nic.sum <- resDE.GXE.nic.sum %>% mutate(Significance.25=
             ifelse(Significance.25.y %in% c("Up", "Down") &  Significance.25.x %in% c("Up", "Down"), "both",
             ifelse(Significance.25.y %in% c("Up", "Down") & (Significance.25.x %in% c("Not Significant", "NA") | is.na(Significance.25.x)), "GXE",
             ifelse((Significance.25.y %in% c("Not Significant", "NA") | is.na(Significance.25.y)) & Significance.25.x %in% c("Up", "Down"), "sc", "NA"))))
head(resDE.GXE.nic.sum)
nrow(resDE.GXE.nic.sum)
table(resDE.GXE.nic.sum$Significance.25)



resDE.GXE.nic.sum <- resDE.GXE.nic.sum %>% mutate(Significance.0=
             ifelse(Significance.0.y %in% c("Up", "Down") &  Significance.0.x %in% c("Up", "Down"), "both",
             ifelse(Significance.0.y %in% c("Up", "Down") & (Significance.0.x %in% c("Not Significant", "NA") | is.na(Significance.0.x)), "GXE",
             ifelse((Significance.0.y %in% c("Not Significant", "NA") | is.na(Significance.0.y)) & Significance.0.x %in% c("Up", "Down"), "sc", "NA"))))
head(resDE.GXE.nic.sum)
nrow(resDE.GXE.nic.sum)
table(resDE.GXE.nic.sum$Significance.0)



## resDE.GXE.nic.sum <- resDE.GXE.nic.sum %>% dplyr::filter(!is.na(Significance))
## head(resDE.GXE.nic.sum)
## nrow(resDE.GXE.nic.sum)
## table(resDE.GXE.nic.sum$Significance)


## resDE.GXE.caff.sum
## list <- c(resDE.GXE.caff.sum, resDE.GXE.nic.sum)
## for (i in list){
##     print(i)
##     res.2 <- i
##     print(head(res.2))
##     }

#----nic plot----
logFC <- 0

MCls <- unique(resDE.GXE.nic.sum$MCls)
MCls
comb.2 <- MCls

for(i in comb.2){
            res1 <- resDE.GXE.nic.sum %>% dplyr::filter(MCls==i)
            png(paste0("./1_GXE.outs/logFC", logFC, "_all_scatter_nic_", i, ".png"))
            p <- ggplot(res1, aes(x=logFC, y=estimate, color=Significance.0)) +
                scale_colour_manual(values = c("#F8766D", "#00BA38", "darkgrey", "#619CFF")) +
                    geom_point(show.legend = TRUE)+
                    xlab("GXE_data") + ylab("sc_data")+
#                    xlim(-4, 10) + ylim(-4, 5)+
                    ggtitle(paste0("logFC_", logFC, "_all_logFC_nic_", i))+
                    theme_bw()+
                    theme(plot.title = element_text(color = "black", size = 18, face = "bold"),
                          axis.text = element_text(size = 18),
                          axis.title = element_text(size = 18))
            print(p)
            dev.off()
        }
































#----vitA----
## resGXE.vitA.1 <- read.table("../2_Differential/DEG_results/treatments/DP1_DEG_stats_T6C1.txt", header=T)
## #head(resGXE.vitA.1)
## resGXE.vitA.2 <- read.table("../2_Differential/DEG_results/treatments/DP4_DEG_stats_T6C1.txt", header=T)
## resGXE.vitA.3 <- read.table("../2_Differential/DEG_results/treatments/DP6_DEG_stats_T6C1.txt", header=T)
## resGXE.vitA.4 <- read.table("../2_Differential/DEG_results/treatments/DP9_DEG_stats_T6C1.txt", header=T)
## resGXE.vitA.5 <- read.table("../2_Differential/DEG_results/treatments/DP11_DEG_stats_T6C1.txt", header=T)


resGXE.vitA <- read.table("../2_Differential/DEG_results/treatments/DP6_DEG_stats_T6C1.txt", header=T)


resGXE.vitA[!is.na(resGXE.vitA$padj) & resGXE.vitA$padj<0.1 & resGXE.vitA$logFC > 0.5, "Significance"] <- "Up"
resGXE.vitA[!is.na(resGXE.vitA$padj) & resGXE.vitA$padj<0.1 & resGXE.vitA$logFC < -0.5, "Significance"] <- "Down"
resGXE.vitA[!is.na(resGXE.vitA$padj) & resGXE.vitA$padj>=0.1, "Significance"] <- "Not Significant"
resGXE.vitA[is.na(resGXE.vitA$padj),"Significance"] <- "NA"
head(resGXE.vitA)
nrow(resGXE.vitA)
table(resGXE.vitA$Significance)


resGXE.vitA[!is.na(resGXE.vitA$padj) & resGXE.vitA$padj<0.1 & resGXE.vitA$logFC > 0.25, "Significance.25"] <- "Up"
resGXE.vitA[!is.na(resGXE.vitA$padj) & resGXE.vitA$padj<0.1 & resGXE.vitA$logFC < -0.25, "Significance.25"] <- "Down"
resGXE.vitA[!is.na(resGXE.vitA$padj) & resGXE.vitA$padj>=0.1, "Significance.25"] <- "Not Significant"
resGXE.vitA[is.na(resGXE.vitA$padj),"Significance.25"] <- "NA"
head(resGXE.vitA)
nrow(resGXE.vitA)
table(resGXE.vitA$Significance.25)


resGXE.vitA[!is.na(resGXE.vitA$padj) & resGXE.vitA$padj<0.1 & resGXE.vitA$logFC > 0, "Significance.0"] <- "Up"
resGXE.vitA[!is.na(resGXE.vitA$padj) & resGXE.vitA$padj<0.1 & resGXE.vitA$logFC < 0, "Significance.0"] <- "Down"
resGXE.vitA[!is.na(resGXE.vitA$padj) & resGXE.vitA$padj>=0.1, "Significance.0"] <- "Not Significant"
resGXE.vitA[is.na(resGXE.vitA$padj),"Significance.0"] <- "NA"
head(resGXE.vitA)
nrow(resGXE.vitA)
table(resGXE.vitA$Significance.0)


## head(resDE)
## head(all.DE)
## length(all.DE)

## resGXE.vitA.filt <- resGXE.vitA %>% dplyr::filter(g.id %in% all.DE)
## head(resGXE.vitA.filt)
## nrow(resGXE.vitA.filt)
## length(unique(resGXE.vitA.filt$g.id))

## genes_summary <- as.data.frame(table(resGXE.vitA$g.id))
## head(genes_summary)
## nrow(genes_summary)
## #genes_summary <- resGXE.vitA %>% group_by(g.id) %>% summarise(sum = sum(g.id))

## genes_summary <- genes_summary %>% dplyr::filter(Freq==1)
## head(genes_summary)
## nrow(genes_summary)

#resGXE.vitA.sum <- resGXE.vitA %>% group_by(g.id) %>% summarise_each(funs(mean))


resGXE.vitA.sum <- resGXE.vitA %>% arrange(padj)
resGXE.vitA.sum <- resGXE.vitA.sum %>% distinct(ensg, .keep_all = TRUE)
head(resGXE.vitA.sum)
nrow(resGXE.vitA.sum)
table(resGXE.vitA.sum$Significance)
table(resGXE.vitA.sum$Significance.25)
table(resGXE.vitA.sum$Significance.0)


resGXE.vitA.sum <- resGXE.vitA %>% arrange(padj)
resGXE.vitA.sum <- resGXE.vitA.sum %>% distinct(g.id, .keep_all = TRUE)
head(resGXE.vitA.sum)
nrow(resGXE.vitA.sum)
table(resGXE.vitA.sum$Significance)
table(resGXE.vitA.sum$Significance.25)
table(resGXE.vitA.sum$Significance.0)


resGXE.vitA.sum <- resGXE.vitA %>% mutate(sig=ifelse(Significance.25 %in% c("Up", "Down"), "sig", "nonsig"))
#head(resGXE.vitA.sum)
#table(resGXE.vitA.sum$sig)
resGXE.vitA.sum <- resGXE.vitA.sum %>% arrange(desc(logFC))
resGXE.vitA.sum <- resGXE.vitA.sum %>% arrange(padj)
resGXE.vitA.sum <- resGXE.vitA.sum %>% arrange(desc(sig))
resGXE.vitA.sum <- resGXE.vitA.sum %>% distinct(ensg, .keep_all = TRUE)
head(resGXE.vitA.sum)
nrow(resGXE.vitA.sum)
table(resGXE.vitA.sum$sig)
table(resGXE.vitA.sum$Significance)
table(resGXE.vitA.sum$Significance.25)
table(resGXE.vitA.sum$Significance.0)


## newdf <- resGXE.vitA %>% group_by(g.id) %>% summarise(mean = mean(logFC))
## head(newdf)

## resGXE.vitA.sum$Significance <- NULL

## resGXE.vitA.sum[!is.na(resGXE.vitA.sum$padj) & resGXE.vitA.sum$padj<0.1 & resGXE.vitA.sum$logFC > 0.5, "Significance"] <- "Up"
## resGXE.vitA.sum[!is.na(resGXE.vitA.sum$padj) & resGXE.vitA.sum$padj<0.1 & resGXE.vitA.sum$logFC < -0.5, "Significance"] <- "Down"
## resGXE.vitA.sum[!is.na(resGXE.vitA.sum$padj) & resGXE.vitA.sum$padj>=0.1, "Significance"] <- "Not Significant"
## resGXE.vitA.sum[is.na(resGXE.vitA.sum$padj),"Significance"] <- "NA"
## head(resGXE.vitA.sum)
## nrow(resGXE.vitA.sum)
## table(resGXE.vitA.sum$Significance)


resGXE.vitA.sum <- resGXE.vitA.sum %>% mutate(contrast.gene=paste0("vitA.", g.id))
head(resGXE.vitA.sum)
nrow(resGXE.vitA.sum)
table(resGXE.vitA.sum$Significance)



#----joining resDE and resGXE.vitA----
#resDE <- resDE %>% mutate(contrast.gene=paste0(contrast, ".", gene))
table(resDE$contrast)

resDE.vitA <- resDE %>% dplyr::filter(contrast=="vitA")
head(resDE.vitA)
nrow(resDE.vitA)
#nrow(sig_resDE.vitA)
#table(sig_resDE.vitA$contrast)


resDE.GXE.vitA.sum <- left_join(resDE.vitA,resGXE.vitA.sum, by="contrast.gene")
head(resDE.GXE.vitA.sum)
nrow(resDE.GXE.vitA.sum)
table(resDE.GXE.vitA.sum$Significance.x)
table(resDE.GXE.vitA.sum$Significance.y)


resDE.GXE.vitA.sum <- resDE.GXE.vitA.sum %>% filter(!is.na(logFC))
#%>% filter(!is.na(estimate.x))
head(resDE.GXE.vitA.sum)
nrow(resDE.GXE.vitA.sum)
table(resDE.GXE.vitA.sum$Significance.x)
table(resDE.GXE.vitA.sum$Significance.y)


resDE.GXE.vitA.sum <- resDE.GXE.vitA.sum %>% mutate(Significance=
             ifelse(Significance.y %in% c("Up", "Down") &  Significance.x %in% c("Up", "Down"), "both",
             ifelse(Significance.y %in% c("Up", "Down") & (Significance.x %in% c("Not Significant", "NA") | is.na(Significance.x)), "GXE",
             ifelse((Significance.y %in% c("Not Significant", "NA") | is.na(Significance.y)) & Significance.x %in% c("Up", "Down"), "sc", "NA"))))
head(resDE.GXE.vitA.sum)
nrow(resDE.GXE.vitA.sum)
table(resDE.GXE.vitA.sum$Significance)


resDE.GXE.vitA.sum <- resDE.GXE.vitA.sum %>% mutate(Significance.25=
             ifelse(Significance.25.y %in% c("Up", "Down") &  Significance.25.x %in% c("Up", "Down"), "both",
             ifelse(Significance.25.y %in% c("Up", "Down") & (Significance.25.x %in% c("Not Significant", "NA") | is.na(Significance.25.x)), "GXE",
             ifelse((Significance.25.y %in% c("Not Significant", "NA") | is.na(Significance.25.y)) & Significance.25.x %in% c("Up", "Down"), "sc", "NA"))))
head(resDE.GXE.vitA.sum)
nrow(resDE.GXE.vitA.sum)
table(resDE.GXE.vitA.sum$Significance.25)



resDE.GXE.vitA.sum <- resDE.GXE.vitA.sum %>% mutate(Significance.0=
             ifelse(Significance.0.y %in% c("Up", "Down") &  Significance.0.x %in% c("Up", "Down"), "both",
             ifelse(Significance.0.y %in% c("Up", "Down") & (Significance.0.x %in% c("Not Significant", "NA") | is.na(Significance.0.x)), "GXE",
             ifelse((Significance.0.y %in% c("Not Significant", "NA") | is.na(Significance.0.y)) & Significance.0.x %in% c("Up", "Down"), "sc", "NA"))))
head(resDE.GXE.vitA.sum)
nrow(resDE.GXE.vitA.sum)
table(resDE.GXE.vitA.sum$Significance.0)



## resDE.GXE.vitA.sum <- resDE.GXE.vitA.sum %>% dplyr::filter(!is.na(Significance))
## head(resDE.GXE.vitA.sum)
## nrow(resDE.GXE.vitA.sum)
## table(resDE.GXE.vitA.sum$Significance)




#----vitA plot----
logFC <- 0

MCls <- unique(resDE.GXE.vitA.sum$MCls)
MCls

comb.2 <- MCls

for(i in comb.2){
            res1 <- resDE.GXE.vitA.sum %>% dplyr::filter(MCls==i)
            png(paste0("./1_GXE.outs/logFC", logFC, "_all_scatter_vitA_", i, ".png"))
            p <- ggplot(res1, aes(x=logFC, y=estimate, color=Significance.0)) +
                scale_colour_manual(values = c("#F8766D", "#00BA38", "darkgrey", "#619CFF")) +
                    geom_point(show.legend = TRUE)+
                    xlab("GXE_data") + ylab("sc_data")+
#                    xlim(-4, 10) + ylim(-4, 5)+
                    ggtitle(paste0("logFC_", logFC, "_all_logFC_vitA_", i))+
                    theme_bw()+
                    theme(plot.title = element_text(color = "black", size = 18, face = "bold"),
                          axis.text = element_text(size = 18),
                          axis.title = element_text(size = 18))
            print(p)
            dev.off()
        }



















#----vitD----
#resGXE.vitD.1 <- read.table("../2_Differential/DEG_results/treatments/DP4_DEG_stats_T23C1.txt", header=T)

resGXE.vitD <- read.table("../2_Differential/DEG_results/treatments/DP6_DEG_stats_T23C1.txt", header=T)

#resGXE.vitD.3 <- read.table("../2_Differential/DEG_results/treatments/DP9_DEG_stats_T23C1.txt", header=T)
#resGXE.vitD.4 <- read.table("../2_Differential/DEG_results/treatments/DP11_DEG_stats_T23C1.txt", header=T)


resGXE.vitD[!is.na(resGXE.vitD$padj) & resGXE.vitD$padj<0.1 & resGXE.vitD$logFC > 0.5, "Significance"] <- "Up"
resGXE.vitD[!is.na(resGXE.vitD$padj) & resGXE.vitD$padj<0.1 & resGXE.vitD$logFC < -0.5, "Significance"] <- "Down"
resGXE.vitD[!is.na(resGXE.vitD$padj) & resGXE.vitD$padj>=0.1, "Significance"] <- "Not Significant"
resGXE.vitD[is.na(resGXE.vitD$padj),"Significance"] <- "NA"
head(resGXE.vitD)
nrow(resGXE.vitD)
table(resGXE.vitD$Significance)


resGXE.vitD[!is.na(resGXE.vitD$padj) & resGXE.vitD$padj<0.1 & resGXE.vitD$logFC > 0.25, "Significance.25"] <- "Up"
resGXE.vitD[!is.na(resGXE.vitD$padj) & resGXE.vitD$padj<0.1 & resGXE.vitD$logFC < -0.25, "Significance.25"] <- "Down"
resGXE.vitD[!is.na(resGXE.vitD$padj) & resGXE.vitD$padj>=0.1, "Significance.25"] <- "Not Significant"
resGXE.vitD[is.na(resGXE.vitD$padj),"Significance.25"] <- "NA"
head(resGXE.vitD)
nrow(resGXE.vitD)
table(resGXE.vitD$Significance.25)


resGXE.vitD[!is.na(resGXE.vitD$padj) & resGXE.vitD$padj<0.1 & resGXE.vitD$logFC > 0, "Significance.0"] <- "Up"
resGXE.vitD[!is.na(resGXE.vitD$padj) & resGXE.vitD$padj<0.1 & resGXE.vitD$logFC < 0, "Significance.0"] <- "Down"
resGXE.vitD[!is.na(resGXE.vitD$padj) & resGXE.vitD$padj>=0.1, "Significance.0"] <- "Not Significant"
resGXE.vitD[is.na(resGXE.vitD$padj),"Significance.0"] <- "NA"
head(resGXE.vitD)
nrow(resGXE.vitD)
table(resGXE.vitD$Significance.0)


## head(resDE)
## head(all.DE)
## length(all.DE)

## resGXE.vitD.filt <- resGXE.vitD %>% dplyr::filter(g.id %in% all.DE)
## head(resGXE.vitD.filt)
## nrow(resGXE.vitD.filt)
## length(unique(resGXE.vitD.filt$g.id))

## genes_summary <- as.data.frame(table(resGXE.vitD$g.id))
## head(genes_summary)
## nrow(genes_summary)
## #genes_summary <- resGXE.vitD %>% group_by(g.id) %>% summarise(sum = sum(g.id))

## genes_summary <- genes_summary %>% dplyr::filter(Freq==1)
## head(genes_summary)
## nrow(genes_summary)

#resGXE.vitD.sum <- resGXE.vitD %>% group_by(g.id) %>% summarise_each(funs(mean))

resGXE.vitD.sum <- resGXE.vitD %>% arrange(padj)
resGXE.vitD.sum <- resGXE.vitD.sum %>% distinct(ensg, .keep_all = TRUE)
head(resGXE.vitD.sum)
nrow(resGXE.vitD.sum)
table(resGXE.vitD.sum$Significance)
table(resGXE.vitD.sum$Significance.25)
table(resGXE.vitD.sum$Significance.0)

resGXE.vitD.sum <- resGXE.vitD %>% arrange(padj)
resGXE.vitD.sum <- resGXE.vitD.sum %>% distinct(g.id, .keep_all = TRUE)
head(resGXE.vitD.sum)
nrow(resGXE.vitD.sum)
table(resGXE.vitD.sum$Significance)
table(resGXE.vitD.sum$Significance.25)
table(resGXE.vitD.sum$Significance.0)


resGXE.vitD.sum <- resGXE.vitD %>% mutate(sig=ifelse(Significance.25 %in% c("Up", "Down"), "sig", "nonsig"))
#head(resGXE.vitD.sum)
#table(resGXE.vitD.sum$sig)
resGXE.vitD.sum <- resGXE.vitD.sum %>% arrange(desc(logFC))
resGXE.vitD.sum <- resGXE.vitD.sum %>% arrange(padj)
resGXE.vitD.sum <- resGXE.vitD.sum %>% arrange(desc(sig))
resGXE.vitD.sum <- resGXE.vitD.sum %>% distinct(ensg, .keep_all = TRUE)
head(resGXE.vitD.sum)
nrow(resGXE.vitD.sum)
table(resGXE.vitD.sum$sig)
table(resGXE.vitD.sum$Significance)
table(resGXE.vitD.sum$Significance.25)
table(resGXE.vitD.sum$Significance.0)


## newdf <- resGXE.vitD %>% group_by(g.id) %>% summarise(mean = mean(logFC))
## head(newdf)

## resGXE.vitD.sum$Significance <- NULL

## resGXE.vitD.sum[!is.na(resGXE.vitD.sum$padj) & resGXE.vitD.sum$padj<0.1 & resGXE.vitD.sum$logFC > 0.5, "Significance"] <- "Up"
## resGXE.vitD.sum[!is.na(resGXE.vitD.sum$padj) & resGXE.vitD.sum$padj<0.1 & resGXE.vitD.sum$logFC < -0.5, "Significance"] <- "Down"
## resGXE.vitD.sum[!is.na(resGXE.vitD.sum$padj) & resGXE.vitD.sum$padj>=0.1, "Significance"] <- "Not Significant"
## resGXE.vitD.sum[is.na(resGXE.vitD.sum$padj),"Significance"] <- "NA"
## head(resGXE.vitD.sum)
## nrow(resGXE.vitD.sum)
## table(resGXE.vitD.sum$Significance)


resGXE.vitD.sum <- resGXE.vitD.sum %>% mutate(contrast.gene=paste0("vitD.", g.id))
head(resGXE.vitD.sum)
nrow(resGXE.vitD.sum)
table(resGXE.vitD.sum$Significance)



#----joining resDE and resGXE.vitD----
#resDE <- resDE %>% mutate(contrast.gene=paste0(contrast, ".", gene))
table(resDE$contrast)

resDE.vitD <- resDE %>% dplyr::filter(contrast=="vitD")
head(resDE.vitD)
nrow(resDE.vitD)
#nrow(sig_resDE.vitD)
#table(sig_resDE.vitD$contrast)


resDE.GXE.vitD.sum <- left_join(resDE.vitD,resGXE.vitD.sum, by="contrast.gene")
head(resDE.GXE.vitD.sum)
nrow(resDE.GXE.vitD.sum)
table(resDE.GXE.vitD.sum$Significance.x)
table(resDE.GXE.vitD.sum$Significance.y)


resDE.GXE.vitD.sum <- resDE.GXE.vitD.sum %>% filter(!is.na(logFC))
#%>% filter(!is.na(estimate.x))
head(resDE.GXE.vitD.sum)
nrow(resDE.GXE.vitD.sum)
table(resDE.GXE.vitD.sum$Significance.x)
table(resDE.GXE.vitD.sum$Significance.y)

resDE.GXE.vitD.sum <- resDE.GXE.vitD.sum %>% mutate(Significance=
             ifelse(Significance.y %in% c("Up", "Down") &  Significance.x %in% c("Up", "Down"), "both",
             ifelse(Significance.y %in% c("Up", "Down") & (Significance.x %in% c("Not Significant", "NA") | is.na(Significance.x)), "GXE",
             ifelse((Significance.y %in% c("Not Significant", "NA") | is.na(Significance.y)) & Significance.x %in% c("Up", "Down"), "sc", "NA"))))
head(resDE.GXE.vitD.sum)
nrow(resDE.GXE.vitD.sum)
table(resDE.GXE.vitD.sum$Significance)


resDE.GXE.vitD.sum <- resDE.GXE.vitD.sum %>% mutate(Significance.25=
             ifelse(Significance.25.y %in% c("Up", "Down") &  Significance.25.x %in% c("Up", "Down"), "both",
             ifelse(Significance.25.y %in% c("Up", "Down") & (Significance.25.x %in% c("Not Significant", "NA") | is.na(Significance.25.x)), "GXE",
             ifelse((Significance.25.y %in% c("Not Significant", "NA") | is.na(Significance.25.y)) & Significance.25.x %in% c("Up", "Down"), "sc", "NA"))))
head(resDE.GXE.vitD.sum)
nrow(resDE.GXE.vitD.sum)
table(resDE.GXE.vitD.sum$Significance.25)



resDE.GXE.vitD.sum <- resDE.GXE.vitD.sum %>% mutate(Significance.0=
             ifelse(Significance.0.y %in% c("Up", "Down") &  Significance.0.x %in% c("Up", "Down"), "both",
             ifelse(Significance.0.y %in% c("Up", "Down") & (Significance.0.x %in% c("Not Significant", "NA") | is.na(Significance.0.x)), "GXE",
             ifelse((Significance.0.y %in% c("Not Significant", "NA") | is.na(Significance.0.y)) & Significance.0.x %in% c("Up", "Down"), "sc", "NA"))))
head(resDE.GXE.vitD.sum)
nrow(resDE.GXE.vitD.sum)
table(resDE.GXE.vitD.sum$Significance.0)



## resDE.GXE.vitD.sum <- resDE.GXE.vitD.sum %>% dplyr::filter(!is.na(Significance))
## head(resDE.GXE.vitD.sum)
## nrow(resDE.GXE.vitD.sum)
## table(resDE.GXE.vitD.sum$Significance)




#----vitD plot----
logFC <- 0

MCls <- unique(resDE.GXE.vitD.sum$MCls)
MCls

comb.2 <- MCls

for(i in comb.2){
            res1 <- resDE.GXE.vitD.sum %>% dplyr::filter(MCls==i)
            png(paste0("./1_GXE.outs/logFC_", logFC, "_all_scatter_vitD_", i, ".png"))
            p <- ggplot(res1, aes(x=logFC, y=estimate, color=Significance.0)) +
                scale_colour_manual(values = c("#F8766D", "#00BA38", "darkgrey", "#619CFF")) +
                    geom_point(show.legend = TRUE)+
                    xlab("GXE_data") + ylab("sc_data")+
#                    xlim(-4, 10) + ylim(-4, 5)+
                    ggtitle(paste0("logFC_", logFC, "_all_logFC_vitD_", i))+
                    theme_bw()+
                    theme(plot.title = element_text(color = "black", size = 18, face = "bold"),
                          axis.text = element_text(size = 18),
                          axis.title = element_text(size = 18))
            print(p)
            dev.off()
        }


















#----vitE----
resGXE.vitE <- read.table("../2_Differential/DEG_results/treatments/DP6_DEG_stats_T7C1.txt", header=T)

nrow(resGXE.vitE)

#resGXE.vitE.2 <- read.table("../2_Differential/DEG_results/treatments/DP11_DEG_stats_T7C1.txt", header=T)


resGXE.vitE[!is.na(resGXE.vitE$padj) & resGXE.vitE$padj<0.1 & resGXE.vitE$logFC > 0.5, "Significance"] <- "Up"
resGXE.vitE[!is.na(resGXE.vitE$padj) & resGXE.vitE$padj<0.1 & resGXE.vitE$logFC < -0.5, "Significance"] <- "Down"
resGXE.vitE[!is.na(resGXE.vitE$padj) & resGXE.vitE$padj>=0.1, "Significance"] <- "Not Significant"
resGXE.vitE[is.na(resGXE.vitE$padj),"Significance"] <- "NA"
head(resGXE.vitE)
nrow(resGXE.vitE)
table(resGXE.vitE$Significance)


resGXE.vitE[!is.na(resGXE.vitE$padj) & resGXE.vitE$padj<0.1 & resGXE.vitE$logFC > 0.25, "Significance.25"] <- "Up"
resGXE.vitE[!is.na(resGXE.vitE$padj) & resGXE.vitE$padj<0.1 & resGXE.vitE$logFC < -0.25, "Significance.25"] <- "Down"
resGXE.vitE[!is.na(resGXE.vitE$padj) & resGXE.vitE$padj>=0.1, "Significance.25"] <- "Not Significant"
resGXE.vitE[is.na(resGXE.vitE$padj),"Significance.25"] <- "NA"
head(resGXE.vitE)
nrow(resGXE.vitE)
table(resGXE.vitE$Significance.25)


resGXE.vitE[!is.na(resGXE.vitE$padj) & resGXE.vitE$padj<0.1 & resGXE.vitE$logFC > 0, "Significance.0"] <- "Up"
resGXE.vitE[!is.na(resGXE.vitE$padj) & resGXE.vitE$padj<0.1 & resGXE.vitE$logFC < 0, "Significance.0"] <- "Down"
resGXE.vitE[!is.na(resGXE.vitE$padj) & resGXE.vitE$padj>=0.1, "Significance.0"] <- "Not Significant"
resGXE.vitE[is.na(resGXE.vitE$padj),"Significance.0"] <- "NA"
head(resGXE.vitE)
nrow(resGXE.vitE)
table(resGXE.vitE$Significance.0)


## head(resDE)
## head(all.DE)
## length(all.DE)

## resGXE.vitE.filt <- resGXE.vitE %>% dplyr::filter(g.id %in% all.DE)
## head(resGXE.vitE.filt)
## nrow(resGXE.vitE.filt)
## length(unique(resGXE.vitE.filt$g.id))

## genes_summary <- as.data.frame(table(resGXE.vitE$g.id))
## head(genes_summary)
## nrow(genes_summary)
## #genes_summary <- resGXE.vitE %>% group_by(g.id) %>% summarise(sum = sum(g.id))

## genes_summary <- genes_summary %>% dplyr::filter(Freq==1)
## head(genes_summary)
## nrow(genes_summary)

#resGXE.vitE.sum <- resGXE.vitE %>% group_by(g.id) %>% summarise_each(funs(mean))

resGXE.vitE.sum <- resGXE.vitE %>% arrange(padj)
resGXE.vitE.sum <- resGXE.vitE.sum %>% distinct(ensg, .keep_all = TRUE)
head(resGXE.vitE.sum)
nrow(resGXE.vitE.sum)
table(resGXE.vitE.sum$Significance)
table(resGXE.vitE.sum$Significance.25)
table(resGXE.vitE.sum$Significance.0)

resGXE.vitE.sum <- resGXE.vitE %>% arrange(padj)
resGXE.vitE.sum <- resGXE.vitE.sum %>% distinct(g.id, .keep_all = TRUE)
head(resGXE.vitE.sum)
nrow(resGXE.vitE.sum)
table(resGXE.vitE.sum$Significance)
table(resGXE.vitE.sum$Significance.25)
table(resGXE.vitE.sum$Significance.0)


resGXE.vitE.sum <- resGXE.vitE %>% mutate(sig=ifelse(Significance.25 %in% c("Up", "Down"), "sig", "nonsig"))
#head(resGXE.vitE.sum)
#table(resGXE.vitE.sum$sig)
resGXE.vitE.sum <- resGXE.vitE.sum %>% arrange(desc(logFC))
resGXE.vitE.sum <- resGXE.vitE.sum %>% arrange(padj)
resGXE.vitE.sum <- resGXE.vitE.sum %>% arrange(desc(sig))
resGXE.vitE.sum <- resGXE.vitE.sum %>% distinct(ensg, .keep_all = TRUE)
head(resGXE.vitE.sum)
nrow(resGXE.vitE.sum)
table(resGXE.vitE.sum$sig)
table(resGXE.vitE.sum$Significance)
table(resGXE.vitE.sum$Significance.25)
table(resGXE.vitE.sum$Significance.0)


## newdf <- resGXE.vitE %>% group_by(g.id) %>% summarise(mean = mean(logFC))
## head(newdf)

## resGXE.vitE.sum$Significance <- NULL

## resGXE.vitE.sum[!is.na(resGXE.vitE.sum$padj) & resGXE.vitE.sum$padj<0.1 & resGXE.vitE.sum$logFC > 0.5, "Significance"] <- "Up"
## resGXE.vitE.sum[!is.na(resGXE.vitE.sum$padj) & resGXE.vitE.sum$padj<0.1 & resGXE.vitE.sum$logFC < -0.5, "Significance"] <- "Down"
## resGXE.vitE.sum[!is.na(resGXE.vitE.sum$padj) & resGXE.vitE.sum$padj>=0.1, "Significance"] <- "Not Significant"
## resGXE.vitE.sum[is.na(resGXE.vitE.sum$padj),"Significance"] <- "NA"
## head(resGXE.vitE.sum)
## nrow(resGXE.vitE.sum)
## table(resGXE.vitE.sum$Significance)


resGXE.vitE.sum <- resGXE.vitE.sum %>% mutate(contrast.gene=paste0("vitE.", g.id))
head(resGXE.vitE.sum)
nrow(resGXE.vitE.sum)
table(resGXE.vitE.sum$Significance)



#----joining resDE and resGXE.vitE----
resDE <- resDE %>% mutate(contrast.gene=paste0(contrast, ".", gene))
table(resDE$contrast)

resDE.vitE <- resDE %>% dplyr::filter(contrast=="vitE")
head(resDE.vitE)
nrow(resDE.vitE)
#nrow(sig_resDE.vitE)
#table(sig_resDE.vitE$contrast)


resDE.GXE.vitE.sum <- left_join(resDE.vitE,resGXE.vitE.sum, by="contrast.gene")
head(resDE.GXE.vitE.sum)
nrow(resDE.GXE.vitE.sum)
table(resDE.GXE.vitE.sum$Significance.x)
table(resDE.GXE.vitE.sum$Significance.y)


resDE.GXE.vitE.sum <- resDE.GXE.vitE.sum %>% filter(!is.na(logFC))
#%>% filter(!is.na(estimate.x))
head(resDE.GXE.vitE.sum)
nrow(resDE.GXE.vitE.sum)
table(resDE.GXE.vitE.sum$Significance.x)
table(resDE.GXE.vitE.sum$Significance.y)

resDE.GXE.vitE.sum <- resDE.GXE.vitE.sum %>% mutate(Significance=
             ifelse(Significance.y %in% c("Up", "Down") &  Significance.x %in% c("Up", "Down"), "both",
             ifelse(Significance.y %in% c("Up", "Down") & (Significance.x %in% c("Not Significant", "NA") | is.na(Significance.x)), "GXE",
             ifelse((Significance.y %in% c("Not Significant", "NA") | is.na(Significance.y)) & Significance.x %in% c("Up", "Down"), "sc", "NA"))))
head(resDE.GXE.vitE.sum)
nrow(resDE.GXE.vitE.sum)
table(resDE.GXE.vitE.sum$Significance)


resDE.GXE.vitE.sum <- resDE.GXE.vitE.sum %>% mutate(Significance.25=
             ifelse(Significance.25.y %in% c("Up", "Down") &  Significance.25.x %in% c("Up", "Down"), "both",
             ifelse(Significance.25.y %in% c("Up", "Down") & (Significance.25.x %in% c("Not Significant", "NA") | is.na(Significance.25.x)), "GXE",
             ifelse((Significance.25.y %in% c("Not Significant", "NA") | is.na(Significance.25.y)) & Significance.25.x %in% c("Up", "Down"), "sc", "NA"))))
head(resDE.GXE.vitE.sum)
nrow(resDE.GXE.vitE.sum)
table(resDE.GXE.vitE.sum$Significance.25)



resDE.GXE.vitE.sum <- resDE.GXE.vitE.sum %>% mutate(Significance.0=
             ifelse(Significance.0.y %in% c("Up", "Down") &  Significance.0.x %in% c("Up", "Down"), "both",
             ifelse(Significance.0.y %in% c("Up", "Down") & (Significance.0.x %in% c("Not Significant", "NA") | is.na(Significance.0.x)), "GXE",
             ifelse((Significance.0.y %in% c("Not Significant", "NA") | is.na(Significance.0.y)) & Significance.0.x %in% c("Up", "Down"), "sc", "NA"))))
head(resDE.GXE.vitE.sum)
nrow(resDE.GXE.vitE.sum)
table(resDE.GXE.vitE.sum$Significance.0)



## resDE.GXE.vitE.sum <- resDE.GXE.vitE.sum %>% dplyr::filter(!is.na(Significance))
## head(resDE.GXE.vitE.sum)
## nrow(resDE.GXE.vitE.sum)
## table(resDE.GXE.vitE.sum$Significance)




#----vitE plot----
logFC <- 0

MCls <- unique(resDE.GXE.vitE.sum$MCls)
MCls

comb.2 <- MCls


sig.color=c("both"="#F8766D",
            "GXE"="#00BA38",
            "NA"="darkgrey",
            "sc"="#619CFF")


for(i in comb.2){
            res1 <- resDE.GXE.vitE.sum %>% dplyr::filter(MCls==i)
            png(paste0("./1_GXE.outs/logFC_", logFC, "_all_scatter_vitE_", i, ".png"))
#            p <- ggplot(res1, aes(x=logFC, y=estimate, color=Significance)) +
#                scale_colour_manual(values = c("#F8766D", "#00BA38", "darkgrey", "#619CFF")) +
            p <- ggplot(res1, aes(x=logFC, y=estimate, colour=factor(Significance.0))) +
            scale_colour_manual(values = sig.color) +
                    geom_point(show.legend = TRUE)+
                    xlab("GXE_data") + ylab("sc_data")+
#                    xlim(-4, 10) + ylim(-4, 5)+
                    ggtitle(paste0("logFC_", logFC, "_all_logFC_vitE_", i))+
                    theme_bw()+
                    theme(plot.title = element_text(color = "black", size = 18, face = "bold"),
                          axis.text = element_text(size = 18),
                          axis.title = element_text(size = 18))
            print(p)
            dev.off()
        }























#----zinc----
resGXE.zinc <- read.table("../2_Differential/DEG_results/treatments/DP6_DEG_stats_T20C1.txt", header=T)

head(resGXE.zinc)
nrow(resGXE.zinc)
#head(resGXE.zinc %>% dplyr::filter(g.id=="LINC01409"))
#nrow(resGXE.zinc %>% dplyr::filter(g.id=="LINC01409"))
length(unique(resGXE.zinc$g.id))


resGXE.zinc.2 <- resGXE.zinc %>% mutate(Significance.check=
             ifelse(padj < 0.1 & logFC > 0.5, "Up",
             ifelse(padj < 0.1 & logFC < -0.5, "Down",
             ifelse(padj >= 0.1, "Not Signicant", "NA"))))
table(resGXE.zinc.2$Significance.check)






resGXE.zinc[!is.na(resGXE.zinc$padj) & resGXE.zinc$padj<0.1 & resGXE.zinc$logFC > 0.5, "Significance"] <- "Up"
resGXE.zinc[!is.na(resGXE.zinc$padj) & resGXE.zinc$padj<0.1 & resGXE.zinc$logFC < -0.5, "Significance"] <- "Down"
resGXE.zinc[!is.na(resGXE.zinc$padj) & resGXE.zinc$padj>=0.1, "Significance"] <- "Not Significant"
resGXE.zinc[is.na(resGXE.zinc$padj),"Significance"] <- "NA"
head(resGXE.zinc)
nrow(resGXE.zinc)
table(resGXE.zinc$Significance)


resGXE.zinc[!is.na(resGXE.zinc$padj) & resGXE.zinc$padj<0.1 & resGXE.zinc$logFC > 0.25, "Significance.25"] <- "Up"
resGXE.zinc[!is.na(resGXE.zinc$padj) & resGXE.zinc$padj<0.1 & resGXE.zinc$logFC < -0.25, "Significance.25"] <- "Down"
resGXE.zinc[!is.na(resGXE.zinc$padj) & resGXE.zinc$padj>=0.1, "Significance.25"] <- "Not Significant"
resGXE.zinc[is.na(resGXE.zinc$padj),"Significance.25"] <- "NA"
head(resGXE.zinc)
nrow(resGXE.zinc)
table(resGXE.zinc$Significance.25)

resGXE.zinc[!is.na(resGXE.zinc$padj) & resGXE.zinc$padj<0.1 & resGXE.zinc$logFC > 0, "Significance.0"] <- "Up"
resGXE.zinc[!is.na(resGXE.zinc$padj) & resGXE.zinc$padj<0.1 & resGXE.zinc$logFC < 0, "Significance.0"] <- "Down"
resGXE.zinc[!is.na(resGXE.zinc$padj) & resGXE.zinc$padj>=0.1, "Significance.0"] <- "Not Significant"
resGXE.zinc[is.na(resGXE.zinc$padj),"Significance.0"] <- "NA"
head(resGXE.zinc)
nrow(resGXE.zinc)
table(resGXE.zinc$Significance.0)



## head(resDE)
## head(all.DE)
## length(all.DE)

## resGXE.zinc.filt <- resGXE.zinc %>% dplyr::filter(g.id %in% all.DE)
## head(resGXE.zinc.filt)
## nrow(resGXE.zinc.filt)
## length(unique(resGXE.zinc.filt$g.id))

## genes_summary <- as.data.frame(table(resGXE.zinc$g.id))
## head(genes_summary)
## nrow(genes_summary)
## #genes_summary <- resGXE.zinc %>% group_by(g.id) %>% summarise(sum = sum(g.id))

## genes_summary <- genes_summary %>% dplyr::filter(Freq==1)
## head(genes_summary)
## nrow(genes_summary)

#resGXE.zinc.sum <- resGXE.zinc %>% group_by(g.id) %>% summarise_each(funs(mean))

resGXE.zinc.sum <- resGXE.zinc %>% arrange(padj)
resGXE.zinc.sum <- resGXE.zinc.sum %>% distinct(ensg, .keep_all = TRUE)
head(resGXE.zinc.sum)
nrow(resGXE.zinc.sum)
table(resGXE.zinc.sum$Significance)
table(resGXE.zinc.sum$Significance.25)
table(resGXE.zinc.sum$Significance.0)

resGXE.zinc.sum <- resGXE.zinc %>% arrange(padj)
resGXE.zinc.sum <- resGXE.zinc.sum %>% distinct(g.id, .keep_all = TRUE)
head(resGXE.zinc.sum)
nrow(resGXE.zinc.sum)
table(resGXE.zinc.sum$Significance)
table(resGXE.zinc.sum$Significance.25)
table(resGXE.zinc.sum$Significance.0)


resGXE.zinc.sum <- resGXE.zinc %>% mutate(sig=ifelse(Significance.25 %in% c("Up", "Down"), "sig", "nonsig"))
#head(resGXE.zinc.sum)
#table(resGXE.zinc.sum$sig)
resGXE.zinc.sum <- resGXE.zinc.sum %>% arrange(desc(logFC))
resGXE.zinc.sum <- resGXE.zinc.sum %>% arrange(padj)
resGXE.zinc.sum <- resGXE.zinc.sum %>% arrange(desc(sig))
resGXE.zinc.sum <- resGXE.zinc.sum %>% distinct(ensg, .keep_all = TRUE)
head(resGXE.zinc.sum)
nrow(resGXE.zinc.sum)
table(resGXE.zinc.sum$sig)
table(resGXE.zinc.sum$Significance)
table(resGXE.zinc.sum$Significance.25)
table(resGXE.zinc.sum$Significance.0)

## newdf <- resGXE.zinc %>% group_by(g.id) %>% summarise(mean = mean(logFC))
## head(newdf)

## resGXE.zinc.sum[!is.na(resGXE.zinc.sum$padj) & resGXE.zinc.sum$padj<0.1 & resGXE.zinc.sum$logFC > 0.5, "Significance"] <- "Up"
## resGXE.zinc.sum[!is.na(resGXE.zinc.sum$padj) & resGXE.zinc.sum$padj<0.1 & resGXE.zinc.sum$logFC < -0.5, "Significance"] <- "Down"
## resGXE.zinc.sum[!is.na(resGXE.zinc.sum$padj) & resGXE.zinc.sum$padj>=0.1, "Significance"] <- "Not Significant"
## resGXE.zinc.sum[is.na(resGXE.zinc.sum$padj),"Significance"] <- "NA"
## head(resGXE.zinc.sum)
## nrow(resGXE.zinc.sum)
## table(resGXE.zinc.sum$Significance)


resGXE.zinc.sum <- resGXE.zinc.sum %>% mutate(contrast.gene=paste0("zinc.", g.id))
head(resGXE.zinc.sum)
nrow(resGXE.zinc.sum)
table(resGXE.zinc.sum$Significance)



#----joining resDE and resGXE.zinc----
#resDE <- resDE %>% mutate(contrast.gene=paste0(contrast, ".", gene))
resDE.zinc <- resDE %>% dplyr::filter(contrast=="zinc")
head(resDE.zinc)
nrow(resDE.zinc)
#nrow(sig_resDE.zinc)
#table(sig_resDE.zinc$contrast)


resDE.GXE.zinc.sum <- left_join(resDE.zinc,resGXE.zinc.sum, by="contrast.gene")
head(resDE.GXE.zinc.sum)
nrow(resDE.GXE.zinc.sum)
table(resDE.GXE.zinc.sum$Significance.x)
table(resDE.GXE.zinc.sum$Significance.y)


resDE.GXE.zinc.sum <- resDE.GXE.zinc.sum %>% filter(!is.na(logFC))
#%>% filter(!is.na(estimate.x))
head(resDE.GXE.zinc.sum)
nrow(resDE.GXE.zinc.sum)
table(resDE.GXE.zinc.sum$Significance.x)
table(resDE.GXE.zinc.sum$Significance.y)


resDE.GXE.zinc.sum <- resDE.GXE.zinc.sum %>% mutate(Significance=
             ifelse(Significance.y %in% c("Up", "Down") &  Significance.x %in% c("Up", "Down"), "both",
             ifelse(Significance.y %in% c("Up", "Down") & (Significance.x %in% c("Not Significant", "NA") | is.na(Significance.x)), "GXE",
             ifelse((Significance.y %in% c("Not Significant", "NA") | is.na(Significance.y)) & Significance.x %in% c("Up", "Down"), "sc", "NA"))))
head(resDE.GXE.zinc.sum)
nrow(resDE.GXE.zinc.sum)
table(resDE.GXE.zinc.sum$Significance)


resDE.GXE.zinc.sum <- resDE.GXE.zinc.sum %>% mutate(Significance.25=
             ifelse(Significance.25.y %in% c("Up", "Down") &  Significance.25.x %in% c("Up", "Down"), "both",
             ifelse(Significance.25.y %in% c("Up", "Down") & (Significance.25.x %in% c("Not Significant", "NA") | is.na(Significance.25.x)), "GXE",
             ifelse((Significance.25.y %in% c("Not Significant", "NA") | is.na(Significance.25.y)) & Significance.25.x %in% c("Up", "Down"), "sc", "NA"))))
head(resDE.GXE.zinc.sum)
nrow(resDE.GXE.zinc.sum)
table(resDE.GXE.zinc.sum$Significance.25)



resDE.GXE.zinc.sum <- resDE.GXE.zinc.sum %>% mutate(Significance.0=
             ifelse(Significance.0.y %in% c("Up", "Down") &  Significance.0.x %in% c("Up", "Down"), "both",
             ifelse(Significance.0.y %in% c("Up", "Down") & (Significance.0.x %in% c("Not Significant", "NA") | is.na(Significance.0.x)), "GXE",
             ifelse((Significance.0.y %in% c("Not Significant", "NA") | is.na(Significance.0.y)) & Significance.0.x %in% c("Up", "Down"), "sc", "NA"))))
head(resDE.GXE.zinc.sum)
nrow(resDE.GXE.zinc.sum)
table(resDE.GXE.zinc.sum$Significance.0)


## resDE.GXE.zinc.sum <- resDE.GXE.zinc.sum %>% dplyr::filter(!is.na(Significance))
## head(resDE.GXE.zinc.sum)
## nrow(resDE.GXE.zinc.sum)
## table(resDE.GXE.zinc.sum$Significance)




#----zinc plot----
logFC <- 0.25

MCls <- unique(resDE.GXE.zinc.sum$MCls)
MCls

comb.2 <- MCls

for(i in comb.2){
            res1 <- resDE.GXE.zinc.sum %>% dplyr::filter(MCls==i)
            png(paste0("./1_GXE.outs/logFC_", logFC, "_all_scatter_zinc_", i, ".png"))
            p <- ggplot(res1, aes(x=logFC, y=estimate, color=Significance.25)) +
                scale_colour_manual(values = c("#F8766D", "#00BA38", "darkgrey", "#619CFF")) +
                    geom_point(show.legend = TRUE)+
                    xlab("GXE_data") + ylab("sc_data")+
#                    xlim(-4, 10) + ylim(-4, 5)+
                    ggtitle(paste0("logFC_", logFC, "_all_logFC_zinc_", i))+
                    theme_bw()+
                    theme(plot.title = element_text(color = "black", size = 18, face = "bold"),
                          axis.text = element_text(size = 18),
                          axis.title = element_text(size = 18))
            print(p)
            dev.off()
        }







#----spearman_corr----
scorr_caff <- c()
scorr_nic <- c()
scorr_vitA <- c()
scorr_vitD <- c()
scorr_vitE <- c()
scorr_zinc <- c()

#scorr_caff

#head(resDE.GXE.caff.sum)

MCls <- unique(resDE.GXE.caff.sum$MCls)
#MCls
comb.2 <- MCls
comb.2



conditions <- c("both", "sc", "GXE")
include <- "union"

conditions <- c("both")
include <- "both"

logFC <- 0



#as.numeric(cor.test(res1$z_score.y, res1$z_score.x, method = "spearman")$estimate)

#----
for(i in comb.2){
            res1 <- resDE.GXE.caff.sum %>% dplyr::filter(MCls==i)                                 %>% dplyr::filter(Significance.0 %in% conditions)
            #scorr <- cor(res1$z_score.y, res1$z_score.x, method = "spearman")
            scorr <- as.numeric(cor.test(res1$logFC, res1$estimate, method = "spearman")$estimate)
            scorr_caff <- c(scorr_caff, scorr)
        }


#----
for(i in comb.2){
            res1 <- resDE.GXE.nic.sum %>% dplyr::filter(MCls==i)                                 %>% dplyr::filter(Significance.0 %in% conditions)
            #scorr <- cor(res1$z_score.y, res1$z_score.x, method = "spearman")
            scorr <- as.numeric(cor.test(res1$logFC, res1$estimate, method = "spearman")$estimate)
            scorr_nic <- c(scorr_nic, scorr)
        }



#----
for(i in comb.2){
            res1 <- resDE.GXE.vitA.sum %>% dplyr::filter(MCls==i)                                 %>% dplyr::filter(Significance.0 %in% conditions)
            #scorr <- cor(res1$z_score.y, res1$z_score.x, method = "spearman")
            scorr <- as.numeric(cor.test(res1$logFC, res1$estimate, method = "spearman")$estimate)
            scorr_vitA <- c(scorr_vitA, scorr)
        }


#----
for(i in comb.2){
            res1 <- resDE.GXE.vitD.sum %>% dplyr::filter(MCls==i)                                 %>% dplyr::filter(Significance.0 %in% conditions)
            #scorr <- cor(res1$z_score.y, res1$z_score.x, method = "spearman")
            scorr <- as.numeric(cor.test(res1$logFC, res1$estimate, method = "spearman")$estimate)
            scorr_vitD <- c(scorr_vitD, scorr)
        }


#----
for(i in comb.2){
            res1 <- resDE.GXE.vitE.sum %>% dplyr::filter(MCls==i)                                 %>% dplyr::filter(Significance.0 %in% conditions)
            #scorr <- cor(res1$z_score.y, res1$z_score.x, method = "spearman")
            scorr <- as.numeric(cor.test(res1$logFC, res1$estimate, method = "spearman")$estimate)
            scorr_vitE <- c(scorr_vitE, scorr)
        }


#----
for(i in comb.2){
            res1 <- resDE.GXE.zinc.sum %>% dplyr::filter(MCls==i)                                 %>% dplyr::filter(Significance.0 %in% conditions)
            #scorr <- cor(res1$z_score.y, res1$z_score.x, method = "spearman")
            scorr <- as.numeric(cor.test(res1$logFC, res1$estimate, method = "spearman")$estimate)
            scorr_zinc <- c(scorr_zinc, scorr)
        }




#scorr_caff

scorr_data <- data.frame(MCls,
                         scorr_caff, scorr_nic,
                         scorr_vitA, scorr_vitD)
                        # scorr_vitE, scorr_zinc)

scorr_data <- data.frame(MCls,
                         scorr_caff, scorr_nic,
                         scorr_vitA, scorr_vitD,
                         scorr_vitE, scorr_zinc)


is.num <- sapply(scorr_data, is.numeric)
scorr_data[is.num] <- lapply(scorr_data[is.num], round, 2)

scorr_data


library(reshape2)
#sig.matrix <- as.matrix(scorr_data)
#sig.matrix

sig.melt <- melt(scorr_data)
sig.melt

sig.p <- ggplot(sig.melt, aes(x = variable, y = MCls)) +
          geom_tile(aes(fill=value), color="black")+
          scale_fill_gradient(low = "white", high = "white")+
          geom_text(aes(label = value), color = "black", size = 5) +
          ggtitle(paste0("logFC_", logFC, "_", include, "_logFC_corr_GXE_vs_sc")) +
          theme(axis.text.x = element_text(angle = 90),
                          axis.title.x = element_blank(),
                          axis.title.y = element_blank(),
                          legend.position="none")

library("ggpubr")

figure <- ggarrange(sig.p +
                    font("x.text", size = 14))
                   # ncol = 2, nrow = 2)

png(paste0("./1_GXE.outs/logFC_", logFC, "_", include, "_logFC_corr.png"), width=500, height=500, pointsize=16, res=125)
print(figure)
dev.off()







#----heatmap----
mat <- scorr_data
rownames(mat) <- mat$MCls
mat <- mat %>% dplyr::select(!MCls)
mat

mat2 <- as.matrix(mat)
mat2 <- mat2[nrow(mat2):1, ] 
mat2


fig1 <- pheatmap(mat2, border_color="NA",
                    cluster_rows=F, cluster_cols=F,
#                    annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend=T,
                    show_colnames=T,
                    show_rownames=T,
                    na_col="white",
                    fontsize_col=20,
                    fontsize_row=10)


figfn <- "./1_GXE.outs/logFC_corr_heatmap.png"
png(figfn, width=1800, height=2100,res=225)
print(fig1)
dev.off() 

















































