###
###
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


##library(ComplexHeatmap)
##library(circlize)
library(openxlsx)
library(cowplot)
library(ggrepel)
##

rm(list=ls())

###
###



######################
### autosome genes ###
######################

autosome <- as.character(1:22) 
grch38_unq <- grch38%>%
    dplyr::filter(chr%in%autosome, grepl("protein", biotype))%>%
    distinct(ensgene, chr, .keep_all=T)%>%dplyr::select(gene=ensgene, chr, biotype, symbol)


#########################
### read twas results ###
#########################

###

for (ii in c("eQTL_base", "eQTL_conditions", "eQTL_union")){

##ii <- "eQTL_conditions"
outdir2 <- paste0("./2_summary_twas/", ii, "/")
if ( !file.exists(outdir2)) dir.create(outdir2, showWarnings=F, recursive=T)


traits <- unique(read.table("traits_ls.txt")$V1)
for ( trait in traits){
    
fn <- paste0("./1_SMR_output/", ii, "/", trait, "_topPIP_twas.txt.gz") 
res <- fread(fn, header=T, data.table=F)

###    
res <- res%>%dplyr::select(gene=gene2, id_b38, pval_eqtl, PIP, pval_gwas)    
res <- res%>%dplyr::filter(gene%in%grch38_unq$gene)    
pval <- res$pval_gwas
fdr <- qvalue(pval)

## if ( class(fdr)=="try-error"){
##    fdr <- qvalue(pval, pi0=1) 
## }
    
res$FDR <- fdr$qvalues

###    
grch38_unq2 <- grch38_unq[,c("gene", "symbol")]     
res <- res%>%arrange(pval_gwas)%>%left_join(grch38_unq2, by="gene")

### output-1 all test genes    
opfn <- paste(outdir2, trait, "_twas.txt", sep="")
write.table(res, opfn, row.names=F, col.names=T, quote=F)

### output-2, twas genes with FDR<0.1    
res2 <- res%>%dplyr::filter(FDR<0.1)
opfn2 <- paste(outdir2, trait, "_twas2_sig.xlsx", sep="")
write.xlsx(res2, file=opfn2, overwrite=T)

cat(ii, trait, nrow(res2), "\n")    
} ### trait
} ### conditions


## fn <- "./2_summary_twas/eQTL_base/Asthma_twas.txt" 
## x1 <- read.table(fn, header=T)%>%dplyr::select(gene, x1=pval_gwas)

## fn <- "./analysis0_no/2_summary_twas/eQTL_base/Asthma_twas.txt"
## x2 <- read.table(fn, header=T)%>%dplyr::select(gene, x2=pval_gwas)

## plotDF <- x1%>%left_join(x2, by="gene")
 
## p <- ggplot(plotDF, aes(x=-log10(x1), y=-log10(x2)))+
##    geom_point()+
##    xlab(bquote(~log[10]~"pval"~" from impute"))+
##    ylab(bquote(~log[10]~"pval"~" from no-impute"))+ 
##    theme_bw()

## figfn <- "./2_summary_twas/Figue0_asthma.png"
## png(figfn, width=380, height=380, res=120)
## print(p)
## dev.off()


#########################################################
### histogram distribution of p-value across traits  ###
########################################################


col2 <- c("eQTL_base"="#868686", "eQTL_conditions"="#d01c8b", "eQTL_union"="#e66101")    

###
traits <- unique(read.table("traits_ls.txt")$V1)
ypos_vec <- c(800, 500, 2100, 1300, 2000)
names(ypos_vec) <- traits

###
for (trait in traits){

cat(trait, "\n")
    
### plot data    
dapdir <- c("eQTL_base", "eQTL_conditions", "eQTL_union")
plotDF <- map_dfr(dapdir, function(ii){
    ###
    fn <- paste("./2_summary_twas/", ii, "/", trait, "_twas.txt",  sep="")
    res <- read.table(fn, header=T) 
    res <- res%>%dplyr::select(gene, pval_gwas)
    res$dap <- ii
    res
})


#### annotation text    
annoDF <- map_dfr(dapdir, function(ii){
    ###
    fn <- paste("./2_summary_twas/", ii, "/", trait, "_twas.txt",  sep="")
    res <- read.table(fn, header=T)
    pval <- res$pval_gwas
    fdr <- qvalue(res$pval_gwas)
    ##pi0 <- propTrueNull(pval)   
    ngene <- nrow(res)
    
    pi0 <- round(fdr$pi0, digits=3)
    eq2 <- as.expression(bquote(~pi==.(pi0)))

    cat(ii, pi0, "\n")
    ##
    df2 <- tibble(dap=ii, pi=pi0, eq=eq2, yline=ngene*pi0)
    df2
})

ypos_i <- ypos_vec[trait]    
annoDF <- annoDF%>%mutate(xpos=0.75, ypos=ypos_i, yline2=yline/40)    

p2 <- ggplot(plotDF)+
   geom_histogram(aes(x=pval_gwas, color=factor(dap)), fill=NA, bins=40)+
   geom_text(data=annoDF, aes(x=xpos, y=ypos, label=eq), size=5, parse=T)+
   geom_hline(data=annoDF, aes(yintercept=yline2, color=factor(dap)), linetype="dashed")+ 
   xlab("P value of TWAS")+
   scale_y_continuous("Number of genes")+
   facet_wrap(~dap, ncol=3, scales="fixed")+
   scale_color_manual(values=col2)+ 
   theme_bw()+
   theme(## legend.title=element_blank(),
         ## legend.text=element_text(size=10),
         ## legend.key.size=grid::unit(0.8, "lines"),
         legend.position="none",
         axis.text=element_text(size=10.5),
         axis.title=element_text(size=12),
         strip.text=element_text(size=12))
     
 
figfn <- paste("./2_summary_twas/Figure1_",  trait, "_pval_hist.png", sep="")
png(figfn, width=850, height=350, res=120)
print(p2)
dev.off()    
 
}
    

###
###
trait <- traits[4]
twas_gene <- lapply(c("eQTL_base", "eQTL_conditions", "eQTL_union"), function(ii){
   ## 
   fn <- paste("./2_summary_twas/", ii, "/", trait, "_twas.txt",  sep="")
   res <- read.table(fn, header=T)
   res2 <- res%>%filter(FDR<0.1)
   res2$gene
})

lengths(twas_gene)

x12 <- intersect(twas_gene[[1]], twas_gene[[2]])
cat(length(x12), "\n")

x13 <- intersect(twas_gene[[1]], twas_gene[[3]])
cat(length(x13), "\n")

x23 <- intersect(twas_gene[[2]], twas_gene[[3]])
cat(length(x23), "\n")



##################################
#### Manhanttan plots for twas ###
##################################

traits <- sort(read.table("traits_ls.txt")$V1)

###
for (trait in traits){
###
### read twas data
fn <- paste("./2_summary_twas/eQTL_conditions/", trait, "_twas.txt", sep="")
res <- read.table(fn, header=T)
res <- res%>%mutate(is_sig=ifelse(FDR<0.1, 1, 0))

###
map_df <- as.data.frame(str_split(res$id_b38, "_", simplify=T))[,1:2]
names(map_df) <- c("chr", "pos")
res <- cbind(res, map_df)

### plot data
res <- res%>%mutate(chr=as.integer(gsub("chr", "", chr)), pos=as.numeric(pos))
res <- res%>%mutate(log10p=-log10(pval_gwas), gr2=factor(chr%%2))

## threshold line
sig <- res%>%filter(is_sig==1)%>%pull(log10p)%>%min()


### aes position 
ngene <- nrow(res)
nchr <- max(res$chr)
res$BPcum <- NA
s <- 0
for (i in 1:nchr){
   x <- res[res$chr==i, "pos"]+s
   res[res$chr==i,"BPcum"] <- x
   s <- max(x)
}

axis.set <- res%>%
   group_by(chr)%>%
   summarize(midpoint=(max(BPcum)+min(BPcum))/2)


###
### annotation data
fn <- paste("./3_examples_output/", trait, "/", trait, "_intact_sig.xlsx", sep="")
PCG <- read.xlsx(fn)%>%pull(Gene)

annoDF <- res%>%filter(gene%in%PCG) ## chr%in%c(2, 6, 17))    


###
### Manhattan plots
ymax <- max(res$log10p)+8

p2 <- ggplot(res) +
   geom_point(aes(x=BPcum, y=log10p, color=factor(gr2), size=log10p))+
   geom_text_repel(data=annoDF,
       aes(x=BPcum, y=log10p, label=symbol), box.padding=0.2,
            max.overlaps=20, fontface="italic", size=2.5)+      
   geom_hline(yintercept=sig, color="red", linetype="dashed")+ 
   scale_x_continuous("chromosome", label=axis.set$chr, breaks=axis.set$midpoint,
      limits=c(min(res$BPcum), max(res$BPcum)) )+
   scale_y_continuous(bquote(-log[10]~"("~italic(p)~")"), limits=c(0,ymax))+
   scale_color_manual(values=c("0"="#2171b5", "1"="#6baed6"))+
   scale_size_continuous(range=c(0.1,1.5)) +
   theme_bw()+
   theme(legend.position="none",
         axis.title=element_text(size=12),
         axis.text.x=element_text(size=10, angle=90, vjust=0.5, hjust=1),
         axis.text.y=element_text(size=12),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank())

###
figfn <- paste("./2_summary_twas/", "Figure2_", trait, "_manhattan.png", sep="")
png(figfn, width=800, height=500, pointsize=12, res=120)
print(p2)
dev.off()

cat(trait, "\n")    
}
