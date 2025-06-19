##
library(Matrix)
library(tidyverse)
##library(clusterProfiler)
library(data.table)
library(qvalue)

library(annotables)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

##
library(GenomicRanges)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(Signac)  ##, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(SeuratWrappers)

## library(cicero, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
## library(monocle3)


##library(ComplexHeatmap)
##library(circlize)
library(openxlsx)
library(cowplot)
library(ggrepel)

## library(plotgardener)
## library(S4Vectors)


## library(plotgardener, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
## install.packages("BiocManager", lib="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
## library(BiocManager, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
## BiocManager::install("plotgardener", lib="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
## BiocManager::install("S4Vectors", lib="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")

rm(list=ls())






###
### input data

fn <- "1_candidate_genes.xlsx"
gene_df <- read.xlsx(fn)
tmp <- gene_df%>%filter(`6_Monocyte_vitD_d`>0)
 
fn <- "traits_name_df.xlsx"
traits_df <- read.xlsx(fn)

fn <- "./vcf/gene_infor.txt"
gene_infor <- read.table(fn, header=F)
names(gene_infor) <- c("Gene", "symbol", "chr", "s0", "s1")


### fast-eqtl,
fn <- "../../eQTL_results/Whole_Blood.allpairs.txt.gz"
fast <- fread(fn, header=T, data.table=F)%>%mutate(ensgene=gsub("\\..*", "", gene_id))


### atac project
## fn <- "/nfs/rprdata/julong/sc_multiome/analysis_multiome_JW_2023-06-03/sc_multiome_data/1_processing/5.1_reCallPeak.outs/3_scATAC.annot.mitofilt.ndup.rds"
## sc <- read_rds(fn)
## sc2 <- subset(sc, subset=EXP%in%c("scGxE_1-1", "scGxE_1-2"))
## opfn <- "sc_multiome_small.rds"
## write_rds(sc2, file=opfn)

sc <- read_rds("sc_multiome_small.rds")

 
########################################################
### local zoom plots 
########################################################

 
### loop for traits 
for (k in 1:11){
###
trait <- traits_df$traits_long[k]
trait2 <- traits_df$traits_short[k]
    
outdir2 <- paste("./plots_outs/", trait2, "/", sep="")
if (!file.exists(outdir2)) dir.create(outdir2, showWarnings=F, recursive=T)
    
gene_df2 <- gene_df%>%filter(traits==trait2)
if ( nrow(gene_df2)==0 ) next


### input data, gwas
fn <- paste("../../gwas_data/", trait, "_impute_gwas.txt.gz", sep="")    
gwas <- fread(fn, header=T, data.table=F)


###
### loop for genes
ngene <- nrow(gene_df2)
for ( i in 1:ngene){

    
##dtss <- 2.5e+04 ## default 5e+04, 10e+04
 
dtss1 <- 5e+04
dtss2 <- 5e+04

###    
ens <- gene_df2$Gene[i]
symbol <- gene_df2$symbol[i]

### snps    
snp_i <- gene_df2$id_b38[i]    
snpinfor <- unlist(str_split(snp_i, "_"))    
chr_i <- snpinfor[1]
pos_i <- as.integer(snpinfor[2])    
pos_i1 <- pos_i-dtss1
pos_i2 <- pos_i+dtss2    

    
cat(trait2, i, ens, symbol, "\n")

    
###############################################    
### 0, calculate LD files
###############################################    
ind <- read.table("./vcf/sampleID.txt", header=F)$V1    

###
snpfn <- paste("./vcf/Example_genes_vcf/", ens, "_", symbol, "_",  chr_i, ".SNPid.txt.gz", sep="")    
snp_df <- fread(snpfn, header=F, data.table=F)

###
### vcf file

genfn <- paste("./vcf/Example_genes_vcf/", ens, "_", symbol, "_", chr_i, ".vcf.gz", sep="")
gen <- fread(genfn, header=T, data.table=F)
colnames(gen) <- gsub(".*\\[.*\\]", "", colnames(gen))

### snp_df-snp infor    
snp_df <- gen[,1:4]%>%mutate(snp_id=paste(CHROM, POS, REF, ALT, "b38", sep="_"))

### gen2-genotype dosage    
nind <- ncol(gen)    
gen2 <- gen[,5:nind]
rownames(gen2) <- snp_df$snp_id 

    
###    
### get ld infor  
snp_df2 <- snp_df%>%filter(POS>(pos_i1-1), POS<pos_i2+1, !snp_id%in%snp_i)
 
g0 <- as.numeric(gen2[snp_i,])    
g2 <- t(as.matrix(gen2[snp_df2$snp_id,]))
    
r2 <- t(cor(g0, g2)^2)[,1]    
snp_df2 <- snp_df2%>%
    mutate("r2"=r2[snp_id],
           gr=case_when(r2>=0.8~"1", r2>=0.6&r2<0.8~"2", r2>=0.4&r2<0.6~"3",
                        r2>=0.2&r2<0.4~"4", TRUE~"5"))
snp_df3 <- snp_df2%>%dplyr::select(id_b38=snp_id, r2, gr)    


    
##########################################
### plot 1, gwas 
##########################################


###    
gwas2 <- gwas%>%
   dplyr::filter(chr==chr_i, pos>pos_i1, pos<pos_i2)%>%    
   mutate(pos2=pos/1e+03, log10p=-log10(pval))
   
gwas2 <- gwas2%>%left_join(snp_df3, by="id_b38")%>%
     mutate(gr2=ifelse(is.na(r2), 6, gr))%>%mutate(gr2=ifelse(pos==pos_i, "0", gr2))
    
##
anno <- gwas2%>%filter(pos==pos_i)
anno$label_rs <- snp_i
    
p1 <- ggplot(gwas2, aes(x=pos2, y=log10p))+
   geom_point(aes(colour=factor(gr2), shape=factor(gr2), size=factor(gr2) ))+
   scale_size_manual("",
       values=c("0"=3.5, "1"=1.5,"2"=2, "3"=2, "4"=2, "5"=2, "6"=2), guide="none")+
   ## scale_color_manual("",
   ##     values=c("1"="#D43F3A", "2"="#EEA236", "3"="#5CB85C", "4"="#46B8DA",
   ##              "5"="#357EBD", "6"="#B8B8B8", "0"="#9632B8"), guide="none")+      
   scale_color_manual("",
       values=c("1"="#D43F3A", "2"="#EEA236", "3"="#5CB85C", "4"="#46B8DA",
                "5"="#357EBD", "6"="#B8B8B8", "0"="#9632B8"),
       breaks=c(1, 2, 3, 4, 5),
       labels=c("1"="0.8~1", "2"="0.6~0.8", "3"="0.4~0.6",
                "4"="0.2~0.4", "5"="0~0.2"),
       guide=guide_legend(override.aes=list(size=2)))+
   scale_shape_manual("",
       values=c("0"=18, "1"=19, "2"=19, "3"=19, "4"=19, "5"=19, "6"=19), guide="none")+    
   geom_text_repel(data=anno, aes(x=pos2, y=log10p, label=label_rs),
                   color="#9632B8", size=3.5, box.padding=0.8)+
   ## annotate("point", x=as.numeric(pos_i/1e+03), y=0, shape=18, size=3.5, color="blue")+ 
   xlab(bquote("position (kb)"~"("~.(chr_i)~")"))+
   ylab(bquote(~"GWAS-"~.(trait2)~-log[10]~"("~italic(p)~")"))+
   theme(plot.title=element_text(hjust=0.5, size=14),
      axis.text=element_text(size=10),
      axis.title=element_text(size=10),
      legend.position="none",
      legend.key=element_blank(),
      legend.key.size=grid::unit(0.4, "cm"),
      legend.background=element_blank(),
      legend.box.background=element_blank(),
      legend.text=element_text(size=8),
      plot.margin=unit(c(5.5, 12, 5.5, 5.5), "pt"),
       panel.background=element_blank(),
       panel.border=element_rect(fill=NA,color="black"),
       panel.grid.major=element_blank(),
       panel.grid.minor=element_blank())

## figfn <- paste(outdir2, "Figure", i, ".1_", ens, "_", symbol, "_gwas.png", sep="")       
## png(figfn, width=580, height=380, res=120)
## print(p1)
## dev.off()



#############################################
### plot2, FastQTL results 
#############################################

    
##    
fast2 <- fast%>% ## mutate(ensgene=gsub("\\..*", "", gene_id))%>%
    filter(ensgene==ens)%>%
    dplyr::select(gene_id, ensgene, id_b38=variant_id, pval=pval_nominal) 

    
x <- str_split(fast2$id_b38, "_", simplify=T)
fast2 <- fast2%>%mutate(chr=x[,1], pos=as.integer(x[,2]))%>%
   filter(pos>pos_i1, pos<pos_i2)%>%    
   mutate(pos2=pos/1e+03, log10p=-log10(pval))
     
ens2 <- fast2$gene_id[1]

###    
fast2 <- fast2%>%left_join(snp_df3, by="id_b38")%>%
    mutate(gr2=ifelse(is.na(r2), 6, gr), gr2=ifelse(pos==pos_i, 0, gr2))    

    
##
anno <- fast2%>%dplyr::filter(pos==pos_i)
anno$label_rs <- snp_i
 
p2 <- ggplot(fast2, aes(x=pos2, y=log10p))+
   geom_point(aes(colour=factor(gr2), size=factor(gr2),shape=factor(gr2)))+
   scale_size_manual("",
       values=c("0"=3.5, "1"=1.5,"2"=2, "3"=2, "4"=2, "5"=2, "6"=2), guide="none")+
   ## scale_color_manual("",
   ##     values=c("1"="#D43F3A", "2"="#EEA236", "3"="#5CB85C", "4"="#46B8DA",
   ##              "5"="#357EBD", "6"="#B8B8B8", "0"="#9632B8"), guide="none")+      
   scale_color_manual("",
       values=c("1"="#D43F3A", "2"="#EEA236", "3"="#5CB85C", "4"="#46B8DA",
                "5"="#357EBD", "6"="#B8B8B8", "0"="#9632B8"),
       breaks=c(1, 2, 3, 4, 5),
       labels=c("1"="0.8~1", "2"="0.6~0.8", "3"="0.4~0.6",
                "4"="0.2~0.4", "5"="0~0.2"),
       guide=guide_legend(override.aes=list(size=2)))+
   scale_shape_manual("",
       values=c("0"=18, "1"=19, "2"=19, "3"=19, "4"=19, "5"=19, "6"=19), guide="none")+        
   geom_text_repel(data=anno, aes(x=pos2, y=log10p, label=label_rs),
                   color="#9632B8", size=3.5, box.padding=0.8)+
   xlab(bquote("position (kb)"~"("~.(chr_i)~")"))+
   ylab(bquote(~"FastQTL-WBL"~log[10]~"("~italic(p)~")"))+
   ##ggtitle(bquote(~italic(.(symbol))~"expression"))+ 
   theme(plot.title=element_text(hjust=0.5, size=14),
      axis.text=element_text(size=10),
      axis.title=element_text(size=10),
      plot.margin=unit(c(5.5, 12, 5.5, 5.5), "pt"),
      legend.position="none",
       panel.background=element_blank(),
       panel.border=element_rect(fill=NA,color="black"),
       panel.grid.major=element_blank(),
       panel.grid.minor=element_blank())

## figfn <- paste(outdir2, "Figure", i, ".3_", ens, "_", symbol, "_eQTL.png", sep="")       
## png(figfn, width=450, height=380, res=120)
## print(p3)
## dev.off()


    

############################################
### plot 3, PIP 
#############################################
 
fn <- paste("./dap-g_outs/Whole_Blood/", ens2,  ".SNP.out", sep="")    
dap <- read.table(fn)%>%dplyr::select(id_b38=V2, PIP=V3, cluster=V5)

x <- str_split(dap$id_b38, "_", simplify=T)    
dap2 <- dap%>%mutate(chr=x[,1], pos=as.integer(x[,2]), pos2=pos/1e+03)%>%
   dplyr::filter(pos>pos_i1, pos<pos_i2)
     
###    
dap2 <- dap2%>%left_join(snp_df3, by="id_b38")%>%
    mutate(gr2=ifelse(is.na(r2), 6, gr))%>%mutate(gr2=ifelse(pos==pos_i, 0, gr2))    

    
##
anno <- dap2%>%dplyr::filter(pos==pos_i)
anno$label_rs <- snp_i
    
     
### plot  
p3 <- ggplot(dap2, aes(x=pos2, y=PIP))+
   geom_point(aes(colour=factor(gr2), size=factor(gr2),shape=factor(gr2)))+
   scale_size_manual("",
       values=c("0"=3, "1"=1.5,"2"=2, "3"=2, "4"=2, "5"=2, "6"=2), guide="none")+
   scale_color_manual("",
       values=c("1"="#D43F3A", "2"="#EEA236", "3"="#5CB85C", "4"="#46B8DA",
                "5"="#357EBD", "6"="#B8B8B8", "0"="#9632B8"),
       breaks=c(1, 2, 3, 4, 5),
       labels=c("1"="0.8~1", "2"="0.6~0.8", "3"="0.4~0.6",
                "4"="0.2~0.4", "5"="0~0.2"),
       guide=guide_legend(override.aes=list(size=2)))+
   scale_shape_manual("",
       values=c("0"=18, "1"=19, "2"=19, "3"=19, "4"=19, "5"=19, "6"=19), guide="none")+   
   geom_text_repel(data=anno, aes(x=pos2, y=PIP, label=label_rs),
                   color="#9632B8", size=3.5, box.padding=0.8)+
   xlab(bquote("position (kb)"~"("~.(chr_i)~")"))+
   ylab(bquote("Fine-map WBL ("~italic(PIP)~")"))+ylim(0,1)+
   ##ggtitle(bquote(~italic(.(symbol))))+ 
   theme(plot.title=element_text(hjust=0.5, size=14),
      axis.text=element_text(size=10),
      axis.title=element_text(size=10),
      plot.margin=unit(c(5.5, 12, 5.5, 5.5), "pt"),
      legend.position=c(0.85,0.8),
      legend.key=element_blank(),
      legend.key.size=grid::unit(0.4, "cm"),
      legend.background=element_blank(),
      legend.box.background=element_blank(),
      legend.text=element_text(size=8),
      panel.background=element_blank(),
       panel.border=element_rect(fill=NA,color="black"),
       panel.grid.major=element_blank(),
       panel.grid.minor=element_blank())



############################################
#### #4, gene track 
############################################
 
region2 <- paste(gsub("chr", "", chr_i), pos_i1, pos_i2, sep="-")     
p4 <- AnnotationPlot(sc, region=region2)&
   xlab(bquote("position (kb)"~"("~.(chr_i)~")"))&    
   theme(axis.title.x=element_text(size=10),
         axis.title.y=element_text(size=12),
         axis.text=element_blank(),       
         axis.ticks.x=element_blank())

ll <- length(p4$layers)

p4$layers[[ll]]$aes_params$size <- 0
p4$layers[[ll]]$aes_params$fontface <- "italic"    
## p4$layers[[5]]$aes_params$fonts <- "italic"    
     
###
comb <- plot_grid(p1, p2, p3, p4, align="hv", axis="lr", nrow=4, rel_heights=c(0.6, 0.6, 0.6, 0.4))
figfn <- paste0(outdir2, "Figure_", ens, "_", symbol, ".comb.png")
ggsave(figfn, comb, device="png", width=750, height=950, units="px", dpi=120)        


    
} ## END genes

} ## END traits
    






#################################################################
### three examples genes for final plots 
#################################################################


rm(list=ls())


trait <- traits_df$traits_long[k]
trait2 <- traits_df$traits_short[k]
    
outdir2 <- "./plots_outs/Example_genes/"
if (!file.exists(outdir2)) dir.create(outdir2, showWarnings=F, recursive=T)



fn <- "1.2_candidate_genes.xlsx"
gene_df <- read.xlsx(fn)
 
fn <- "traits_name_df.xlsx"
traits_df <- read.xlsx(fn)

fn <- "./vcf/gene_infor.txt"
gene_infor <- read.table(fn, header=F)
names(gene_infor) <- c("Gene", "symbol", "chr", "s0", "s1")


### fast-eqtl,
fn <- "../../eQTL_results/Whole_Blood.allpairs.txt.gz"
fast <- fread(fn, header=T, data.table=F)%>%mutate(ensgene=gsub("\\..*", "", gene_id))


### atac project
## fn <- "/nfs/rprdata/julong/sc_multiome/analysis_multiome_JW_2023-06-03/sc_multiome_data/1_processing/5.1_reCallPeak.outs/3_scATAC.annot.mitofilt.ndup.rds"
## sc <- read_rds(fn)
## sc2 <- subset(sc, subset=EXP%in%c("scGxE_1-1", "scGxE_1-2"))
## opfn <- "sc_multiome_small.rds"
## write_rds(sc2, file=opfn)

sc <- read_rds("sc_multiome_small.rds")


###
####
geneSel <- c("FOSL2", "IFNGR2", "IRF5")
traitSel <- c("Eczema", "Crohns", "RA")
gene_df2 <- gene_df%>%filter(symbol%in%geneSel, traits%in%traitSel)
 
###
### loop for genes
ngene <- nrow(gene_df2)
for ( i in 1:ngene){

    
trait2 <- gene_df2$traits[i]    
trait <- traits_df%>%filter(traits_short==trait2)%>%pull(traits_long)

    
### input data, gwas
fn <- paste("../../gwas_data/", trait, "_impute_gwas.txt.gz", sep="")    
gwas <- fread(fn, header=T, data.table=F)

    
##dtss <- 2.5e+04 ## default 5e+04, 10e+04
 
dtss1 <- 5e+04
dtss2 <- 5e+04

###    
ens <- gene_df2$Gene[i]
symbol <- gene_df2$symbol[i]

### snps    
snp_i <- gene_df2$id_b38[i]    
snpinfor <- unlist(str_split(snp_i, "_"))    
chr_i <- snpinfor[1]
pos_i <- as.integer(snpinfor[2])    
pos_i1 <- pos_i-dtss1
pos_i2 <- pos_i+dtss2    

snp_rs <- gene_df2$refsnp_id[i]
    
cat(trait2, i, ens, symbol, "\n")

    
###############################################    
### 0, calculate LD files
###############################################    
ind <- read.table("./vcf/sampleID.txt", header=F)$V1    

###
snpfn <- paste("./vcf/Example_genes_vcf/", ens, "_", symbol, "_",  chr_i, ".SNPid.txt.gz", sep="")    
snp_df <- fread(snpfn, header=F, data.table=F)

###
### vcf file

genfn <- paste("./vcf/Example_genes_vcf/", ens, "_", symbol, "_", chr_i, ".vcf.gz", sep="")
gen <- fread(genfn, header=T, data.table=F)
colnames(gen) <- gsub(".*\\[.*\\]", "", colnames(gen))

### snp_df-snp infor    
snp_df <- gen[,1:4]%>%mutate(snp_id=paste(CHROM, POS, REF, ALT, "b38", sep="_"))

### gen2-genotype dosage    
nind <- ncol(gen)    
gen2 <- gen[,5:nind]
rownames(gen2) <- snp_df$snp_id 

    
###    
### get ld infor  
snp_df2 <- snp_df%>%filter(POS>(pos_i1-1), POS<pos_i2+1, !snp_id%in%snp_i)
 
g0 <- as.numeric(gen2[snp_i,])    
g2 <- t(as.matrix(gen2[snp_df2$snp_id,]))
    
r2 <- t(cor(g0, g2)^2)[,1]    
snp_df2 <- snp_df2%>%
    mutate("r2"=r2[snp_id],
           gr=case_when(r2>=0.8~"1", r2>=0.6&r2<0.8~"2", r2>=0.4&r2<0.6~"3",
                        r2>=0.2&r2<0.4~"4", TRUE~"5"))
snp_df3 <- snp_df2%>%dplyr::select(id_b38=snp_id, r2, gr)    


    
##########################################
### plot 1, gwas 
##########################################


###    
gwas2 <- gwas%>%
   dplyr::filter(chr==chr_i, pos>pos_i1, pos<pos_i2)%>%    
   mutate(pos2=pos/1e+03, log10p=-log10(pval))
   
gwas2 <- gwas2%>%left_join(snp_df3, by="id_b38")%>%
     mutate(gr2=ifelse(is.na(r2), 6, gr))%>%mutate(gr2=ifelse(pos==pos_i, "0", gr2))
    
##
anno <- gwas2%>%filter(pos==pos_i)
anno$label_rs <- snp_rs
    
p1 <- ggplot(gwas2, aes(x=pos2, y=log10p))+
   geom_point(aes(colour=factor(gr2), shape=factor(gr2), size=factor(gr2) ))+
   scale_size_manual("",
       values=c("0"=3.5, "1"=1.5,"2"=2, "3"=2, "4"=2, "5"=2, "6"=2), guide="none")+
   ## scale_color_manual("",
   ##     values=c("1"="#D43F3A", "2"="#EEA236", "3"="#5CB85C", "4"="#46B8DA",
   ##              "5"="#357EBD", "6"="#B8B8B8", "0"="#9632B8"), guide="none")+      
   scale_color_manual("",
       values=c("1"="#D43F3A", "2"="#EEA236", "3"="#5CB85C", "4"="#46B8DA",
                "5"="#357EBD", "6"="#B8B8B8", "0"="#9632B8"),
       breaks=c(1, 2, 3, 4, 5),
       labels=c("1"="0.8~1", "2"="0.6~0.8", "3"="0.4~0.6",
                "4"="0.2~0.4", "5"="0~0.2"),
       guide=guide_legend(override.aes=list(size=2)))+
   scale_shape_manual("",
       values=c("0"=18, "1"=19, "2"=19, "3"=19, "4"=19, "5"=19, "6"=19), guide="none")+    
   geom_text_repel(data=anno, aes(x=pos2, y=log10p, label=label_rs),
                   color="#9632B8", size=3.5, box.padding=0.8,
                   min.segment.length=0, point.padding=0.1,
                   segment.curvature=1e-20, arrow=arrow(length=unit(0.015, "npc"))
                  )+
   ## annotate("point", x=as.numeric(pos_i/1e+03), y=0, shape=18, size=3.5, color="blue")+ 
   xlab(bquote("position (kb)"~"("~.(chr_i)~")"))+
   ylab(bquote(~"GWAS-"~-log[10]~"("~italic(p)~")"))+
   theme(plot.title=element_text(hjust=0.5, size=14),
      axis.text=element_text(size=10),
      axis.title=element_text(size=10),
      legend.position="none",
      legend.key=element_blank(),
      legend.key.size=grid::unit(0.4, "cm"),
      legend.background=element_blank(),
      legend.box.background=element_blank(),
      legend.text=element_text(size=8),
      plot.margin=unit(c(5.5, 12, 5.5, 5.5), "pt"),
       panel.background=element_blank(),
       panel.border=element_rect(fill=NA,color="black"),
       panel.grid.major=element_blank(),
       panel.grid.minor=element_blank())

## figfn <- paste(outdir2, "Figure", i, ".1_", ens, "_", symbol, "_gwas.png", sep="")       
## png(figfn, width=580, height=380, res=120)
## print(p1)
## dev.off()



#############################################
### plot2, FastQTL results 
#############################################

    
##    
fast2 <- fast%>% ## mutate(ensgene=gsub("\\..*", "", gene_id))%>%
    filter(ensgene==ens)%>%
    dplyr::select(gene_id, ensgene, id_b38=variant_id, pval=pval_nominal) 

    
x <- str_split(fast2$id_b38, "_", simplify=T)
fast2 <- fast2%>%mutate(chr=x[,1], pos=as.integer(x[,2]))%>%
   filter(pos>pos_i1, pos<pos_i2)%>%    
   mutate(pos2=pos/1e+03, log10p=-log10(pval))
     
ens2 <- fast2$gene_id[1]

###    
fast2 <- fast2%>%left_join(snp_df3, by="id_b38")%>%
    mutate(gr2=ifelse(is.na(r2), 6, gr), gr2=ifelse(pos==pos_i, 0, gr2))    

    
##
anno <- fast2%>%dplyr::filter(pos==pos_i)
anno$label_rs <- snp_rs
 
p2 <- ggplot(fast2, aes(x=pos2, y=log10p))+
   geom_point(aes(colour=factor(gr2), size=factor(gr2),shape=factor(gr2)))+
   scale_size_manual("",
       values=c("0"=3.5, "1"=1.5,"2"=2, "3"=2, "4"=2, "5"=2, "6"=2), guide="none")+
   ## scale_color_manual("",
   ##     values=c("1"="#D43F3A", "2"="#EEA236", "3"="#5CB85C", "4"="#46B8DA",
   ##              "5"="#357EBD", "6"="#B8B8B8", "0"="#9632B8"), guide="none")+      
   scale_color_manual("",
       values=c("1"="#D43F3A", "2"="#EEA236", "3"="#5CB85C", "4"="#46B8DA",
                "5"="#357EBD", "6"="#B8B8B8", "0"="#9632B8"),
       breaks=c(1, 2, 3, 4, 5),
       labels=c("1"="0.8~1", "2"="0.6~0.8", "3"="0.4~0.6",
                "4"="0.2~0.4", "5"="0~0.2"),
       guide=guide_legend(override.aes=list(size=2)))+
   scale_shape_manual("",
       values=c("0"=18, "1"=19, "2"=19, "3"=19, "4"=19, "5"=19, "6"=19), guide="none")+        
   geom_text_repel(data=anno, aes(x=pos2, y=log10p, label=label_rs),
                   color="#9632B8", size=3.5, box.padding=0.8,
                   min.segment.length=0, point.padding=0.1,
                   segment.curvature=1e-20, arrow=arrow(length=unit(0.015, "npc"))
                   )+
   xlab(bquote("position (kb)"~"("~.(chr_i)~")"))+
   ylab(bquote(~"FastQTL-WBL"~-log[10]~"("~italic(p)~")"))+
   ##ggtitle(bquote(~italic(.(symbol))~"expression"))+ 
   theme(plot.title=element_text(hjust=0.5, size=14),
      axis.text=element_text(size=10),
      axis.title=element_text(size=10),
      plot.margin=unit(c(5.5, 12, 5.5, 5.5), "pt"),
      legend.position="none",
       panel.background=element_blank(),
       panel.border=element_rect(fill=NA,color="black"),
       panel.grid.major=element_blank(),
       panel.grid.minor=element_blank())

## figfn <- paste(outdir2, "Figure", i, ".3_", ens, "_", symbol, "_eQTL.png", sep="")       
## png(figfn, width=450, height=380, res=120)
## print(p3)
## dev.off()


    

############################################
### plot 3, PIP 
#############################################
 
fn <- paste("./dap-g_outs/Whole_Blood/", ens2,  ".SNP.out", sep="")    
dap <- read.table(fn)%>%dplyr::select(id_b38=V2, PIP=V3, cluster=V5)

x <- str_split(dap$id_b38, "_", simplify=T)    
dap2 <- dap%>%mutate(chr=x[,1], pos=as.integer(x[,2]), pos2=pos/1e+03)%>%
   dplyr::filter(pos>pos_i1, pos<pos_i2)
     
###    
dap2 <- dap2%>%left_join(snp_df3, by="id_b38")%>%
    mutate(gr2=ifelse(is.na(r2), 6, gr))%>%mutate(gr2=ifelse(pos==pos_i, 0, gr2))    

    
##
anno <- dap2%>%dplyr::filter(pos==pos_i)
anno$label_rs <- snp_rs
    
     
### plot  
p3 <- ggplot(dap2, aes(x=pos2, y=PIP))+
   geom_point(aes(colour=factor(gr2), size=factor(gr2),shape=factor(gr2)))+
   scale_size_manual("",
       values=c("0"=3, "1"=1.5,"2"=2, "3"=2, "4"=2, "5"=2, "6"=2), guide="none")+
   scale_color_manual("",
       values=c("1"="#D43F3A", "2"="#EEA236", "3"="#5CB85C", "4"="#46B8DA",
                "5"="#357EBD", "6"="#B8B8B8", "0"="#9632B8"),
       breaks=c(1, 2, 3, 4, 5),
       labels=c("1"="0.8~1", "2"="0.6~0.8", "3"="0.4~0.6",
                "4"="0.2~0.4", "5"="0~0.2"),
       guide=guide_legend(override.aes=list(size=2)))+
   scale_shape_manual("",
       values=c("0"=18, "1"=19, "2"=19, "3"=19, "4"=19, "5"=19, "6"=19), guide="none")+   
   geom_text_repel(data=anno, aes(x=pos2, y=PIP, label=label_rs),
                   color="#9632B8", size=3.5, box.padding=0.8,
                   min.segment.length=0, point.padding=0.1,
                   segment.curvature=1e-20, arrow=arrow(length=unit(0.015, "npc"))
                   )+
   xlab(bquote("position (kb)"~"("~.(chr_i)~")"))+
   ylab(bquote("Fine-map WBL "~italic(PIP)))+ylim(0,1)+
   ##ggtitle(bquote(~italic(.(symbol))))+ 
   theme(plot.title=element_text(hjust=0.5, size=14),
      axis.text=element_text(size=10),
      axis.title=element_text(size=10),
      plot.margin=unit(c(5.5, 12, 5.5, 5.5), "pt"),
      legend.position=c(0.85,0.8),
      legend.key=element_blank(),
      legend.key.size=grid::unit(0.4, "cm"),
      legend.background=element_blank(),
      legend.box.background=element_blank(),
      legend.text=element_text(size=8),
      panel.background=element_blank(),
       panel.border=element_rect(fill=NA,color="black"),
       panel.grid.major=element_blank(),
       panel.grid.minor=element_blank())



############################################
#### #4, gene track 
############################################
 
region2 <- paste(gsub("chr", "", chr_i), pos_i1, pos_i2, sep="-")     
p4 <- AnnotationPlot(sc, region=region2)&
   xlab(bquote("position (kb)"~"("~.(chr_i)~")"))&    
   theme(axis.title.x=element_text(size=10),
         axis.title.y=element_text(size=12),
         axis.text=element_blank(),       
         axis.ticks.x=element_blank())

ll <- length(p4$layers)

p4$layers[[ll]]$aes_params$size <- 0
p4$layers[[ll]]$aes_params$fontface <- "italic"    
## p4$layers[[5]]$aes_params$fonts <- "italic"    
     
###
comb <- plot_grid(p1, p2, p3, p4, align="hv", axis="lr", nrow=4, rel_heights=c(0.6, 0.6, 0.6, 0.4))
figfn <- paste0(outdir2, "Figure_", trait2, "_", ens, "_", symbol, ".comb.png")
ggsave(figfn, comb, device="png", width=750, height=950, units="px", dpi=120)        

    
} ## END genes






     
######################################
### gene track  
######################################


fn <- "traits_name_df.xlsx"
traits_df <- read.xlsx(fn)


###
### gwas 
k <- 7
trait <- traits_df$traits_long[k]
trait2 <- traits_df$traits_short[k]

fn <- paste("../../gwas_data/", trait, "_impute_gwas.txt.gz", sep="")    
gwas <- fread(fn, header=T, data.table=F)



###
### snp information 
fn <- "SCAIP_final_bed.gz"
df_snp <- fread(fn, header=F, data.table=F)
names(df_snp) <- c("chr", "pos", "chr_pos_grch37", "snp_id", "rs", "chr_pos_grch38")


###
### motif snps
gen_fn <- "/nfs/rprdata/julong/sc_multiome/genetic_variant_annot/annot_jaspar2022/zzz_allmotif.bed.gz"
motif_snp <- fread(gen_fn, header=F, data.table=F)


###
### Response motif
motif_fn <- "/nfs/rprdata/julong/sc_multiome/genetic_variant_annot/4_SNPAnnot_correct/2.2_response_motif_th0.1.txt"
res_motif <- read.table(motif_fn, header=T)




fn <- "1.2_candidate_genes.xlsx"
gene_df <- read.xlsx(fn)


###
df0 <- gene_df%>%filter(symbol=="IRF5")


### fine-mapping 
### gene
ens <- df0$Gene[1] 
fns <- list.files("./dap-g_outs/Whole_Blood/", ".SNP.out")
fn0 <- fns[grepl(ens, fns)]

dap_fn <- paste("./dap-g_outs/Whole_Blood/", fn0, sep="")    
dap <- read.table(dap_fn)%>%dplyr::select(id_b38=V2, PIP=V3, cluster=V5)


###
### get snp rs
snp0 <- unlist(str_split(gsub("chr", "", dap$id_b38[1]), "_"))
snp0 <- paste(snp0[1], snp0[2], sep="_")

df_snp%>%filter(chr_pos_grch38==snp0)

 
## snps <- getBM(attributes=c('refsnp_id', 'chr_name', 'chrom_start', 'chrom_end'),
##     filters=c('chromosome_name', "start", "end"),
##      values=list(2, snp0[2], snp0[2]),
##      mart=ensembl)


gwas%>%filter(id_b38%in%dap$id_b38[1])
 
###
### motif infor
snp0 <- gsub("_b38$", "", df0$id_b38[1]) 
motif_snp2 <- motif_snp%>%filter(V4==snp0)
res_motif2 <- res_motif%>%filter(motif_ID%in%motif_snp2$V5)


###
###  differential peaks

fn <- "/nfs/rprdata/julong/sc_multiome/analysis_multiome_JW_2023-06-03/sc_multiome_data/2_Differential/1.2_DiffPeak.outs/2.0_DESeq.results_0.1_cn_sampleID+new_treat.rds"
res_DP <- read_rds(fn)

sigs <- res_DP%>%filter(p.adjusted<0.1, abs(estimate)>0.5)%>%pull(gene)%>%unique() 
 
pos_df <- str_split(unique(res_DP$gene), "-", simplify=T)%>%as.data.frame()
names(pos_df) <- c("chr", "pos1", "pos2")
pos_df <- pos_df%>%mutate(peak=paste(chr, pos1, pos2, sep="-"))

                          
snp0 <- unlist(str_split(dap$id_b38[1], "_"))
chr0 <- snp0[1]
pos0 <- as.numeric(snp0[2])

pos_df2 <- pos_df%>%filter(chr==chr0, as.numeric(pos1)<pos0, as.numeric(pos2)>pos0)
 
peakSel <- pos_df2$peak

res_DP2 <- res_DP%>%filter(gene%in%peakSel)%>%filter(p.adjusted<0.1)%>%arrange(desc(abs(estimate)))
