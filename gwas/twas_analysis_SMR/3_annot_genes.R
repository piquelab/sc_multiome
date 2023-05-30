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


##library(ComplexHeatmap)
##library(circlize)
library(openxlsx)
library(cowplot)
library(ggrepel)
##


##options(scipen=200)

rm(list=ls())




####
#### annotation files 
fn <- "../../gtex_v8/torus/torus_input/zzz_torus.annot.gz"
annot <- fread(fn, header=T, data.table=F)


traits <- read.table("traits_ls.txt")$V1

for (trait in traits){

###    
outdir <- paste("./3_examples_output/", trait, "/", sep="")
if ( !file.exists(outdir)) dir.create(outdir, recursive=T)

    
### twas results
fn <- paste("./2_summary_twas/eQTL_conditions/", trait, "_twas.txt", sep="")
res <- read.table(fn, header=T)
res2 <- res%>%dplyr::select(Gene=gene, symbol, id_b38, pval_eqtl, PIP_eqtl=PIP, pval_gwas, FDR_twas=FDR)%>%
   arrange(FDR_twas) ### rank by FDR_twas
    
### colocalization FDR_intact<0.1   
fn <- paste("../enloc_analysis/INTACT_output/", trait, "_intact.txt", sep="")
intact <- read.table(fn, header=T)
intact2 <- intact%>%dplyr::filter(FDR<0.1)%>%
   dplyr::select(Gene, zscore_twas=zscore, FDR_intact=FDR, GRCP, GLCP, PCG)

### combine    
res_comb <- res2%>%inner_join(intact2, by="Gene")    


### annotation     
annot2 <- annot%>%filter(SNP%in%res_comb$id_b38)
x <- annot2[,-1]
rnz <- rowSums(x)
DF <- data.frame(id_b38=annot2$SNP, nannot=rnz)

    
### output all the genes FDR_intact<0.1
res_comb <- res_comb%>%left_join(DF, by="id_b38")
opfn <- paste(outdir, trait, "_intact_sig.xlsx", sep="")
write.xlsx(res_comb, opfn)


### output the genes with FDR_intact<0.1 and also annotated
res_comb2 <- res_comb%>%filter(nannot>1)
opfn2 <- paste(outdir, trait, "_intact_sig_annoted.xlsx", sep="")
write.xlsx(res_comb2, opfn2)

cat(trait, "total:", nrow(res_comb), "annotated:", nrow(res_comb2), "\n")
    
}



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


    



###
### summary traits
summ <- map_dfr(traits, function(trait){
   ##
   fn <- paste("./2_summary_twas/eQTL_conditions/", trait, "_twas.txt", sep="")
   gene1 <- read.table(fn, header=T)%>%dplyr::filter(FDR<0.1)%>%pull(gene)%>%unique()

   ##
   fn <- paste("./3_examples_output/", trait, "/", trait, "_intact_sig.xlsx", sep="") 
   x <- read.xlsx(fn)
   gene2 <- unique(x$Gene)

   gene3 <- x%>%dplyr::filter(nannot>1)%>%pull(Gene)%>%unique()

   df <- data.frame(traits=trait, ntwas=length(gene1), nPCG=length(gene2), n_annoted=length(gene3))
   df
})



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
