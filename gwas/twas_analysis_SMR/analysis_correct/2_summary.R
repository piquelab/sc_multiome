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


outdir <- "./2_summary_twas/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)



######################
### autosome genes ###
######################

autosome <- as.character(1:22) 
grch38_unq <- grch38%>%
    dplyr::filter(chr%in%autosome, grepl("protein", biotype))%>%
    distinct(ensgene, chr, .keep_all=T)%>%dplyr::select(gene=ensgene, chr, biotype, symbol)
 

 
## traits <- read.table("traits_of_interest.txt")$V1
## traits <- sort(traits)

## traitDF <- data.frame(traits_long=traits,
##     traits_short=c("CAD", "Crohns", "IBD", "UC", "MS", "SLE", "RA", "HTN",
##                    "Asthma", "Eczema", "Psoriasis", "Height"))

## opfn <- "traits_name_df.xlsx"
## write.xlsx(traitDF, opfn, overwrite=T)


######################################################
### Extract protein coding genes and calcualte FDR ###
######################################################


### use all gene to  calculate FDR then extract the protein coding genes


traits <- read.table("../traits_of_interest.txt")$V1
traits <- sort(traits)


###
### loop by approach to selecting lead variants
for (ii in c("eQTL_fastQTL", "eQTL_dap_base", "eQTL_dap_conditions")){

##ii <- "eQTL_conditions"
outdir2 <- paste0("./2_summary_twas/", ii, "/")
if ( !file.exists(outdir2)) dir.create(outdir2, showWarnings=F, recursive=T)

### loop by  traits 
for ( trait in traits){

if ( grepl("dap", ii)){    
   fn <- paste0("./1_SMR_output/", ii, "/", trait, "_topPIP_twas.txt.gz") 
}else{
   fn <- paste0("./1_SMR_output/", ii, "/", trait, "_minP_twas.txt.gz")
}   
res <- fread(fn, header=T, data.table=F)

###    
res <- res%>%mutate(gene2=gsub("\\..*", "", gene), FDR=qvalue(pval_gwas)$qvalues)
res <- res%>%dplyr::filter(gene2%in%grch38_unq$gene)    

## pval <- res$pval_gwas
## fdr <- qvalue(pval)
## res$FDR <- fdr$qvalues

###    
grch38_unq2 <- grch38_unq[,c("gene", "symbol")]     
res <- res%>%arrange(pval_gwas)%>%left_join(grch38_unq2, by=c("gene2"="gene"))
    
 
### output-1 all test genes    
opfn <- paste(outdir2, trait, "_twas.txt", sep="")
write.table(res, opfn, row.names=F, col.names=T, quote=F)

cat(ii, trait, sum(res$FDR<0.1), "\n")    
} ### trait
} ### conditions




###
### summary number of significant genes

traitDF <- read.xlsx("../traits_name_df.xlsx")
nt <- nrow(traitDF)

conditions <- c("eQTL_fastQTL", "eQTL_dap_base", "eQTL_dap_conditions")
summ <- map_dfr(1:nt, function(k){
    trait <- traitDF$traits_long[k]
    trait2 <- traitDF$traits_short[k]
    ##
    ngene <- sapply(conditions, function(ii){
       ## 
       fn <- paste0("./2_summary_twas/", ii, "/", trait, "_twas.txt", sep="")
       x <- read.table(fn, fill=T, header=T)
       sigs <- x%>%filter(FDR<0.05)%>%pull(gene2)%>%unique()
       sigs <- unique(gsub("\\..*", "", sigs)) 
       nn <- length(sigs)
       nn
     })
    ##
    df2 <- data.frame(traits=trait, trait_short=trait2)
    df2 <- cbind(df2, t(ngene))
    df2
})    

opfn <- "./2_summary_twas/1_summary_ngene_FDR0.05.xlsx"
write.xlsx(summ, opfn, overwrite=T)
 


###
### proportion of p<0.01
traitDF <- read.xlsx("../traits_name_df.xlsx")
nt <- nrow(traitDF)

conditions <- c("eQTL_fastQTL", "eQTL_dap_base", "eQTL_dap_conditions")
summ <- map_dfr(1:nt, function(i){
    trait <- traitDF$traits_long[i]
    trait2 <- traitDF$traits_short[i]
    ##
    ngene <- sapply(conditions, function(ii){
       ##
       if ( grepl("dap", ii)){    
           fn <- paste0("./1_SMR_output/", ii, "/", trait, "_topPIP_twas.txt.gz") 
       }else{
           fn <- paste0("./1_SMR_output/", ii, "/", trait, "_minP_twas.txt.gz")
       }   
       x <- fread(fn, header=T, data.table=F)
       prop <- round(mean(x$pval_gwas<0.01), digits=3)
       prop 
    })
    ##
    df2 <- data.frame(traits=trait, traits_short=trait2)
    df2 <- cbind(df2, t(ngene))
    df2
})    

opfn <- "./2_summary_twas/1.2_summary_prop_sig0.01.xlsx"
write.xlsx(summ, opfn, overwrite=T)

## opfn <- "./2_summary_twas/1.3_summary_prop_sig0.01.xlsx"
## write.xlsx(summ, opfn, overwrite=T)







#########################################################
### histogram distribution of p-value across traits  ###
########################################################


outdir2 <- "./2_summary_twas/twa/"
if ( !file.exists(outdir2)) dir.create(outdir2, showWarnings=F, recursive=T)


col2 <- c("fastQTL"="#868686", "dap_base"="#e66101", "dap_conditions"="#d01c8b")    
dap_value <- c("fastQTL"=1, "dap_base"=2, "dap_conditions"=3)
facet_lab <- as_labeller(c("fastQTL"="fastQTL", "dap_base"="dap_base", "dap_conditions"="dap_scMultiome"))

###
traits <- read.table("traits_of_interest.txt")$V1
traits <- sort(traits)
 
## ypos_vec <- c(800, 500, 2100, 1300, 2000)
## names(ypos_vec) <- traits
 
###
for (trait in traits){

cat(trait, "\n")
    
### plot data    
dapdir <- c("eQTL_fastQTL", "eQTL_dap_base", "eQTL_dap_conditions")
plotDF <- map_dfr(dapdir, function(ii){
   ###
   if ( grepl("dap", ii)){    
      fn <- paste0("./1_SMR_output/", ii, "/", trait, "_topPIP_twas.txt.gz") 
   }else{
      fn <- paste0("./1_SMR_output/", ii, "/", trait, "_minP_twas.txt.gz")
   }   
   res <- read.table(fn, header=T)
   res <- res%>%dplyr::select(gene, pval_gwas)
   res$dap <- gsub("eQTL_", "", ii)
   res
})

plotDF <- plotDF%>%mutate(dap_val=as.numeric(dap_value[dap]), dap2=fct_reorder(dap, dap_val))
    
        
#### annotation text    
annoDF <- map_dfr(dapdir, function(ii){
   ###
   if ( grepl("dap", ii)){    
      fn <- paste0("./1_SMR_output/", ii, "/", trait, "_topPIP_twas.txt.gz") 
   }else{
      fn <- paste0("./1_SMR_output/", ii, "/", trait, "_minP_twas.txt.gz")
   }
   res <- read.table(fn, header=T) ##%>%dplyr::select(gene, pval_gwas)
   pval <- res$pval_gwas
   fdr <- qvalue(res$pval_gwas)
    ##pi0 <- propTrueNull(pval)   
    ngene <- nrow(res)

    ## FDR
    pi0 <- round(fdr$pi0, digits=3) ## from FDR
    eq2 <- as.expression(bquote(~pi==.(pi0)))

    ## lambda
    ## pi0 <- (2*sum(pval<=0.5))/length(pval)
    ## pi0 <- round(pi0, digits=3)
    ## eq2 <- as.expression(bquote(~pi<=.(pi0)))
    
    cat(ii, pi0, "\n")
    ii2 <- gsub("eQTL_", "", ii)
    ##
    df2 <- tibble(dap=ii2, pi=pi0, eq=eq2, yline=ngene*pi0)
    df2
})
annoDF <- annoDF%>%mutate(dap_val=as.numeric(dap_value[dap]), dap2=fct_reorder(dap, dap_val))
    
    
ypos_i <- 800   
annoDF <- annoDF%>%mutate(xpos=0.75, ypos=ypos_i, yline2=yline/40)    
    
p2 <- ggplot(plotDF)+
   geom_histogram(aes(x=pval_gwas, color=factor(dap)), fill=NA, bins=40)+
   geom_text(data=annoDF, aes(x=xpos, y=ypos, label=eq), size=5, parse=T)+
   geom_hline(data=annoDF, aes(yintercept=yline2, color=factor(dap)), linetype="dashed")+ 
   xlab("P value of TWAS")+
   scale_y_continuous("Number of genes", expand=expansion(mult=c(0, 0.3)))+
   facet_wrap(~dap2, ncol=3, scales="fixed", labeller=facet_lab)+
   scale_color_manual(values=col2)+
   ggtitle(trait)+ 
   theme_bw()+
   theme(## legend.title=element_blank(),
         ## legend.text=element_text(size=10),
         ## legend.key.size=grid::unit(0.8, "lines"),
         legend.position="none",
         axis.text=element_text(size=10),
         axis.title=element_text(size=12),
         strip.text=element_text(size=10),
         plot.title=element_text(hjust=0.5, size=12))
     
 
figfn <- paste(outdir2, "Figure1_",  trait, "_pval_hist.png", sep="")
ggsave(figfn, p2, device="png", width=850, height=380, units="px", dpi=120)
## png(figfn, width=850, height=350, res=120)
## print(p2)
## dev.off()     
}
    

###
###
## trait <- traits[4]
## twas_gene <- lapply(c("eQTL_base", "eQTL_conditions", "eQTL_union"), function(ii){
##    ## 
##    fn <- paste("./2_summary_twas/", ii, "/", trait, "_twas.txt",  sep="")
##    res <- read.table(fn, header=T)
##    res2 <- res%>%filter(FDR<0.1)
##    res2$gene
## })

## lengths(twas_gene)

## x12 <- intersect(twas_gene[[1]], twas_gene[[2]])
## cat(length(x12), "\n")

## x13 <- intersect(twas_gene[[1]], twas_gene[[3]])
## cat(length(x13), "\n")

## x23 <- intersect(twas_gene[[2]], twas_gene[[3]])
## cat(length(x23), "\n")



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
