####################################################################
# gene_enrichment_analysis with ClusterProfiler gseGO
####################################################################
# https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/
# BiocManager::install("clusterProfiler", version = "3.8")
# BiocManager::install("clusterProfiler")
# BiocManager::install("pathview")
# BiocManager::install("enrichplot")
library(clusterProfiler)
library(enrichplot)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)
#library(Seurat)
#library(SeuratDisk)
#library(SeuratData)
#library(Signac)
library(ggnewscale)
library(ggridges)
library(europepmc)
library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(SeuratObject)
library(Signac)
library(SeuratWrappers)
library(SeuratData)
library(qqman)
library(JASPAR2020)
library(JASPAR2022)
library(TFBSTools)
## library(JASPAR2022, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
## library(TFBSTools, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
library(BSgenome.Hsapiens.1000genomes.hs37d5, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(ChIPseeker, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
###
library(ggplot2)
library(cowplot) ##lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(grid)
library(gridExtra)
library(ggExtra)
library(RColorBrewer)
library(ggsci)
library(viridis)
library(ComplexHeatmap)
library(circlize)

## install.packages("ggnewscale")
## install.packages("ggridges")
## install.packages("europepmc")

# SET THE DESIRED ORGANISM HERE
# http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
#organism = "org.Dm.eg.db"


# use without quotes in the function, but with quotes in library
#organism = get("org.Hs.eg.db")
organism = org.Hs.eg.db 
#organism = EnsDb.Hsapiens.v86

BiocManager::install(organism, character.only = TRUE)

library("org.Hs.eg.db", character.only = TRUE)


###################################################
# reading in data from deseq2
###################################################
# df = read.csv("drosphila_example_de.csv", header=TRUE)
#df = caff_up

#----resDE----
resDE <-  read_rds("../2_Differential/1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn.rds") %>% as.data.frame()
head(resDE)

class(resDE$gene)

resDE$gene <- as.factor(resDE$gene)

resDE.caff <- resDE %>% dplyr::filter(contrast=="caffeine")
resDE.nic <- resDE %>% dplyr::filter(contrast=="nicotine")
resDE.vitA <- resDE %>% dplyr::filter(contrast=="vitA")
resDE.vitD <- resDE %>% dplyr::filter(contrast=="vitD")
resDE.vitE <- resDE %>% dplyr::filter(contrast=="vitE")
resDE.zinc <- resDE %>% dplyr::filter(contrast=="zinc")


sig_resDE <- resDE %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)
head(sig_resDE)

#df
#df[,1]
#rownames(df)
#colnames(df)

sig_resDE.caff <- sig_resDE %>% dplyr::filter(contrast=="caffeine")
sig_resDE.nic <- sig_resDE %>% dplyr::filter(contrast=="nicotine")
sig_resDE.vitA <- sig_resDE %>% dplyr::filter(contrast=="vitA")
sig_resDE.vitD <- sig_resDE %>% dplyr::filter(contrast=="vitD")
sig_resDE.vitE <- sig_resDE %>% dplyr::filter(contrast=="vitE")
sig_resDE.zinc <- sig_resDE %>% dplyr::filter(contrast=="zinc")




#----gene data----
# we want the log2 fold change
#original_gene_list <- df$log2FoldChange
original_gene_list <- sig_resDE$estimate
head(original_gene_list)

# name the vector
# names(original_gene_list) <- df$X
names(original_gene_list) <- sig_resDE$gene
head(original_gene_list)

# omit any NA values
gene_list<-na.omit(original_gene_list)
head(gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)
head(gene_list)

#keytypes(org.Hs.eg.db)
#keytypes(EnsDb.Hsapiens.v86)


#----gseGO function----
gse <- gseGO(geneList=gene_list,
                          ont ="BP",
                          keyType = "SYMBOL",
                          nPerm = 10000,
                          minGSSize = 1,
                          maxGSSize = 800,
                          pvalueCutoff = 0.1,
                          verbose = TRUE,
                          OrgDb = organism,
                          pAdjustMethod = "none")

# gse <- gseGO(geneList=gene_list,
#              ont ="ALL",
#              keyType = "ENSEMBL",
#              pvalueCutoff = 0.1,
#              verbose = TRUE,
#              OrgDb = organism,
#              pAdjustMethod = "BH")

gse

require(DOSE)

#dotplot(gse)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)






#----function_loop----
require(DOSE)

# for all treatments
#head(resDE)
head(sig_resDE)

treats <- unique(sig_resDE$contrast)

treats <- "nicotine"

treats

celltypes <- unique(sig_resDE$MCls)
celltypes

#keytypes(org.Hs.eg.db)
#keytypes(EnsDb.Hsapiens.v86)

for (i in dev.list()[1]:dev.list()[length(dev.list())]) {
       dev.off()
       }

dev.list()


for(i in treats){
    for (j in celltypes){
        tryCatch({
        #resDE2 <- resDE%>%dplyr::filter(contrast==i)
        #resDE2 <- sig_resDE%>%dplyr::filter(contrast==i)
        resDE2 <- sig_resDE%>%dplyr::filter(contrast==i, MCls==j)
        #print(head(resDE2))
        original_gene_list <- resDE2$estimate
        #head(original_gene_list)
        names(original_gene_list) <- resDE2$gene
        #head(original_gene_list)
        gene_list<-na.omit(original_gene_list)
        gene_list = sort(gene_list, decreasing = TRUE)
        #head(gene_list)
        #----gseGO function----
        gse <- gseGO(geneList=gene_list,
                     ont ="BP",
                     keyType = "SYMBOL",
                     nPerm = 10000,
                     minGSSize = 1,
                     maxGSSize = 800,
                     pvalueCutoff = 0.1,
                     verbose = TRUE,
                     OrgDb = organism,
                     pAdjustMethod = "none")
                    # pAdjustMethod = "BH")
        # gse <- gseGO(geneList=gene_list,
        #              ont ="ALL",
        #              keyType = "ENSEMBL",
        #              pvalueCutoff = 0.1,
        #              verbose = TRUE,
        #              OrgDb = organism,
        #              pAdjustMethod = "BH")
        #gse
        #dotplot(gse)
        ##
        ##
        #png(paste0("./1_GOenrichment.outs/GO/sig_gseGO_", i, ".png"), width=1800, height=2100,res=125)
        #p <- dotplot(gse, showCategory=10, title = paste0("Gene_expression_", i), split=".sign") + facet_grid(.~.sign)
        ##
        ##
        png(paste0("./1_GOenrichment.outs/GO/all_combinations/sig_gseGO_", i, "_", j, ".png"), width=1800, height=2100,res=125)
        p <- dotplot(gse, showCategory=10, title = paste0("Gene_expression_", i, "_", j), split=".sign") + facet_grid(.~.sign)
        print(p)
        dev.off()
        },  error=function(e){cat("ERROR : no term enriched under specific pvalueCutoff for",i,"_",j, "\n")})
    }}







## #----other plots----
## emapplot(gse, showCategory = 10)

## # categorySize can be either 'pvalue' or 'geneNum'
## cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 10)


## ridgeplot(gse) + labs(x = "enrichment distribution")

## # Use the `Gene Set` param for the index in the title, and as the value for geneSetId
## gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)


## terms <- gse$Description[1:3]
## pmcplot(terms, 2010:2018, proportion=FALSE)

















####################################################################
# gene_enrichment_analysis with ClusterProfiler enrichGO function
####################################################################
enrichGO(gene, OrgDb, keytype = "ENTREZID", ont = "MF", pvalueCutoff = 0.05, pAdjustMethod = "BH", universe, qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = FALSE)













######################################################################
###### KEGG Gene Set Enrichment Analysis with ClusterProfiler
######################################################################
## # Convert gene IDs for gseKEGG function
## # We will lose some genes here because not all IDs will be converted
## resDE2 <- sig_resDE %>% dplyr::filter(contrast == "caffeine")
## head(resDE2)
## nrow(resDE2)

## original_gene_list <- resDE2$estimate
## #head(original_gene_list)

## names(original_gene_list) <- resDE2$gene
## #head(original_gene_list)

## #gene_list<-na.omit(original_gene_list)
## #gene_list = sort(gene_list, decreasing = TRUE)
    
## original_gene_list
## length(original_gene_list)

## #ids<-bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=organism)
## ids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
## head(ids)
## nrow(ids)


## # remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
## #dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]
## dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]
## colnames(dedup_ids)[1] <- "gene"
## head(dedup_ids)
## nrow(dedup_ids)


## # Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
## #df2 = df[df$X %in% dedup_ids$SYMBOL,]
## df2 <- resDE2 %>% dplyr::filter(gene %in% dedup_ids$gene)
## head(df2)
## nrow(df2)


## # Create a new column in df2 with the corresponding ENTREZ IDs
## #df2$Y = dedup_ids$ENTREZID
## df2 <- left_join(df2, dedup_ids, by="gene")
## head(df2)
## nrow(df2)

## # Create a vector of the gene unuiverse
## #kegg_gene_list <- df2$log2FoldChange
## kegg_gene_list <- df2$estimate
## head(kegg_gene_list)

## # Name vector with ENTREZ ids
## #names(kegg_gene_list) <- df2$Y
## names(kegg_gene_list) <- df2$ENTREZID
## head(kegg_gene_list)

## # omit any NA values
## kegg_gene_list<-na.omit(kegg_gene_list)
## head(kegg_gene_list)

## # sort the list in decreasing order (required for clusterProfiler)
## kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
## head(kegg_gene_list)

## kegg_organism = "hsa"

## kk2 <- gseKEGG(geneList     = kegg_gene_list,
##                organism     = kegg_organism,
##                nPerm        = 10000,
##                minGSSize    = 1,
##                maxGSSize    = 800,
##                pvalueCutoff = 0.1,
##                pAdjustMethod = "none",
##                keyType       = "ncbi-geneid")

## png(paste0("./1_GOenrichment.outs/KEGG/sig_gseKEGG_", "caffeine", ".png"), width=1800, height=2100,res=125)
## p <- dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
## print(p)
## dev.off()



#----function_loop----
head(sig_resDE)

treats <- unique(sig_resDE$contrast)
treats

celltypes <- unique(sig_resDE$MCls)
celltypes

kegg_organism = "hsa"


for (i in dev.list()[1]:dev.list()[length(dev.list())]) {
       dev.off()
       }

dev.list()


for(i in treats){
    for (j in celltypes){
        tryCatch({
        #resDE2 <- resDE%>%dplyr::filter(contrast==i)
        resDE2 <- sig_resDE%>%dplyr::filter(contrast==i)
        #resDE2 <- sig_resDE%>%dplyr::filter(contrast==i, MCls==j)
        #print(head(resDE2))
        original_gene_list <- resDE2$estimate
        #head(original_gene_list)
        names(original_gene_list) <- resDE2$gene
        #head(original_gene_list)
        ids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
        dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]
        colnames(dedup_ids)[1] <- "gene"
        df2 <- resDE2 %>% dplyr::filter(gene %in% dedup_ids$gene)
        df2 <- left_join(df2, dedup_ids, by="gene")
        kegg_gene_list <- df2$estimate
        names(kegg_gene_list) <- df2$ENTREZID
        kegg_gene_list<-na.omit(kegg_gene_list)
        kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
        #----gseKEGG function----
        kk2 <- gseKEGG(geneList     = kegg_gene_list,
                       organism     = kegg_organism,
                       nPerm        = 10000,
                       minGSSize    = 1,
                       maxGSSize    = 800,
                       pvalueCutoff = 0.1,
                       pAdjustMethod = "none",
                       keyType       = "ncbi-geneid")
        ##
        ##
        #png(paste0("./1_GOenrichment.outs/KEGG/sig_gseKEGG_", i, ".png"), width=1800, height=2100,res=125)
        #p <- dotplot(kk2, showCategory=10, title = paste0("KEGG_Gene_expression_", i), split=".sign") + facet_grid(.~.sign)
        ##
        ##
        png(paste0("./1_GOenrichment.outs/KEGG/all_combinations/sig_gseKEGG_", i, "_", j,  ".png"), width=1800, height=2100,res=125)
        p <- dotplot(kk2, showCategory=10, title = paste0("KEGG_Gene_expression_", i, "_", j), split=".sign") + facet_grid(.~.sign)
        print(p)
        dev.off()
        },  error=function(e){cat("ERROR : no term enriched under specific pvalueCutoff for",i,"_",j, "\n")})
    }
}

























#################################################################################
#################################################################################
# reading in data for resDP
#################################################################################
#################################################################################
# df = read.csv("drosphila_example_de.csv", header=TRUE)
#df = caff_up

#----resDP----
resDP <-  read_rds("../2_Differential/1.2_DiffPeak.outs/2.0_DESeq.results_0.1_cn.rds") %>% as.data.frame()

resDP <- resDP %>% dplyr::rename("peak"="gene") 

head(resDP)





#----peakAnno----
peakAnno <- read_rds("../2_Differential/2.2_compareRNAandATAC.outs/1_annot.signac_macs2_0.1_cn.rds")
x1 <- peakAnno
head(x1)

x1.2 <- x1 %>% dplyr::rename("peak"="query_region")%>% mutate(peak=gsub(":", "-", peak)) %>% dplyr::select(peak, gene_id, gene_name, dtss=distance) %>% mutate(dtss=abs(dtss))
head(x1.2)
nrow(x1.2)





#----resDP.join----
resDP.join <- resDP%>%left_join(x1.2, by="peak")
head(resDP.join)
nrow(resDP.join)

colnames(resDP.join)

resDP.join.filt <- resDP.join %>% dplyr::filter(dtss <= 100000)
head(resDP.join.filt)
nrow(resDP.join.filt)



resDP <- resDP.join.filt %>% dplyr::rename("gene"="gene_name")
head(resDP)

resDP.caff <- resDP %>% dplyr::filter(contrast=="caffeine")
resDP.nic <- resDP %>% dplyr::filter(contrast=="nicotine")
resDP.vitA <- resDP %>% dplyr::filter(contrast=="vitA")
resDP.vitD <- resDP %>% dplyr::filter(contrast=="vitD")
resDP.vitE <- resDP %>% dplyr::filter(contrast=="vitE")
resDP.zinc <- resDP %>% dplyr::filter(contrast=="zinc")


sig_resDP <- resDP %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)
head(sig_resDP)

#df
#df[,1]
#rownames(df)
#colnames(df)

sig_resDP.caff <- sig_resDP %>% dplyr::filter(contrast=="caffeine")
sig_resDP.nic <- sig_resDP %>% dplyr::filter(contrast=="nicotine")
sig_resDP.vitA <- sig_resDP %>% dplyr::filter(contrast=="vitA")
sig_resDP.vitD <- sig_resDP %>% dplyr::filter(contrast=="vitD")
sig_resDP.vitE <- sig_resDP %>% dplyr::filter(contrast=="vitE")
sig_resDP.zinc <- sig_resDP %>% dplyr::filter(contrast=="zinc")




#----function_loop----
require(DOSE)

# for all treatments
#head(resDP)

comb <- unique(resDP$contrast)

comb <- unique(sig_resDP$contrast)

comb

for(i in comb){
#    resDP2 <- resDP%>%dplyr::filter(contrast==i)
    resDP2 <- sig_resDP%>%dplyr::filter(contrast==i)
    original_gene_list <- resDP2$estimate
    #head(original_gene_list)
    names(original_gene_list) <- resDP2$gene
    #head(original_gene_list)
    gene_list<-na.omit(original_gene_list)
    gene_list = sort(gene_list, decreasing = TRUE)
    #head(gene_list)
    #----gseGO function----
    gse <- gseGO(geneList=gene_list,
                          ont ="BP",
                          keyType = "SYMBOL",
                          nPerm = 10000,
                          minGSSize = 1,
                          maxGSSize = 800,
                          pvalueCutoff = 0.1,
                          verbose = TRUE,
                          OrgDb = organism,
                          pAdjustMethod = "none")
#                          pAdjustMethod = "BH")
    # gse <- gseGO(geneList=gene_list,
    #              ont ="ALL",
    #              keyType = "ENSEMBL",
    #              pvalueCutoff = 0.1,
    #              verbose = TRUE,
    #              OrgDb = organism,
    #              pAdjustMethod = "BH")
    #gse
    #dotplot(gse)
#    png(paste0("./1_GOenrichment.outs/resDP_gseGO_BH_", i, ".png"), width=1800, height=2100,res=125)
    png(paste0("./1_GOenrichment.outs/resDP_sig_gseGO_", i, ".png"), width=1800, height=2100,res=125)
    p <- dotplot(gse, showCategory=10, title = paste0("Accessibility_", i), split=".sign") + facet_grid(.~.sign)
    print(p)
    dev.off()
    }










































#################################################################################
#################################################################################
# reading in data for motif
#################################################################################
#################################################################################
# df = read.csv("drosphila_example_de.csv", header=TRUE)
#df = caff_up

#----motif----
df <- read_rds("../3_motif/1_motif.outs/3_motif.enrich.direction_macs2_0.1_cn.rds")
head(df)

#----resDP----
resDP <-  read_rds("../2_Differential/1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn.rds")
head(resDP)

resDP.caff <- resDP %>% dplyr::filter(contrast=="caffeine")
resDP.nic <- resDP %>% dplyr::filter(contrast=="nicotine")
resDP.vitA <- resDP %>% dplyr::filter(contrast=="vitA")
resDP.vitD <- resDP %>% dplyr::filter(contrast=="vitD")
resDP.vitE <- resDP %>% dplyr::filter(contrast=="vitE")
resDP.zinc <- resDP %>% dplyr::filter(contrast=="zinc")




sig_resDP <- resDP %>% drop_na(p.adjusted) %>% dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)
head(sig_resDP)

#df
#df[,1]
#rownames(df)
#colnames(df)

sig_resDP.caff <- sig_resDP %>% dplyr::filter(contrast=="caffeine")
sig_resDP.nic <- sig_resDP %>% dplyr::filter(contrast=="nicotine")
sig_resDP.vitA <- sig_resDP %>% dplyr::filter(contrast=="vitA")
sig_resDP.vitD <- sig_resDP %>% dplyr::filter(contrast=="vitD")
sig_resDP.vitE <- sig_resDP %>% dplyr::filter(contrast=="vitE")
sig_resDP.zinc <- sig_resDP %>% dplyr::filter(contrast=="zinc")













