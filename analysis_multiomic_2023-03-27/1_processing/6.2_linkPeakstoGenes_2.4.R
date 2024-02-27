###
library(tidyverse)
## library(parallel)
library(data.table)
## library(purrr)
##library(GenomicRanges)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(Signac)
library(SeuratWrappers)
library(SeuratObject)
##library(EnsDb.Hsapiens.v75)
## library(cicero, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
## library(monocle3)
###
##library(EnsDb.Hsapiens.v75)
##library(EnsDb.Hsapiens.v86)
##library(BSgenome.Hsapiens.UCSC.hg38)
###
##library(ggplot2)
##library(cowplot)
##library(RColorBrewer)
##library(viridis)
##library(ggrastr)
##theme_set(theme_grey())


rm(list=ls())

#outdir <- "./3_Clustering.outs/"
#if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


##############################################
##  sys_args ##
##############################################
### passing argument
args=commandArgs(trailingOnly=T)

#condition <- args[1]

if (length(args)>0){
    condition <- args[1]
###dir <- args[2]
}else{
    condition <- "caffeine.Bcell.AL-002"
###dir <- "SCAIP_results"
}

condition <- gsub("\r", "", condition)

print(condition)








## ## #########################################################
## ## ############## list - ind  ##########################
## ## #########################################################
## ##library("rlist")

## ##list <- read.table("./comb_list_platelet.txt") %>% as.data.frame()
## list <- read.table(paste0("./", condition)) %>% as.data.frame()
## ##%>% mutate(V1=gsub(".", "-", V1)
## head(list)


## bti <- list$V1
## ##gsub("\\.", "-", bti) 
## ##head(bti)

## split <- str_split(bti, "\\.", simplify=T) %>% as.data.frame()
## ##head(split)


## MCls <- unique(split$V2)
## #MCls <- unique(split[,2])
## #MCls <- c("Bcell", "Platelet")
## ##MCls

## #results <- list()

## columns <- c("celltype",
##              "treatment",
##              "individual1",
##              "individual2",
##              "cor",
##              "pvalue",
##              "z_cor",
##              "z_pvalue")

## df <- data.frame(matrix(nrow = 0, ncol = length(columns)))
## colnames(df) <- columns
## ##head(df)

## for(i in 1:length(MCls)){
##     split_1 <- split %>% dplyr::filter(V2==MCls[i])
##     treats <- unique(split_1$V1)
##     ##treats <- c("vitE", "water")
##     ##treats
##     for(j in 1:length(treats)){
##         split_2 <- split_1 %>% dplyr::filter(V1==treats[j])
##         ind <- unique(split_2$V3)
##         ##ind <- c("AL_002", "AL_018")
##         ##ind[1]
##         for(k in 1:(length(ind)-1)){
##             for(l in (k+1):length(ind)){
##                 tryCatch({    
##                     dat1 <- readRDS(paste0("./6.2_linkPeakstoGenes_2.2.outs/links_", treats[j], ".", MCls[i], ".", ind[k], "_pval1_distdefault_score0.rds"))
##                     dat2 <- readRDS(paste0("./6.2_linkPeakstoGenes_2.2.outs/links_", treats[j], ".", MCls[i], ".", ind[l], "_pval1_distdefault_score0.rds"))
##                     dat1 <- dat1 %>% mutate(peak.gene = paste0(peak, ".", gene))
##                     dat2 <- dat2 %>% mutate(peak.gene = paste0(peak, ".", gene))
##                     dat <- full_join(dat1, dat2, by="peak.gene")
##                     dat[is.na(dat)] <- 0
##                     scorr <- as.numeric(cor.test(dat$score.x, dat$score.y, method = "spearman")$estimate)
##                     pval <- as.numeric(cor.test(dat$score.x, dat$score.y, method = "spearman")$p.value)
##                     scorr_z <- as.numeric(cor.test(dat$zscore.x, dat$zscore.y, method = "spearman")$estimate)
##                     pval_z <- as.numeric(cor.test(dat$zscore.x, dat$zscore.y, method = "spearman")$p.value)
##                     ## df1 <- data.frame(celltype = MCls[i],
##                     ##                  treatment = treats[j],
##                     ##                  individual1=ind[k],
##                     ##                  individual2=ind[l],
##                     ##                  cor=scorr,
##                     ##                  pvalue=pval,
##                     ##                  z_cor=scorr_z,
##                     ##                  z_pavlue=pval_z)
##                     df1 <- c(MCls[i],
##                              treats[j],
##                              ind[k],
##                              ind[l],
##                              scorr,
##                              pval,
##                              scorr_z,
##                              pval_z)
##                     df[nrow(df)+1, ] <- df1
##                     ## results <- list.append(results, op)
##                                         #print(scorr)
##                                         #print(pval)
##                                         #print(scorr_z)
##                                         #print(pval_z)
##                 }, error=function(e){
##                     ##message('An Error Occurred') 
##                     print('An Error Occurred') 
##                     ##print(e)
##                     cat(treats[j], MCls[i], ind[k], ind[l])
##                 })  
##             }
##         }
##     }
## }

## #head(results)    
## head(df)
## nrow(df)

## ##write_rds(df, paste0("./6.2_linkPeakstoGenes_2.4.outs/df_ind_platelet.rds")
## write_rds(df, paste0("./6.2_linkPeakstoGenes_2.4.outs/df_ind_", MCls[1], ".rds"))








## #########################################################
## ############## list - treat  ##########################
## #########################################################
##library("rlist")

##list <- read.table("./comb_list_platelet.txt") %>% as.data.frame()
list <- read.table(paste0("./", condition)) %>% as.data.frame()
##%>% mutate(V1=gsub(".", "-", V1)
head(list)


bti <- list$V1
##gsub("\\.", "-", bti) 
##head(bti)

split <- str_split(bti, "\\.", simplify=T) %>% as.data.frame()
##head(split)


MCls <- unique(split$V2)
#MCls <- unique(split[,2])
#MCls <- c("Bcell", "Platelet")
##MCls

#results <- list()

columns <- c("celltype",
             "individual",
             "treatment1",
             "treatment2",
             "cor",
             "pvalue",
             "z_cor",
             "z_pvalue")

df <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(df) <- columns
##head(df)

for(i in 1:length(MCls)){
    split_1 <- split %>% dplyr::filter(V2==MCls[i])
    ind <- unique(split_1$V3)
    ##ind <- c("AL_002", "AL_018")
    ##ind[1]
    for(j in 1:length(ind)){
        split_2 <- split_1 %>% dplyr::filter(V3==ind[j])
        treats <- unique(split_2$V1)
        ##treats <- c("vitE", "water")
        ##treats
        for(k in 1:(length(treats)-1)){
            for(l in (k+1):length(treats)){
                tryCatch({    
                    dat1 <- readRDS(paste0("./6.2_linkPeakstoGenes_2.2.outs/links_", treats[k], ".", MCls[i], ".", ind[j], "_pval1_distdefault_score0.rds"))
                    dat2 <- readRDS(paste0("./6.2_linkPeakstoGenes_2.2.outs/links_", treats[l], ".", MCls[i], ".", ind[j], "_pval1_distdefault_score0.rds"))
                    dat1 <- dat1 %>% mutate(peak.gene = paste0(peak, ".", gene))
                    dat2 <- dat2 %>% mutate(peak.gene = paste0(peak, ".", gene))
                    dat <- full_join(dat1, dat2, by="peak.gene")
                    dat[is.na(dat)] <- 0
                    scorr <- as.numeric(cor.test(dat$score.x, dat$score.y, method = "spearman")$estimate)
                    pval <- as.numeric(cor.test(dat$score.x, dat$score.y, method = "spearman")$p.value)
                    scorr_z <- as.numeric(cor.test(dat$zscore.x, dat$zscore.y, method = "spearman")$estimate)
                    pval_z <- as.numeric(cor.test(dat$zscore.x, dat$zscore.y, method = "spearman")$p.value)
                    ## df1 <- data.frame(celltype = MCls[i],
                    ##                  treatment = treats[j],
                    ##                  individual1=ind[k],
                    ##                  individual2=ind[l],
                    ##                  cor=scorr,
                    ##                  pvalue=pval,
                    ##                  z_cor=scorr_z,
                    ##                  z_pavlue=pval_z)
                    df1 <- c(MCls[i],
                             ind[j],
                             treats[k],
                             treats[l],
                             scorr,
                             pval,
                             scorr_z,
                             pval_z)
                    df[nrow(df)+1, ] <- df1
                    ## results <- list.append(results, op)
                                        #print(scorr)
                                        #print(pval)
                                        #print(scorr_z)
                                        #print(pval_z)
                }, error=function(e){
                    ##message('An Error Occurred') 
                    print('An Error Occurred') 
                    ##print(e)
                    cat(ind[j], MCls[i], treats[k], treats[l])
                })  
            }
        }
    }
}

#head(results)    
head(df)
nrow(df)

##write_rds(df, paste0("./6.2_linkPeakstoGenes_2.4.outs/df_treats_platelet.rds")
write_rds(df, paste0("./6.2_linkPeakstoGenes_2.4.outs/df_treats_", MCls[1], ".rds"))
