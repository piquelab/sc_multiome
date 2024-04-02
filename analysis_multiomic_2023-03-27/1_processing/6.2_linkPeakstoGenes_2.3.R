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


## #########################################################
## ############## list  as in 6.2_2.4.txt ##################
## #########################################################
##library("rlist")

list <- read.table("./comb_list_test.txt") %>% as.data.frame()
#%>% mutate(V1=gsub(".", "-", V1)
head(list)


bti <- list$V1
#gsub("\\.", "-", bti) 
#head(bti)

split <- str_split(bti, "\\.", simplify=T) %>% as.data.frame()
head(split)


MCls <- unique(split$V2)
#MCls <- unique(split[,2])
#MCls <- c("Bcell", "Platelet")
MCls

#results <- list()

columns <- c("celltype",
             "treatment",
             "individual1",
             "individual2",
             "cor",
             "pvalue",
             "z_cor",
             "z_pvalue")

df <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(df) <- columns
head(df)

for(i in 1:length(MCls)){
    split_1 <- split %>% dplyr::filter(V2==MCls[i])
    treats <- unique(split_1$V1)
    ##treats <- c("vitE", "water")
    ##treats
    for(j in 1:length(treats)){
        split_2 <- split_1 %>% dplyr::filter(V1==treats[j])
        ind <- unique(split_2$V3)
        ##ind <- c("AL_002", "AL_018")
        ##ind[1]
        for(k in 1:(length(ind)-1)){
            for(l in (k+1):length(ind)){
                tryCatch({    
                    dat1 <- readRDS(paste0("./6.2_linkPeakstoGenes_2.2.outs/links_", treats[j], ".", MCls[i], ".", ind[k], "_pval1_distdefault_score0.rds"))
                    dat2 <- readRDS(paste0("./6.2_linkPeakstoGenes_2.2.outs/links_", treats[j], ".", MCls[i], ".", ind[l], "_pval1_distdefault_score0.rds"))
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
                             treats[j],
                             ind[k],
                             ind[l],
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
                    cat(treats[j], MCls[i], ind[k], ind[l])
                })  
            }
        }
    }
}

#head(results)    
head(df)

write_rds(df, "./6.2_linkPeakstoGenes_2.3.outs/df_ind.rds")
            
    


## #########################################################
## ############## dat -  2.2.outs ##########################
## #########################################################
dat1 <- readRDS(paste0("./6.2_linkPeakstoGenes_2.2.outs/links_", "caffeine.NKcell.AL_045", "_pval1_distdefault_score0.rds"))
dat1 <- dat1 %>% mutate(peak.gene = paste0(peak, ".", gene))

dat2 <- readRDS(paste0("./6.2_linkPeakstoGenes_2.2.outs/links_", "nicotine.NKcell.AL_045", "_pval1_distdefault_score0.rds"))
dat2 <- dat2 %>% mutate(peak.gene = paste0(peak, ".", gene))


head(dat1)
head(dat2)

nrow(dat1 %>% dplyr::filter(seqnames=="chrY"))

nrow(dat1)
nrow(dat2)

#write.csv(head(dat1), paste0("./6.2_linkPeakstoGenes_2.3.outs/links_", "caffeine.NKcell.AL_045", "_pval1_distdefault_score0.csv"), row.names=FALSE)

#write.csv(head(dat2), paste0("./6.2_linkPeakstoGenes_2.3.outs/links_", "nicotine.NKcell.AL_045", "_pval1_distdefault_score0.csv"), row.names=FALSE)


dat <- full_join(dat1, dat2, by="peak.gene")
dat[is.na(dat)] <- 0

head(dat)
nrow(dat)

#write.csv(head(dat), paste0("./6.2_linkPeakstoGenes_2.3.outs/links_", "caffeine.nicotine.NKcell.AL_045", "_pval1_distdefault_score0.csv"), row.names=FALSE)


dat.2 <- dat %>% dplyr::filter(score.x != 0, score.y != 0)
head(dat.2)
nrow(dat.2)

write.csv(head(dat.2), paste0("./6.2_linkPeakstoGenes_2.3.outs/links_", "caffeine.nicotine.2.NKcell.AL_045", "_pval1_distdefault_score0.csv"), row.names=FALSE)


scorr <- as.numeric(cor.test(dat$score.x, dat$score.y, method = "spearman")$estimate)
scorr
pval <- as.numeric(cor.test(dat$score.x, dat$score.y, method = "spearman")$p.value)
pval


scorr <- as.numeric(cor.test(dat.2$score.x, dat.2$score.y, method = "spearman")$estimate)
scorr
pval <- as.numeric(cor.test(dat.2$score.x, dat.2$score.y, method = "spearman")$p.value)
pval


scorr_z <- as.numeric(cor.test(dat$zscore.x, dat$zscore.y, method = "spearman")$estimate)
scorr_z
pval_z <- as.numeric(cor.test(dat$zscore.x, dat$zscore.y, method = "spearman")$p.value)
pval_z



head(dat)
nrow(dat)

