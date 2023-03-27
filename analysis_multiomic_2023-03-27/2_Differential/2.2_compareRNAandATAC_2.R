###
library(Matrix)
library(tidyverse)
library(qqman)
library(qvalue)
##
library(DESeq2)
library(biobroom)
library(ashr)
library(GenomicRanges)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(Signac)
### annotation required package
library(clusterProfiler)
library(org.Hs.eg.db)
library(annotables)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(EnsDb.Hsapiens.v75)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(EnsDb.Hsapiens.v86)
#BiocManager::install("ChIPseeker")
library(ChIPseeker)
library(ChIPseeker, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
###
library(ggplot2)
library(cowplot)
#library(cowplot, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(grid)
library(gridExtra)
library(ggExtra)
library(RColorBrewer)
library(ggsci)
library(RColorBrewer)
library(viridis)
theme_set(theme_grey())


outdir <- "./2.2_compareRNAandATAC_2.outs/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)


#----resDP----
resDP <- read_rds("../2_Differential/1.2_DiffPeak.outs/2.0_DESeq.results_0.1_cn.rds") %>% as.data.frame()
resDP <- resDP %>% mutate(MCls=gsub("_", "-", MCls)) %>% mutate(comb=paste(MCls, contrast, sep="_"))%>% dplyr::rename("peak"="gene") %>% drop_na(p.value)
#head(resDP)
#nrow(resDP)

sig_resDP <- resDP %>% dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)
#head(sig_resDP)
#nrow(sig_resDP)

peakAll <- unique(sig_resDP$peak)
#head(peakAll)
#length(peakAll)



########################
#----peakAnno----
########################
peakAnno <- read_rds("../2_Differential/2.2_compareRNAandATAC.outs/1_annot.signac_macs2_0.1_cn.rds")
x1 <- peakAnno

head(x1)
nrow(x1)
length(unique(x1$query_region))

x1.2 <- x1 %>% dplyr::rename("peak"="query_region")%>% mutate(peak=gsub(":", "-", peak)) %>% dplyr::select(peak, gene_id, gene_name, dtss=distance) %>% mutate(dtss=abs(dtss))

head(x1.2)
nrow(x1.2)

nrow(x1.2 %>% dplyr::filter(dtss<100000))

length(unique(x1.2$peak))


peakAnno2 <- x1.2%>%dplyr::filter(peak%in%peakAll)
head(peakAnno2)
nrow(peakAnno2)

###########################
##resDP.join
###########################
resDP.join <- resDP%>%left_join(x1.2, by="peak")

head(resDP.join)
nrow(resDP.join)

resDP.join.filt <- resDP.join %>% dplyr::filter(dtss <= 100000)
#head(resDP.join.filt)
nrow(resDP.join.filt)

sig_resDP.join <- resDP.join %>% dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)
#head(sig_resDP.join)
#table(sig_resDP.join$dtss)

nrow(sig_resDP.join)
#unique(sig_resDP.join$comb)

sig_resDP.join.filt <- sig_resDP.join %>% dplyr::filter(dtss <= 100000)
nrow(sig_resDP.join.filt)





#----links----
dat <- read_rds("../1_processing/6.1_linkPeakstoGenes.outs/links_dist100000.rds")
head(dat)
nrow(dat)
max(dat$pvalue)


dat2 <- dat%>%dplyr::filter(peak%in%peakAll)
head(dat2)
nrow(dat2)

length(unique(dat2$peak))


resDP.join.links <- resDP%>%left_join(dat, by="peak")
head(resDP.join.links)
nrow(resDP.join.links)

resDP.join.links.2 <- resDP.join.links %>% dplyr::filter(!is.na(gene))
head(resDP.join.links.2)
nrow(resDP.join.links.2)

link.peaks <- unique(resDP.join.links.2$peak)
head(link.peaks)
length(link.peaks)

sig_resDP.join.links <- resDP.join.links.2 %>% dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5) %>% mutate(comb.gene=paste0(comb, ".", gene))

head(sig_resDP.join.links)
#table(sig_resDP.join$dtss)
nrow(sig_resDP.join.links)
#unique(sig_resDP.join$comb)

sig.link.peaks <- unique(sig_resDP.join.links$peak)
head(sig.link.peaks)
length(sig.link.peaks)








## #----TRIAL----
## sig_resDP.join.links.CD4.caff <- sig_resDP.join.links %>% dplyr::filter(MCls=="0-CD4Naive", contrast=="caffeine")
## nrow(sig_resDP.join.links.CD4.caff)
## length(unique(sig_resDP.join.links.CD4.caff$peak))

## peaks.trial <- unique(sig_resDP.join.links.CD4.caff$peak)

## links.trial <- map_dfr(peaks.trial, function(ii){
##     if(nrow(sig_resDP.join.links.CD4.caff %>% dplyr::filter(peak==ii))==1){
##         res <- sig_resDP.join.links.CD4.caff %>% dplyr::filter(peak==ii)
##         } else{
##             res <- NULL
##             }
##     res
##     })

## head(links.trial)
## nrow(links.trial)


## sig_resDP.join.filt.2.CD4.caff <- sig_resDP.join.filt.2 %>% dplyr::filter(MCls=="0-CD4Naive", contrast=="caffeine")
## nrow(sig_resDP.join.filt.2.CD4.caff)
## length(unique(sig_resDP.join.filt.2.CD4.caff$peak))

## filt.2.trial <- sig_resDP.join.filt.2.CD4.caff %>% dplyr::filter(peak %in% links.trial$peak)
## head(filt.2.trial)
## nrow(filt.2.trial)

## trial.join <- left_join(links.trial, filt.2.trial, by="peak") %>% dplyr::select(peak, seqnames, start, end, gene.x, gene.y, score, zscore, comb.gene.y)

## head(trial.join)
## nrow(trial.join)



## head(dat)










##----filtering sig_resDP.join.filt based on link.peaks----
sig_resDP.join.filt.2 <- sig_resDP.join.filt %>% dplyr::filter(peak %in% link.peaks) %>% mutate(gene=gene_name, comb.gene=paste0(comb, ".", gene_name)) 
head(sig_resDP.join.filt.2)
nrow(sig_resDP.join.filt.2)



length(unique(sig_resDP.join.links$peak))
length(unique(sig_resDP.join.filt.2$peak))


links.filt2 <- left_join(sig_resDP.join.links, sig_resDP.join.filt.2, by="gene")


#----resDE----
## resDE <- read_rds("../2_Differential/1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn.rds")%>%as.data.frame()
## resDE <- resDE %>% drop_na(p.value) %>% mutate(MCls=gsub("_", "-", MCls)) %>% mutate(comb=paste(MCls, contrast, sep="_"))
## #head(resDE)
## #nrow(resDE)

## sig_resDE <- resDE %>% dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)
## #head(sig_resDE)
## #nrow(sig_resDE)
## #unique(sig_resDE$comb)


resDE <- read_rds("../2_Differential/3_treatEffect_2_2.outs/resDE.mitofilt.rds")%>%as.data.frame()

resDE <- resDE %>% drop_na(p.value) %>% mutate(MCls=gsub("_", "-", MCls)) %>% mutate(comb=paste(MCls, contrast, sep="_"))
#head(resDE)
#nrow(resDE)

sig_resDE <- resDE %>% dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)
head(sig_resDE)
nrow(sig_resDE)
#unique(sig_resDE$comb)






## #----individual steps DE and DP combined----
## DPG <- unique(resDP.join.filt$gene_name)
## length(DPG)

## DEG <- unique(resDE$gene)
## length(DEG)

## comb <- resDP.join.filt %>% dplyr::filter(gene_name%in%DEG)
## nrow(comb)

## comb <- comb%>%mutate(gene=gene_name)

## comb <- comb%>%left_join(resDE, by="gene")

## head(comb)
## nrow(comb)


## max(comb$estimate.y)
## min(comb$estimate.y)
## max(comb$estimate.x)
## min(comb$estimate.x)

## png("./2.2_compareRNAandATAC_2.outs/Figure1.2_all_scatter__macs2_0.1_cn.png")
## p <- ggplot(comb, aes(x=estimate.y, y=estimate.x, color=comb.x)) +
##         geom_point()+
##         xlab("genes") + ylab("filtered_peaks")+
##     #    xlim(-4, 10) + ylim(-4, 10)+
##         ggtitle("fold.enrichment")
## print(p)
## dev.off()



#----FUNCTION----
comb <- sort(unique(resDE$comb))
comb

geneList <- NULL
geneList <- lapply(comb, function(ii){
    resDP2 <- sig_resDP%>%dplyr::filter(comb==ii)
    resDE2 <- sig_resDE%>%dplyr::filter(comb==ii)
    DP <- unique(resDP2$peak)
    DEG <- unique(resDE2$gene)
    x <- peakAnno2%>%dplyr::filter(peak%in%DP, gene_name%in%DEG)
    if(nrow(x)==0){
        res1 <- NA
    }else{
        res1 <- unique(x$gene_name)
    }
    res <- list(DEG=res1)
})

geneList
length(geneList)


## geneList <- NULL
## geneList  <- mapdfr(lapply(comb, function(ii){
##     resDP2 <- sig_resDP.join%>%dplyr::filter(comb==ii)
##     resDE2 <- sig_resDE%>%dplyr::filter(comb==ii)
##     #DP <- unique(resDP2$peak)
##     DEG <- unique(resDE2$gene)
##     #x <- peakAnno2%>%dplyr::filter(peak%in%DP, gene_name%in%DEG)
##     x <- resDP2 %>% dplyr::filter(gene_name%in%DEG)
##     if(nrow(x)==0){
##         res1 <- NA
##     }else{
##         res1 <- unique(x$gene_name)
##         res1 <- x
##         }
##     res <- list(DEG=res1)
## }))
## geneList
## length(geneList)

## gene1 <- lapply(geneList, function(x) x[[1]])
## gene1

## gene1 <- unique(unlist(gene1[!is.na(gene1)]))
## gene1

DEG <- unique(resDE$gene)
sig_DEG <- unique(sig_resDE$gene)
#DEG
length(gene1)
length(DEG)
length(sig_DEG)


#----FUNCTION directly using sig object----
head(sig_resDP.join.filt)
head(sig_resDE)

comb <- sort(unique(resDE$comb))
comb

head(sig_resDP.join.filt.2)
head(sig_resDP.join.links)


sig_comb.filt <- map_dfr(comb, function(ii){
        resDP2 <- sig_resDP.join.filt%>%dplyr::filter(comb==ii)
        resDE2 <- sig_resDE%>%dplyr::filter(comb==ii)
        #DP <- unique(resDP2$peak)
        DEG <- unique(resDE2$gene)
        #x <- peakAnno2%>%dplyr::filter(peak%in%DP, gene_name%in%DEG)
        x <- resDP2 %>% dplyr::filter(gene_name%in%DEG)
        if(nrow(x)==0){
            res1 <- NULL
            #next
        }else{
            #res1 <- unique(x$gene_name)
            #res1 <- x
            res1 <- x%>%mutate(gene=gene_name)
            res1 <- res1%>%left_join(resDE2, by="gene")
        }
        #res <- list(DEG=res1)
        res1
})

head(sig_comb.filt)
nrow(sig_comb.filt)
length(unique(sig_comb.filt$gene))
table(sig_comb.filt$comb.x)
max(sig_comb.filt$estimate.y)
min(sig_comb.filt$estimate.y)
max(sig_comb.filt$estimate.x)
min(sig_comb.filt$estimate.x)



sig_comb.filt.2 <- map_dfr(comb, function(ii){
        resDP2 <- sig_resDP.join.filt.2%>%dplyr::filter(comb==ii)
        resDE2 <- sig_resDE%>%dplyr::filter(comb==ii)
        #DP <- unique(resDP2$peak)
        DEG <- unique(resDE2$gene)
        #x <- peakAnno2%>%dplyr::filter(peak%in%DP, gene_name%in%DEG)
        x <- resDP2 %>% dplyr::filter(gene_name%in%DEG)
        if(nrow(x)==0){
            res1 <- NULL
            #next
        }else{
            #res1 <- unique(x$gene_name)
            #res1 <- x
            res1 <- x%>%mutate(gene=gene_name)
            res1 <- res1%>%left_join(resDE2, by="gene")
        }
        #res <- list(DEG=res1)
        res1
})


head(sig_resDE)


sig_comb.links <- map_dfr(comb, function(ii){
        resDP2 <- sig_resDP.join.links%>%dplyr::filter(comb==ii)
        resDE2 <- sig_resDE%>%dplyr::filter(comb==ii)
        ##DP <- unique(resDP2$peak)
        DEG <- unique(resDE2$gene)
        ##x <- peakAnno2%>%dplyr::filter(peak%in%DP, gene_name%in%DEG)
        x <- resDP2 %>% dplyr::filter(gene%in%DEG)
        if(nrow(x)==0){
            res1 <- NULL
            #next
        }else{
            #res1 <- unique(x$gene_name)
            #res1 <- x
            #res1 <- x%>%mutate(gene=gene_name)
            res1 <- x%>%left_join(resDE2, by="gene")
        }
        #res <- list(DEG=res1)
        res1
})

#install.packages("ggpmisc")
library(ggpmisc)

typeofDP <- "links"
comb.2 <- unique(sig_comb.links$contrast.x)
comb.2

for(i in comb.2){
    res1 <- sig_comb.links %>% dplyr::filter(contrast.x==i)
    png(paste0("./2.2_compareRNAandATAC_2.outs/Figure1.1_significant_scatter_macs2_0.1_cn_", typeofDP, "_", i, ".png"))
    p <- ggplot(res1, aes(x=estimate.y, y=estimate.x, color=contrast.x)) +
        geom_point()+
        xlab("sig_genes") + ylab("sig_peaks")+
        xlim(-4, 10) + ylim(-4, 5)+
        stat_poly_line() +
        stat_poly_eq() +
        ggtitle("fold.enrichment")
    print(p)
    dev.off()
}




## comb.2 <- unique(sig_comb$MCls.x)
## comb.2
## for(i in comb.2){
##     res1 <- sig_comb %>% dplyr::filter(MCls.x==i)
##     png(paste0("./2.2_compareRNAandATAC_2.outs/Figure1.1_significant_scatter_macs2_0.1_cn_", i, ".png"))
##     p <- ggplot(res1, aes(x=estimate.y, y=estimate.x, color=MCls.x)) +
##     geom_point()+
##     xlab("sig_genes") + ylab("sig_peaks")+
##     xlim(-4, 10) + ylim(-4, 5)+
##     ggtitle("fold.enrichment")
##     print(p)
##     dev.off()
## }


png("./2.2_compareRNAandATAC_2.outs/Figure1.1_significant_scatter_macs2_0.1_cn.png")
p <- ggplot(sig_comb, aes(x=estimate.y, y=estimate.x, color=comb.x)) +
        geom_point()+
        xlab("sig_genes") + ylab("sig_peaks")+
        xlim(-4, 10) + ylim(-4, 5)+
        ggtitle("fold.enrichment")
print(p)
dev.off()



## sig_resDE_sig_DP.join <- map_dfr(comb, function(ii){
##     resDP2 <- sig_resDP.join%>%dplyr::filter(comb==ii)
##     resDE2 <- sig_resDE%>%dplyr::filter(comb==ii)
##     #DP <- unique(resDP2$peak)
##     DEG <- unique(resDP2$gene_name)
##     #x <- peakAnno2%>%dplyr::filter(peak%in%DP, gene_name%in%DEG)
##     x <- resDE2 %>% dplyr::filter(gene%in%DEG)
##     if(nrow(x)==0){
##         res1 <- NULL
##     }else{
##         #res1 <- unique(x$gene_name)
##         res1 <- x
##         }
##     #res <- list(DEG=res1)
##     res1
## })
## head(sig_resDE_sig_DP.join)
## nrow(sig_resDE_sig_DP.join)
## length(unique(sig_resDE_sig_DP.join$gene))



#----compute the number of differential peaks nearby DEG or DVG genes and distance to DPs of DEG or DVG----
dd <- map_dfr(comb, function(ii){
            resDP2 <- sig_resDP%>%dplyr::filter(comb==ii)
                #    head(resDP2)
                    resDE2 <- sig_resDE%>%dplyr::filter(comb==ii)
                #    head(resDE2)
                    DP <- unique(resDP2$peak)
                #    head(DP)
                    DEG <- unique(resDE2$gene)
                #    head(DEG)
                    d1 <- peakAnno2%>%dplyr::filter(peak%in%DP, gene_name%in%DEG)
                #    head(d1)
                    d1 <- d1%>%group_by(gene_name)%>%
                        summarise(npeaks=n(), dtss=min(dtss), .groups="drop")%>%
                        ungroup()%>%as.data.frame()%>%
                        mutate(conditions=ii, grp="1")
            dd <- d1
            dd
})
head(dd)

## comb <- sort(unique(res.DE$comb))
## comb
## dfNew <- map_dfr(comb, function(ii){
##     res2 <- resDP.join%>%dplyr::filter(comb==ii)
##     DEG <- resDE%>%dplyr::filter(comb==ii)%>%dplyr::pull(gene)
##     res2 <- res2%>%mutate(is_DEG=ifelse(gene_name%in%DEG,1,0))
##     ###
##     dx <- map_dfr(c(0,1),function(i){
##         di <- res2%>%dplyr::filter(is_DEG==i)%>%arrange(p.value)
##         ngene <- nrow(di)
##         di <- di%>%mutate(observed=-log10(p.value), expected=-log10(ppoints(ngene)))
##         di
##         })
##     dx
##     })
## head(dfNew)
## table(dfNew$is_DEG)




## #----fisher test function----
## cal.fisher <- function(df){
##        ###
##      resfisher <- map_dfr(1:nrow(df),function(i){
##                dmat <- matrix(as.numeric(df[i,]), 2, 2)
##                      colnames(dmat) <- c("interest", "not.interest")
##                      rownames(dmat) <- c("in.DAR", "not.DAR")
##                      res <- fisher.test(dmat)
##                      res2 <- data.frame(odds=res$estimate, pval.fisher=res$p.value,
##                                                  lower=res$conf.int[1], upper=res$conf.int[2])
##                      res2
##                   })
##        resfisher
##     }




## #----df.fisher FUNCTION----
## #head(res.DE)
## #head(resDP)
## #head(peakAnno)

## cell <- gsub("[0-9]_.*", "", "0_CD4Naive_caffeine")
## cell


## df.fisher <- function(res, resDP, peakAnno){
##        ###
##        comb <- sort(unique(res$comb))
##        df.ls <- lapply(1:length(comb), function(i){
##            ii <- comb[i]
##            cell <- gsub("_.*", "", ii)
##            contrast <- gsub(".*_", "", ii)
##            res2 <- res%>%dplyr::filter(comb==ii) ## dplyr::filter(qval<0.1, abs(beta)>0.5)
##            resDP2 <- resDP%>%
##                dplyr::filter(comb==ii, p.adjusted<0.1, abs(estimate)>0.5)%>%
##                dplyr::left_join(peakAnno, by="peak")
##            if( nrow(resDP2)>5){
##                sigGene <- res2%>%
##                    dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)%>%
##                    dplyr::pull(gene)
##                ###
##                interest.in.DARs <- sum(sigGene%in%resDP2$SYMBOL)
##                interest.not.DARs <- length(sigGene)-interest.in.DARs
##                ###
##                notSig <- setdiff(res2$gene, sigGene)
##                not.interest.in.DARs <- sum(notSig%in%resDP2$SYMBOL)
##                not.interest.not.DARs <- length(notSig)-not.interest.in.DARs
##                df2 <- data.frame("interest.in.DARs"=interest.in.DARs,
##                                  "interest.not.DARs"=interest.not.DARs,
##                                  "not.interest.in.DARs"=not.interest.in.DARs,
##                                  "not.interest.not.DARs"=not.interest.not.DARs)
##                df2$cell <- cell
##                df2$contrast <- contrast
##                df2$comb <- ii
##                }else{
##                    df2 <- NA
##                }
##                df2
##        })
##        df.ls <- df.ls[!is.na(df.ls)]
##        df <- do.call(rbind, df.ls)
##        res <- cal.fisher(df[,1:4])
##        res2 <- cbind(df, res)
##        ###
##        as.data.frame(res2)
## }



## #----EXECUTE FUNCTION----
## resFisher <- df.fisher(res.DE, resDP, peakAnno2)
## head(resFisher)
## opfn <- "./2.2_compareRNAandATAC.outs/3.1_enrich.DEG_macs2_0.1_cn.csv"
## write.csv(resFisher, opfn, row.names=F)
