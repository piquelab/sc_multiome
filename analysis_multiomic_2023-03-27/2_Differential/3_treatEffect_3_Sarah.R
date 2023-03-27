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
library(reshape2)
library("ggpubr")

rm(list=ls())

##
outdir <- "./3_treatEffect_3_DEG.outs/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=F)



#####################################################################
#### resDE heatmap   ###
#####################################################################
resDE <- read_rds("3_treatEffect_2_2.outs/resDE.mitofilt.rds")
head(resDE)
nrow(resDE)



##all.DE
all.DE  <- resDE %>% dplyr::pull(gene)
all.DE <- as.character(unique(all.DE))
head(all.DE)
length(all.DE)









#####################################################################
#### Sarah's data <- bulk RNA seq ###
#####################################################################
## #################################
## # DATA
## #################################
## resDEG <- read_rds("../2_Differential/1_DEG.outs/2.0_DESeq.results.rds") %>% as.data.frame()
## head(resDEG)
## nrow(resDEG)
## max(resDEG$estimate)
## min(resDEG$estimate)

## resDEG <- resDEG %>% mutate(comb_2=paste(contrast, Significance, sep="_"))
## head(resDEG)
## length(unique(resDEG$gene))


## resDEG <- resDEG %>% mutate(contrast=gsub("Caffeine", "caffeine", contrast)) %>% mutate(contrast=gsub("VitaminA", "vitA", contrast)) %>% mutate(contrast=gsub("Water", "water", contrast)) %>% mutate(contrast=gsub("Zinc", "zinc", contrast))


## #################################
## # DATA Filtering
## #################################
## resDEG[!is.na(resDEG$p.adjusted) & resDEG$p.adjusted<0.1 & resDEG$estimate > 0.5, "Significance"] <- "Up"
## resDEG[!is.na(resDEG$p.adjusted) & resDEG$p.adjusted<0.1 & resDEG$estimate < -0.5, "Significance"] <- "Down"
## resDEG[!is.na(resDEG$p.adjusted) & resDEG$p.adjusted>=0.1, "Significance"] <- "Not Significant"
## resDEG[is.na(resDEG$p.adjusted),"Significance"] <- "NA"
## table(resDEG$Significance)
## ##
## resDEG[!is.na(resDEG$p.adjusted) & resDEG$p.adjusted<0.1 & resDEG$estimate > 0.25, "Significance.25"] <- "Up"
## resDEG[!is.na(resDEG$p.adjusted) & resDEG$p.adjusted<0.1 & resDEG$estimate < -0.25, "Significance.25"] <- "Down"
## resDEG[!is.na(resDEG$p.adjusted) & resDEG$p.adjusted>=0.1, "Significance.25"] <- "Not Significant"
## resDEG[is.na(resDEG$p.adjusted),"Significance.25"] <- "NA"
## table(resDEG$Significance.25)
## ##
## resDEG[!is.na(resDEG$p.adjusted) & resDEG$p.adjusted<0.1 & resDEG$estimate > 0, "Significance.0"] <- "Up"
## resDEG[!is.na(resDEG$p.adjusted) & resDEG$p.adjusted<0.1 & resDEG$estimate < 0, "Significance.0"] <- "Down"
## resDEG[!is.na(resDEG$p.adjusted) & resDEG$p.adjusted>=0.1, "Significance.0"] <- "Not Significant"
## resDEG[is.na(resDEG$p.adjusted),"Significance.0"] <- "NA"
## table(resDEG$Significance.0)


## ##all.DEG
## all.DEG  <- resDEG %>% dplyr::pull(gene)
## length(all.DEG)
## all.DEG <- as.character(unique(all.DEG))
## head(all.DEG)
## length(all.DEG)


## #################################
## # gene name from biomart
## #################################
## #----can get gene names from bioMART----
## library(biomaRt)
## mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
## genes.table <- try(biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description", "gene_biotype", "chromosome_name", "start_position"), mart = mart, useCache = F))
## genes.table.2 <- genes.table%>%dplyr::filter(!chromosome_name %in% c("X", "Y", "MT"), ensembl_gene_id %in% all.DEG) %>% mutate(gene=ensembl_gene_id) %>% dplyr::select(gene, external_gene_name)
## #no_genes <- genes.table.2$external_gene_name
## #check <- no_genes %in% count_genes
## resDEG <- left_join(resDEG, genes.table.2, by = "gene") %>% mutate(ensemble=gene) %>% mutate(gene=external_gene_name)

## resDEG <- resDEG %>% mutate(contrast.gene=paste0(contrast, ".", gene))
## resDEG <- resDEG %>% mutate(z_score=estimate/stderror)

## opfn <- paste0("./3_treatEffect_3_DEG.outs/DEG.rds")
## write_rds(resDEG, opfn)


resDEG <- read_rds("./3_treatEffect_3_DEG.outs/DEG.rds")

head(resDEG)

#################################
# LOOP
#################################
treats <- c(unique(resDEG$contrast))
treats


sig.color=c("both"="#F8766D",
            "DEG"="#00BA38",
            "NA"="darkgrey",
            "sc"="#619CFF")

for(j in treats){
    resDE2 <- resDE %>% dplyr::filter(contrast==j)
    resDEG2 <- resDEG %>% dplyr::filter(contrast==j)
    resDEG2.sum <- resDEG2 %>% arrange(p.adjusted)
    resDEG2.sum <- resDEG2.sum %>% distinct(gene, .keep_all = TRUE)
    ##
    resDE.DEG <- left_join(resDE2, resDEG2.sum, by="contrast.gene")
    resDE.DEG <- resDE.DEG %>% filter(!is.na(estimate.y))
    ##
    resDE.DEG <- resDE.DEG %>% mutate(Significance=ifelse(Significance.y %in% c("Up", "Down") &  Significance.x %in% c("Up", "Down"), "both", ifelse(Significance.y %in% c("Up", "Down") & (Significance.x %in% c("Not Significant", "NA") | is.na(Significance.x)), "DEG", ifelse((Significance.y %in% c("Not Significant", "NA") | is.na(Significance.y)) & Significance.x %in% c("Up", "Down"), "sc", "NA"))))
    print(paste0(j, "_resDE.DEG$Significance"))
    print(table(resDE.DEG$Significance))
    ##
    resDE.DEG <- resDE.DEG %>% mutate(Significance.25=ifelse(Significance.25.y %in% c("Up", "Down") &  Significance.25.x %in% c("Up", "Down"), "both", ifelse(Significance.25.y %in% c("Up", "Down") & (Significance.25.x %in% c("Not Significant", "NA") | is.na(Significance.25.x)), "DEG", ifelse((Significance.25.y %in% c("Not Significant", "NA") | is.na(Significance.25.y)) & Significance.25.x %in% c("Up", "Down"), "sc", "NA"))))
    print(paste0(j, "_resDE.DEG$Significance.25"))
    print(table(resDE.DEG$Significance.25))
    ##
    resDE.DEG <- resDE.DEG %>% mutate(Significance.0=ifelse(Significance.0.y %in% c("Up", "Down") &  Significance.0.x %in% c("Up", "Down"), "both", ifelse(Significance.0.y %in% c("Up", "Down") & (Significance.0.x %in% c("Not Significant", "NA") | is.na(Significance.0.x)), "DEG", ifelse((Significance.0.y %in% c("Not Significant", "NA") | is.na(Significance.0.y)) & Significance.0.x %in% c("Up", "Down"), "sc", "NA"))))
    print(paste0(j, "_resDE.DEG$Significance.0"))
    print(table(resDE.DEG$Significance.0))
    ##
    ##resDE.DEG <- resDE.DEG %>% dplyr::filter(!is.na(Significance))
    ##head(resDE.DEG)
    ##nrow(resDE.DEG)
    ##table(resDE.DEG$Significance)
    ##
    ##----zscorevszscore----
    MCls <- unique(resDE.DEG$MCls)
    print(MCls)
    comb.2 <- MCls
    ##
    for(i in comb.2){
        logFC <- 0.5
        ##print(i)
        res1 <- resDE.DEG %>% dplyr::filter(MCls==i)
        ## png(paste0("./3_treatEffect_3_DEG.outs/logFC_", logFC, "_all_logFCvslogfc_scatter_", j, "_", i, ".png"))
        ## p <- ggplot(res1, aes(x=estimate.y, y=estimate.x, color=Significance)) +
        ##     scale_colour_manual(values = c("#F8766D", "#00BA38", "darkgrey", "#619CFF")) +
        ##     scale_colour_manual(values = sig.color) +
        ##     geom_point(show.legend = TRUE)+
        ##     xlab("DEG_data") + ylab("sc_data")+
        ##     xlim(-4, 10) + ylim(-4, 5)+
        ##     ggtitle(paste0("logFC_", logFC, "_all_", j, "_", i))+
        ##     theme_bw()+
        ##     theme(plot.title = element_text(color = "black", size = 18, face = "bold"), axis.text = element_text(size = 18), axis.title = element_text(size = 18))
        ## print(p)
        ## dev.off()
        png(paste0("./3_treatEffect_3_DEG.outs/logFC_", logFC, "_all_zvsz_scatter_", j, "_", i, ".png"))
        p <- ggplot(res1, aes(x=z_score.y, y=z_score.x, color=Significance)) +
            ##                    scale_colour_manual(values = c("#F8766D", "#00BA38", "darkgrey", "#619CFF")) +
            scale_colour_manual(values = sig.color) +
            geom_point(show.legend = TRUE)+
            xlab("DEG_data") + ylab("sc_data")+
            ##xlim(-4, 10) + ylim(-4, 5)+
            ggtitle(paste0("logFC_", logFC, "_all_", j, "_", i))+
            theme_bw()+
            theme(plot.title = element_text(color = "black", size = 18, face = "bold"), axis.text = element_text(size = 18), axis.title = element_text(size = 18))
        print(p)
        dev.off()
        ##
        logFC <- 0.25
        #print(i)
        res1 <- resDE.DEG %>% dplyr::filter(MCls==i)
        png(paste0("./3_treatEffect_3_DEG.outs/logFC_", logFC, "_all_zvsz_scatter_", j, "_", i, ".png"))
        p <- ggplot(res1, aes(x=z_score.y, y=z_score.x, color=Significance.25)) +
            ##                    scale_colour_manual(values = c("#F8766D", "#00BA38", "darkgrey", "#619CFF")) +
            scale_colour_manual(values = sig.color) +
            geom_point(show.legend = TRUE)+
            xlab("DEG_data") + ylab("sc_data")+
            ##xlim(-4, 10) + ylim(-4, 5)+
            ggtitle(paste0("logFC_", logFC, "_all_", j, "_", i))+
            theme_bw()+
            theme(plot.title = element_text(color = "black", size = 18, face = "bold"), axis.text = element_text(size = 18), axis.title = element_text(size = 18))
        print(p)
        dev.off()
        ##
        logFC <- 0
        #print(i)
        res1 <- resDE.DEG %>% dplyr::filter(MCls==i)
        png(paste0("./3_treatEffect_3_DEG.outs/logFC_", logFC, "_all_zvsz_scatter_", j, "_", i, ".png"))
        p <- ggplot(res1, aes(x=z_score.y, y=z_score.x, color=Significance)) +
            ##                    scale_colour_manual(values = c("#F8766D", "#00BA38", "darkgrey", "#619CFF")) +
            scale_colour_manual(values = sig.color) +
            geom_point(show.legend = TRUE)+
            xlab("DEG_data") + ylab("sc_data")+
            ##xlim(-4, 10) + ylim(-4, 5)+
            ggtitle(paste0("logFC_", logFC, "_all_", j, "_", i))+
            theme_bw()+
            theme(plot.title = element_text(color = "black", size = 18, face = "bold"), axis.text = element_text(size = 18), axis.title = element_text(size = 18))
        print(p)
        dev.off()
    }
}



#################################
# SCORR LOOP 3
#################################
#----scorr----
#----spearman_corr----
MCls <- unique(resDE$MCls)

MCls <- c("0-CD4Naive", "6-Monocyte")

MCls

comb.2 <- MCls

#treats.DE <- unique(resDE$contrast)
#treats.DE

treats <- unique(resDEG$contrast)
treats

#treats <- c("caffeine", "vitA", "zinc")


## #----all----
## logFC <- "NA"
## include <- "all"


## #----Significance.0----
## logFC <- 0
## conditions <- c("both", "sc", "DEG")
## include <- "union"

## #logFC <- 0
## #conditions <- c("both")
## #include <- "both"


## #----Significance.25 ----
## #logFC <- 0.25
## #conditions <- c("both", "sc", "DEG")
## #include <- "union"

## #logFC <- 0.25
## #conditions <- c("both")
## #include <- "both"


## #----Significance----
## #logFC <- 0.5
## #conditions <- c("both", "sc", "DEG")
## #include <- "union"

## #logFC <- 0.5
## #conditions <- c("both")
## #include <- "both"


scorr_caff_z <- c()
scorr_vitA_z <- c()
scorr_zinc_z <- c()
scorr_caff_z.5 <- c()
scorr_vitA_z.5 <- c()
scorr_zinc_z.5 <- c()
scorr_caff_z.25 <- c()
scorr_vitA_z.25 <- c()
scorr_zinc_z.25 <- c()
scorr_caff_z.0 <- c()
scorr_vitA_z.0 <- c()
scorr_zinc_z.0 <- c()


##head(resDE)
##as.numeric(cor.test(resDE$z_score, resDE$estimate, method = "spearman")$p.value)
##as.numeric(cor.test(resDE$z_score, resDE$estimate, method = "spearman")$estimate)


for (j in treats){
    resDE2 <- resDE %>% dplyr::filter(contrast==j)
    resDEG2 <- resDEG %>% dplyr::filter(contrast==j)
    resDEG2.sum <- resDEG2 %>% arrange(p.adjusted)
    resDEG2.sum <- resDEG2.sum %>% distinct(gene, .keep_all = TRUE)
    ##
    resDE.DEG <- left_join(resDE2, resDEG2.sum, by="contrast.gene")
    resDE.DEG <- resDE.DEG %>% filter(!is.na(estimate.y))
    ##
    resDE.DEG <- resDE.DEG %>% mutate(Significance=ifelse(Significance.y %in% c("Up", "Down") &  Significance.x %in% c("Up", "Down"), "both", ifelse(Significance.y %in% c("Up", "Down") & (Significance.x %in% c("Not Significant", "NA") | is.na(Significance.x)), "DEG", ifelse((Significance.y %in% c("Not Significant", "NA") | is.na(Significance.y)) & Significance.x %in% c("Up", "Down"), "sc", "NA"))))
    resDE.DEG <- resDE.DEG %>% mutate(Significance.25=ifelse(Significance.25.y %in% c("Up", "Down") &  Significance.25.x %in% c("Up", "Down"), "both", ifelse(Significance.25.y %in% c("Up", "Down") & (Significance.25.x %in% c("Not Significant", "NA") | is.na(Significance.25.x)), "DEG", ifelse((Significance.25.y %in% c("Not Significant", "NA") | is.na(Significance.25.y)) & Significance.25.x %in% c("Up", "Down"), "sc", "NA"))))
    resDE.DEG <- resDE.DEG %>% mutate(Significance.0=ifelse(Significance.0.y %in% c("Up", "Down") &  Significance.0.x %in% c("Up", "Down"), "both", ifelse(Significance.0.y %in% c("Up", "Down") & (Significance.0.x %in% c("Not Significant", "NA") | is.na(Significance.0.x)), "DEG", ifelse((Significance.0.y %in% c("Not Significant", "NA") | is.na(Significance.0.y)) & Significance.0.x %in% c("Up", "Down"), "sc", "NA"))))
    ##
    ##resDE.DEG <- resDE.DEG %>% dplyr::filter(!is.na(Significance))
    ##head(resDE.DEG)
    ##nrow(resDE.DEG)
    ##table(resDE.DEG$Significance)
    ##
    ####################
    ## scorr_values ##
    ####################
    for (i in comb.2){
        tryCatch({
        logFC <- "NA"
        include <- "all"
        res1 <- resDE.DEG %>% dplyr::filter(contrast.x==j, MCls==i)
        ##scorr <- as.numeric(cor.test(res1$estimate.y, res1$estimate.x, method = "spearman")$p.value)
        scorr_z <- as.numeric(cor.test(res1$z_score.y, res1$z_score.x, method = "spearman")$p.value)
        ##print(scorr)
        print(scorr_z)
        if (is.numeric(scorr_z) & j=="caffeine"){
            ##scorr_caff <- c(scorr_caff, scorr)
            scorr_caff_z <- c(scorr_caff_z, scorr_z)
        } else if (is.numeric(scorr_z) & j=="vitA"){
            ##scorr_vitA <- c(scorr_vitA, scorr)
            scorr_vitA_z <- c(scorr_vitA_z, scorr_z)
        } else if (is.numeric(scorr_z) & j=="zinc"){
            ##scorr_zinc <- c(scorr_zinc, scorr)
            scorr_zinc_z <- c(scorr_zinc_z, scorr_z)
        }
        ##
        ##
        ##
        ##
        logFC <- 0.5
        conditions <- c("both", "sc", "DEG")
        include <- "union"
        res1 <- resDE.DEG %>% dplyr::filter(contrast.x==j, MCls==i) %>% dplyr::filter(Significance %in% conditions)
        ##scorr <- as.numeric(cor.test(res1$estimate.y, res1$estimate.x, method = "spearman")$p.value)
        scorr_z <- as.numeric(cor.test(res1$z_score.y, res1$z_score.x, method = "spearman")$p.value)
        ##print(scorr)
        print(scorr_z)
        if (is.numeric(scorr_z) & j=="caffeine"){
            ##scorr_caff <- c(scorr_caff, scorr)
            scorr_caff_z.5 <- c(scorr_caff_z.5, scorr_z)
        } else if (is.numeric(scorr_z) & j=="vitA"){
            ##scorr_vitA <- c(scorr_vitA, scorr)
            scorr_vitA_z.5 <- c(scorr_vitA_z.5, scorr_z)
        } else if (is.numeric(scorr_z) & j=="zinc"){
            ##scorr_zinc <- c(scorr_zinc, scorr)
            scorr_zinc_z.5 <- c(scorr_zinc_z.5, scorr_z)
        }
        ##
        ##
        ##
        ##
        logFC <- 0.25
        conditions <- c("both", "sc", "DEG")
        include <- "union"
        res1 <- resDE.DEG %>% dplyr::filter(contrast.x==j, MCls==i) %>% dplyr::filter(Significance.25 %in% conditions)
        ##scorr <- as.numeric(cor.test(res1$estimate.y, res1$estimate.x, method = "spearman")$p.value)
        scorr_z <- as.numeric(cor.test(res1$z_score.y, res1$z_score.x, method = "spearman")$p.value)
        ##print(scorr)
        print(scorr_z)
        if (is.numeric(scorr_z) & j=="caffeine"){
            ##scorr_caff <- c(scorr_caff, scorr)
            scorr_caff_z.25 <- c(scorr_caff_z.25, scorr_z)
        } else if (is.numeric(scorr_z) & j=="vitA"){
            ##scorr_vitA <- c(scorr_vitA, scorr)
            scorr_vitA_z.25 <- c(scorr_vitA_z.25, scorr_z)
        } else if (is.numeric(scorr_z) & j=="zinc"){
            ##scorr_zinc <- c(scorr_zinc, scorr)
            scorr_zinc_z.25 <- c(scorr_zinc_z.25, scorr_z)
        }
        ##
        ##
        ##
        ##
        logFC <- 0
        conditions <- c("both", "sc", "DEG")
        include <- "union"
        res1 <- resDE.DEG %>% dplyr::filter(contrast.x==j, MCls==i) %>% dplyr::filter(Significance.0 %in% conditions)
        ##scorr <- as.numeric(cor.test(res1$estimate.y, res1$estimate.x, method = "spearman")$p.value)
        scorr_z <- as.numeric(cor.test(res1$z_score.y, res1$z_score.x, method = "spearman")$p.value)
        ##print(scorr)
        print(scorr_z)
        if (is.numeric(scorr_z) & j=="caffeine"){
            ##scorr_caff <- c(scorr_caff, scorr)
            scorr_caff_z.0 <- c(scorr_caff_z.0, scorr_z)
        } else if (is.numeric(scorr_z) & j=="vitA"){
            ##scorr_vitA <- c(scorr_vitA, scorr)
            scorr_vitA_z.0 <- c(scorr_vitA_z.0, scorr_z)
        } else if (is.numeric(scorr_z) & j=="zinc"){
            ##scorr_zinc <- c(scorr_zinc, scorr)
            scorr_zinc_z.0 <- c(scorr_zinc_z.0, scorr_z)
        } 
        },  error=function(e){cat("ERROR", "\n")})
    }
}

scorr_caff_z
scorr_vitA_z
scorr_zinc_z

scorr_caff_z.5
scorr_vitA_z.5
scorr_zinc_z.5

scorr_caff_z.25
scorr_vitA_z.25
scorr_zinc_z.25

scorr_caff_z.0
scorr_vitA_z.0
scorr_zinc_z.0


#################
## out of loop ##
#################
##############################
# all zscore vs zscore plot ##
##############################
logFC <- "NA"
include <- "all"
scorr_data <- data.frame(MCls, scorr_caff_z,
                         scorr_vitA_z, scorr_zinc_z)
is.num <- sapply(scorr_data, is.numeric)
scorr_data[is.num] <- lapply(scorr_data[is.num], round, 2)
##
sig.melt <- melt(scorr_data)
##
sig.p <- ggplot(sig.melt, aes(x = variable, y = MCls)) +
      geom_tile(aes(fill=value), color="black")+
      scale_fill_gradient(low = "white", high = "white")+
      geom_text(aes(label = value), color = "black", size = 5) +
      ggtitle(paste0("logFC_", logFC, "_", include, "_zscorevszscore_corr_DEG_vs_sc")) +
      theme(axis.text.x = element_text(angle = 90),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    legend.position="none")
##
figure <- ggarrange(sig.p + font("x.text", size = 14))
png(paste0("./3_treatEffect_3_DEG.outs/logFC_", logFC, "_", include, "_zscorevszscore_corr.pvalue.png"), width=500, height=500, pointsize=16, res=125)
print(figure)
dev.off()






logFC <- 0.5
include <- "union"
scorr_data <- data.frame(MCls, scorr_caff_z.5,
                         scorr_vitA_z.5, scorr_zinc_z.5)
is.num <- sapply(scorr_data, is.numeric)
scorr_data[is.num] <- lapply(scorr_data[is.num], round, 2)
##
sig.melt <- melt(scorr_data)
##
sig.p <- ggplot(sig.melt, aes(x = variable, y = MCls)) +
      geom_tile(aes(fill=value), color="black")+
      scale_fill_gradient(low = "white", high = "white")+
      geom_text(aes(label = value), color = "black", size = 5) +
      ggtitle(paste0("logFC_", logFC, "_", include, "_zscorevszscore_corr_DEG_vs_sc")) +
      theme(axis.text.x = element_text(angle = 90),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    legend.position="none")
##
figure <- ggarrange(sig.p + font("x.text", size = 14))
png(paste0("./3_treatEffect_3_DEG.outs/logFC_", logFC, "_", include, "_zscorevszscore_corr.pvalue.png"), width=500, height=500, pointsize=16, res=125)
print(figure)
dev.off()





logFC <- 0.25
include <- "union"
scorr_data <- data.frame(MCls, scorr_caff_z.25,
                         scorr_vitA_z.25, scorr_zinc_z.25)
is.num <- sapply(scorr_data, is.numeric)
scorr_data[is.num] <- lapply(scorr_data[is.num], round, 2)
##
sig.melt <- melt(scorr_data)
##
sig.p <- ggplot(sig.melt, aes(x = variable, y = MCls)) +
      geom_tile(aes(fill=value), color="black")+
      scale_fill_gradient(low = "white", high = "white")+
      geom_text(aes(label = value), color = "black", size = 5) +
      ggtitle(paste0("logFC_", logFC, "_", include, "_zscorevszscore_corr_DEG_vs_sc")) +
      theme(axis.text.x = element_text(angle = 90),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    legend.position="none")
##
figure <- ggarrange(sig.p + font("x.text", size = 14))
png(paste0("./3_treatEffect_3_DEG.outs/logFC_", logFC, "_", include, "_zscorevszscore_corr.pvalue.png"), width=500, height=500, pointsize=16, res=125)
print(figure)
dev.off()







logFC <- 0
include <- "union"
scorr_data <- data.frame(MCls, scorr_caff_z.0,
                         scorr_vitA_z.0, scorr_zinc_z.0)
is.num <- sapply(scorr_data, is.numeric)
scorr_data[is.num] <- lapply(scorr_data[is.num], round, 2)
##
sig.melt <- melt(scorr_data)
##
sig.p <- ggplot(sig.melt, aes(x = variable, y = MCls)) +
      geom_tile(aes(fill=value), color="black")+
      scale_fill_gradient(low = "white", high = "white")+
      geom_text(aes(label = value), color = "black", size = 5) +
      ggtitle(paste0("logFC_", logFC, "_", include, "_zscorevszscore_corr_DEG_vs_sc")) +
      theme(axis.text.x = element_text(angle = 90),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    legend.position="none")
##
figure <- ggarrange(sig.p + font("x.text", size = 14))
png(paste0("./3_treatEffect_3_DEG.outs/logFC_", logFC, "_", include, "_zscorevszscore_corr.pvalue.png"), width=500, height=500, pointsize=16, res=125)
print(figure)
dev.off()








####################################################################
# SCORR LOOP 3 - for all treatments(sc) vs all treatments(Sarah)
####################################################################
#----scorr----
#----spearman_corr----
MCls <- unique(resDE$MCls)

MCls <- c("0-CD4Naive", "6-Monocyte")

MCls

comb.2 <- MCls

treats.DE <- unique(resDE$contrast)
treats.DE

treats <- unique(resDEG$contrast)
treats

#treats <- c("caffeine", "vitA", "zinc")


## scorr_caff_z <- c()
## scorr_vitA_z <- c()
## scorr_zinc_z <- c()
## scorr_caff_z.5 <- c()
## scorr_vitA_z.5 <- c()
## scorr_zinc_z.5 <- c()
## scorr_caff_z.25 <- c()
## scorr_vitA_z.25 <- c()
## scorr_zinc_z.25 <- c()
## scorr_caff_z.0 <- c()
## scorr_vitA_z.0 <- c()
## scorr_zinc_z.0 <- c()


scorr_caff_z <- c()
scorr_nic_z <- c()
scorr_vitA_z <- c()
scorr_vitD_z <- c()
scorr_vitE_z <- c()
scorr_water_z <- c()
scorr_zinc_z <- c()


##head(resDE)
##as.numeric(cor.test(resDE$z_score, resDE$estimate, method = "spearman")$p.value)
##as.numeric(cor.test(resDE$z_score, resDE$estimate, method = "spearman")$estimate)

#head(resDE)
#head(resDEG)

oneMCl = "0-CD4Naive"

oneMCl = "6-Monocyte"

for (j in treats.DE){
    print(j)
    resDE2 <- resDE %>% dplyr::filter(MCls==oneMCl, contrast==j)
    for (k in treats){
        print(k)
        resDEG2 <- resDEG %>% dplyr::filter(contrast==k)
        resDEG2.sum <- resDEG2 %>% arrange(p.adjusted)
        resDEG2.sum <- resDEG2.sum %>% distinct(gene, .keep_all = TRUE)
        ##
        resDE.DEG <- left_join(resDE2, resDEG2.sum, by="gene")
        resDE.DEG <- resDE.DEG %>% filter(!is.na(estimate.y))
        ##
        #resDE.DEG <- resDE.DEG %>% mutate(Significance=ifelse(Significance.y %in% c("Up", "Down") &  Significance.x %in% c("Up", "Down"), "both", ifelse(Significance.y %in% c("Up", "Down") & (Significance.x %in% c("Not Significant", "NA") | is.na(Significance.x)), "DEG", ifelse((Significance.y %in% c("Not Significant", "NA") | is.na(Significance.y)) & Significance.x %in% c("Up", "Down"), "sc", "NA"))))
        #resDE.DEG <- resDE.DEG %>% mutate(Significance.25=ifelse(Significance.25.y %in% c("Up", "Down") &  Significance.25.x %in% c("Up", "Down"), "both", ifelse(Significance.25.y %in% c("Up", "Down") & (Significance.25.x %in% c("Not Significant", "NA") | is.na(Significance.25.x)), "DEG", ifelse((Significance.25.y %in% c("Not Significant", "NA") | is.na(Significance.25.y)) & Significance.25.x %in% c("Up", "Down"), "sc", "NA"))))
        #resDE.DEG <- resDE.DEG %>% mutate(Significance.0=ifelse(Significance.0.y %in% c("Up", "Down") &  Significance.0.x %in% c("Up", "Down"), "both", ifelse(Significance.0.y %in% c("Up", "Down") & (Significance.0.x %in% c("Not Significant", "NA") | is.na(Significance.0.x)), "DEG", ifelse((Significance.0.y %in% c("Not Significant", "NA") | is.na(Significance.0.y)) & Significance.0.x %in% c("Up", "Down"), "sc", "NA"))))
    ##
    ##resDE.DEG <- resDE.DEG %>% dplyr::filter(!is.na(Significance))
    ##head(resDE.DEG)
    ##nrow(resDE.DEG)
    ##table(resDE.DEG$Significance)
    ##
    ####################
    ## scorr_values ##
    ####################
    #for (i in comb.2){
        tryCatch({
        logFC <- "NA"
        include <- "all"
        #res1 <- resDE.DEG %>% dplyr::filter(contrast.x==j, MCls==i)
        res1 <- resDE.DEG
        ##scorr <- as.numeric(cor.test(res1$estimate.y, res1$estimate.x, method = "spearman")$p.value)
        ##print(scorr)
        scorr_z <- as.numeric(cor.test(res1$z_score.y, res1$z_score.x, method = "spearman")$estimate)
         print(scorr_z)
        if (is.numeric(scorr_z) & j=="caffeine"){
            ##scorr_caff <- c(scorr_caff, scorr)
            scorr_caff_z <- c(scorr_caff_z, scorr_z)
        }  else if (is.numeric(scorr_z) & j=="nicotine"){
            ##scorr_nic <- c(scorr_nic, scorr)
            scorr_nic_z <- c(scorr_nic_z, scorr_z)
        } else if (is.numeric(scorr_z) & j=="vitA"){
            ##scorr_vitA <- c(scorr_vitA, scorr)
            scorr_vitA_z <- c(scorr_vitA_z, scorr_z)
        }  else if (is.numeric(scorr_z) & j=="vitD"){
            ##scorr_vitD <- c(scorr_vitD, scorr)
            scorr_vitD_z <- c(scorr_vitD_z, scorr_z)
        } else if (is.numeric(scorr_z) & j=="vitE"){
            ##scorr_vitE <- c(scorr_vitE, scorr)
            scorr_vitE_z <- c(scorr_vitE_z, scorr_z)
        } else if (is.numeric(scorr_z) & j=="water"){
            ##scorr_water <- c(scorr_water, scorr)
            scorr_water_z <- c(scorr_water_z, scorr_z)
        } else if (is.numeric(scorr_z) & j=="zinc"){
            ##scorr_zinc <- c(scorr_zinc, scorr)
            scorr_zinc_z <- c(scorr_zinc_z, scorr_z)
        }
        },  error=function(e){cat("ERROR", "\n")})
    }
}

scorr_caff_z
scorr_vitA_z
scorr_zinc_z
scorr_nic_z
scorr_vitD_z
scorr_vitE_z
scorr_water_z




## scorr_caff_z.5
## scorr_vitA_z.5
## scorr_zinc_z.5
## scorr_caff_z.25
## scorr_vitA_z.25
## scorr_zinc_z.25
## scorr_caff_z.0
## scorr_vitA_z.0
## scorr_zinc_z.0


#################
## out of loop ##
#################
##############################
# all zscore vs zscore plot ##
##############################
logFC <- "NA"
include <- "all"
scorr_data <- data.frame(treats, scorr_caff_z, scorr_nic_z,
                         scorr_vitA_z, scorr_vitD_z,
                         scorr_vitE_z, scorr_water_z, scorr_zinc_z
                         )
scorr_data
is.num <- sapply(scorr_data, is.numeric)
scorr_data[is.num] <- lapply(scorr_data[is.num], round, 2)
##
sig.melt <- melt(scorr_data)
sig.melt

##
sig.p <- ggplot(sig.melt, aes(x = treats, y = variable)) +
      geom_tile(aes(fill=value), color="black")+
      scale_fill_gradient(low = "white", high = "white")+
      geom_text(aes(label = value), color = "black", size = 5) +
      ggtitle(paste0("logFC_", logFC, "_", include, "_zscorevszscore_corr_DEG_vs_sc")) +
      theme(axis.text.x = element_text(angle = 90),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    legend.position="none")
##
figure <- ggarrange(sig.p + font("x.text", size = 14))
png(paste0("./3_treatEffect_3_DEG.outs/logFC_", logFC, "_", include, "_zscorevszscore_corr.monocyte.allcomb.png"), width=500, height=500, pointsize=16, res=125)
print(figure)
dev.off()






logFC <- 0.5
include <- "union"
scorr_data <- data.frame(MCls, scorr_caff_z.5,
                         scorr_vitA_z.5, scorr_zinc_z.5)
is.num <- sapply(scorr_data, is.numeric)
scorr_data[is.num] <- lapply(scorr_data[is.num], round, 2)
##
sig.melt <- melt(scorr_data)
##
sig.p <- ggplot(sig.melt, aes(x = variable, y = MCls)) +
      geom_tile(aes(fill=value), color="black")+
      scale_fill_gradient(low = "white", high = "white")+
      geom_text(aes(label = value), color = "black", size = 5) +
      ggtitle(paste0("logFC_", logFC, "_", include, "_zscorevszscore_corr_DEG_vs_sc")) +
      theme(axis.text.x = element_text(angle = 90),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    legend.position="none")
##
figure <- ggarrange(sig.p + font("x.text", size = 14))
png(paste0("./3_treatEffect_3_DEG.outs/logFC_", logFC, "_", include, "_zscorevszscore_corr.pvalue.png"), width=500, height=500, pointsize=16, res=125)
print(figure)
dev.off()





logFC <- 0.25
include <- "union"
scorr_data <- data.frame(MCls, scorr_caff_z.25,
                         scorr_vitA_z.25, scorr_zinc_z.25)
is.num <- sapply(scorr_data, is.numeric)
scorr_data[is.num] <- lapply(scorr_data[is.num], round, 2)
##
sig.melt <- melt(scorr_data)
##
sig.p <- ggplot(sig.melt, aes(x = variable, y = MCls)) +
      geom_tile(aes(fill=value), color="black")+
      scale_fill_gradient(low = "white", high = "white")+
      geom_text(aes(label = value), color = "black", size = 5) +
      ggtitle(paste0("logFC_", logFC, "_", include, "_zscorevszscore_corr_DEG_vs_sc")) +
      theme(axis.text.x = element_text(angle = 90),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    legend.position="none")
##
figure <- ggarrange(sig.p + font("x.text", size = 14))
png(paste0("./3_treatEffect_3_DEG.outs/logFC_", logFC, "_", include, "_zscorevszscore_corr.pvalue.png"), width=500, height=500, pointsize=16, res=125)
print(figure)
dev.off()







logFC <- 0
include <- "union"
scorr_data <- data.frame(MCls, scorr_caff_z.0,
                         scorr_vitA_z.0, scorr_zinc_z.0)
is.num <- sapply(scorr_data, is.numeric)
scorr_data[is.num] <- lapply(scorr_data[is.num], round, 2)
##
sig.melt <- melt(scorr_data)
##
sig.p <- ggplot(sig.melt, aes(x = variable, y = MCls)) +
      geom_tile(aes(fill=value), color="black")+
      scale_fill_gradient(low = "white", high = "white")+
      geom_text(aes(label = value), color = "black", size = 5) +
      ggtitle(paste0("logFC_", logFC, "_", include, "_zscorevszscore_corr_DEG_vs_sc")) +
      theme(axis.text.x = element_text(angle = 90),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    legend.position="none")
##
figure <- ggarrange(sig.p + font("x.text", size = 14))
png(paste0("./3_treatEffect_3_DEG.outs/logFC_", logFC, "_", include, "_zscorevszscore_corr.pvalue.png"), width=500, height=500, pointsize=16, res=125)
print(figure)
dev.off()

