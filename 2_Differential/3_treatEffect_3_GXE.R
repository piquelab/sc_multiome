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
outdir <- "./3_treatEffect_3_GXE.outs/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=F)



#####################################################################
#### resDE heatmap   ###
#####################################################################
## #----resDE----
## #resDE <- read_rds("../2_Differential/1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn.rds") %>% as.data.frame()
## #resDE <- read_rds("../2_Differential/1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn_sampleID+treat_allvsetOH.rds") %>% as.data.frame()
## resDE <- read_rds("../2_Differential/1_DiffRNA_2.outs/2.0_DESeq.results_clusters_treat&ind_filtered_0.1_cn_sampleID+new_treat.rds") %>% as.data.frame()
## head(resDE)

## #resDE <- resDE %>% mutate(MCls=gsub("_", "-", MCls)) %>% mutate(comb=paste(MCls, contrast, sep="_"))
## #%>% drop_na(p.value)

## resDE <- resDE %>% mutate(MCls=gsub("_", "-", MCls)) %>% mutate(comb=paste(contrast, MCls, sep="_"))
## #%>% drop_na(p.value)
## head(resDE)
## nrow(resDE)
## max(resDE$estimate)
## min(resDE$estimate)


## #resDE$Significance <- NULL
## #head(resDE)

## resDE[!is.na(resDE$p.adjusted) & resDE$p.adjusted<0.1 & resDE$estimate > 0.5, "Significance"] <- "Up"
## resDE[!is.na(resDE$p.adjusted) & resDE$p.adjusted<0.1 & resDE$estimate < -0.5, "Significance"] <- "Down"
## resDE[!is.na(resDE$p.adjusted) & resDE$p.adjusted>=0.1, "Significance"] <- "Not Significant"
## resDE[is.na(resDE$p.adjusted),"Significance"] <- "NA"
## head(resDE)
## nrow(resDE)
## table(resDE$Significance)


## resDE[!is.na(resDE$p.adjusted) & resDE$p.adjusted<0.1 & resDE$estimate > 0.25, "Significance.25"] <- "Up"
## resDE[!is.na(resDE$p.adjusted) & resDE$p.adjusted<0.1 & resDE$estimate < -0.25, "Significance.25"] <- "Down"
## resDE[!is.na(resDE$p.adjusted) & resDE$p.adjusted>=0.1, "Significance.25"] <- "Not Significant"
## resDE[is.na(resDE$p.adjusted),"Significance.25"] <- "NA"
## #head(resDE)
## table(resDE$Significance.25)


## resDE[!is.na(resDE$p.adjusted) & resDE$p.adjusted<0.1 & resDE$estimate > 0, "Significance.0"] <- "Up"
## resDE[!is.na(resDE$p.adjusted) & resDE$p.adjusted<0.1 & resDE$estimate < 0, "Significance.0"] <- "Down"
## resDE[!is.na(resDE$p.adjusted) & resDE$p.adjusted>=0.1, "Significance.0"] <- "Not Significant"
## resDE[is.na(resDE$p.adjusted),"Significance.0"] <- "NA"
## #head(resDE)
## table(resDE$Significance.0)


## resDE <- resDE %>% mutate(comb_2=paste(comb, Significance, sep="_"))
## head(resDE)
## table(resDE$Significance)
## length(unique(resDE$gene))


## resDE <- resDE %>% mutate(contrast.gene=paste0(contrast, ".", gene))
## head(resDE)
## nrow(resDE)


## #dplyr::filter does not work for vector
## #mit.genes <- resDE$gene %>% dplyr::filter(substr(resDE$gene, 1, 2)=="MT")
## #the foll works but dont need as doing it other way
## #mit.genes <- resDE$gene[substr(resDE$gene, 1, 3)=="MT-"]
## #mit.genes

## #this select all the rows with "MT" present anywhere in string
## #library(data.table)
## #mit.resDE <- resDE[resDE$gene %like% "MT", ]
## #mit.resDE$gene


## resDE <- resDE %>% mutate(mito=ifelse(substr(resDE$gene, 1, 3)=="MT-", 1, 0))
## #head(resDE)
## nrow(resDE)
## table(resDE$mito)

## resDE <- resDE %>% mutate(z_score=estimate/stderror)
## head(resDE)
## nrow(resDE)


## #----mitofilt resDE----
## resDE <- resDE %>% dplyr::filter(resDE$mito==0)
## #head(resDE)
## nrow(resDE)
## length(unique(resDE$gene))

## table(resDE$mito)
## table(resDE$Significance)
## table(resDE$Significance.25)
## table(resDE$Significance.0)






## #opfn <- paste0("3_treatEffect_2_2.outs/resDE_contrastetOH.mitofilt.rds")
## #write_rds(resDE, opfn)


#resDE <- read_rds("3_treatEffect_2_2.outs/resDE.mitofilt.rds")
#resDE <- read_rds("3_treatEffect_2_2.outs/resDE_contrastetOH.mitofilt.rds")

resDE <- read_rds("3_treatEffect_2_2.outs/resDE_control.mitofilt.rds")
head(resDE)
nrow(resDE)

#----all.DE----
all.DE  <- resDE %>% dplyr::pull(gene)
all.DE <- as.character(unique(all.DE))
head(all.DE)
length(all.DE)









#####################################################################
#### GxE browser data <- bulk RNA seq ###
#####################################################################
#################################
# DATA
#################################
## resGXE.caff <- read.table("../2_Differential/DEG_results/treatments/DP6_DEG_stats_T13C1.txt", header=T)
## resGXE.caff <- resGXE.caff %>% mutate(contrast="caffeine")
## resGXE.nic <- read.table("../2_Differential/DEG_results/treatments/DP6_DEG_stats_T14C1.txt", header=T)
## resGXE.nic <- resGXE.nic %>% mutate(contrast="nicotine")
## ##
## resGXE.vitA <- read.table("../2_Differential/DEG_results/treatments/DP6_DEG_stats_T6C1.txt", header=T)
## resGXE.vitA <- resGXE.vitA %>% mutate(contrast="vitA")
## ##
## resGXE.vitD <- read.table("../2_Differential/DEG_results/treatments/DP6_DEG_stats_T23C1.txt", header=T)
## resGXE.vitD <- resGXE.vitD %>% mutate(contrast="vitD")
## ##
## resGXE.vitE <- read.table("../2_Differential/DEG_results/treatments/DP6_DEG_stats_T7C1.txt", header=T)
## resGXE.vitE <- resGXE.vitE %>% mutate(contrast="vitE")
## ##
## resGXE.zinc <- read.table("../2_Differential/DEG_results/treatments/DP6_DEG_stats_T20C1.txt", header=T)
## resGXE.zinc <- resGXE.zinc %>% mutate(contrast="zinc")
## resGXE.water <- read.table("../2_Differential/DEG_results/DP6_DEG_stats_CO2.txt", header=T)
## resGXE.water <- resGXE.water %>% mutate(contrast="water")


## head(resGXE.caff)

## nrow(resGXE.caff)
## nrow(resGXE.nic)
## nrow(resGXE.vitA)
## nrow(resGXE.vitD)
## nrow(resGXE.vitE)
## nrow(resGXE.zinc)



## #################################
## # DATA MERGE
## #################################
## resGXE <- rbind(resGXE.caff, resGXE.nic, resGXE.vitA, resGXE.vitD, resGXE.vitE, resGXE.water, resGXE.zinc)
## nrow(resGXE)
## resGXE <- resGXE %>% mutate(contrast.gene=paste0(contrast, ".", g.id))

## head(resGXE)


## #################################
## # DATA Filtering
## #################################
## resGXE[!is.na(resGXE$padj) & resGXE$padj<0.1 & resGXE$logFC > 0.5, "Significance"] <- "Up"
## resGXE[!is.na(resGXE$padj) & resGXE$padj<0.1 & resGXE$logFC < -0.5, "Significance"] <- "Down"
## resGXE[!is.na(resGXE$padj) & resGXE$padj>=0.1, "Significance"] <- "Not Significant"
## resGXE[is.na(resGXE$padj),"Significance"] <- "NA"
## head(resGXE)
## nrow(resGXE)
## table(resGXE$Significance)
## ##
## resGXE[!is.na(resGXE$padj) & resGXE$padj<0.1 & resGXE$logFC > 0.25, "Significance.25"] <- "Up"
## resGXE[!is.na(resGXE$padj) & resGXE$padj<0.1 & resGXE$logFC < -0.25, "Significance.25"] <- "Down"
## resGXE[!is.na(resGXE$padj) & resGXE$padj>=0.1, "Significance.25"] <- "Not Significant"
## resGXE[is.na(resGXE$padj),"Significance.25"] <- "NA"
## head(resGXE)
## nrow(resGXE)
## table(resGXE$Significance.25)
## ##
## resGXE[!is.na(resGXE$padj) & resGXE$padj<0.1 & resGXE$logFC > 0, "Significance.0"] <- "Up"
## resGXE[!is.na(resGXE$padj) & resGXE$padj<0.1 & resGXE$logFC < 0, "Significance.0"] <- "Down"
## resGXE[!is.na(resGXE$padj) & resGXE$padj>=0.1, "Significance.0"] <- "Not Significant"
## resGXE[is.na(resGXE$padj),"Significance.0"] <- "NA"
## head(resGXE)
## nrow(resGXE)
## table(resGXE$Significance.0)


## #length(unique(resGXE$g.id))
## #length(unique(resGXE$ensg))
## #table(resGXE$ensg)
## #resGXE %>% dplyr::filter(ensg=="ENSG00000273492")
## #write.csv(resGXE,"./1_GXE.outs/resGXE.csv", row.names = FALSE)


## write_rds(resGXE, "./3_treatEffect_3_GXE.outs/resGXE.rds")

resGXE <- read_rds("./3_treatEffect_3_GXE.outs/resGXE.rds")
head(resGXE)








#################################
# LOOP
#################################
head(resDE)
head(resGXE)

#treats <- unique(resDE$contrast)
#treats

treats <- c(unique(resGXE$contrast))
treats

head(resDE)

#resDE2 <- resDE %>% dplyr::filter(contrast=="caffeine") 
#head(resDE2)

sig.color=c("both"="#F8766D",
            "GXE"="#00BA38",
            "NA"="darkgrey",
            "sc"="#619CFF")

for(j in treats){
    ##print(j)
    resDE2 <- resDE %>% dplyr::filter(contrast==j)
    ##print(head(resDE2))
    resGXE2 <- resGXE %>% dplyr::filter(contrast==j)
    ##
    ##resGXE2.sum <- resGXE2 %>% arrange(padj)
    ## print(paste0("sc_all_mitofilt_", i))
    ## table(resGXE2.sum$Significance)
    ## table(resGXE2.sum$Significance.25)
    ## table(resGXE2.sum$Significance.0)
    ##
    ##resGXE2.sum <- resGXE2.sum %>% distinct(g.id, .keep_all = TRUE)
    ## print(paste0("GXE_all_", i))
    ## table(resGXE2.sum$Significance)
    ## table(resGXE2.sum$Significance.25)
    ## table(resGXE2.sum$Significance.0)
    ## print(paste0("GXE_unique_gids_", i))
    ## table(resGXE2.sum$Significance)
    ## table(resGXE2.sum$Significance.25)
    ## table(resGXE2.sum$Significance.0)
    ##
    resGXE2.sum <- resGXE2 %>% mutate(sig=ifelse(Significance.25 %in% c("Up", "Down"), "sig", "nonsig"))
    resGXE2.sum <- resGXE2.sum %>% arrange(desc(logFC))
    resGXE2.sum <- resGXE2.sum %>% arrange(padj)
    resGXE2.sum <- resGXE2.sum %>% arrange(desc(sig))
    resGXE2.sum <- resGXE2.sum %>% distinct(ensg, .keep_all = TRUE)
    ##
    resDE.GXE <- left_join(resDE2, resGXE2.sum, by="contrast.gene")
    resDE.GXE <- resDE.GXE %>% filter(!is.na(logFC))
    ##print(head(resDE.GXE))
    ##
    resDE.GXE <- resDE.GXE %>% mutate(Significance=ifelse(Significance.y %in% c("Up", "Down") &  Significance.x %in% c("Up", "Down"), "both", ifelse(Significance.y %in% c("Up", "Down") & (Significance.x %in% c("Not Significant", "NA") | is.na(Significance.x)), "GXE", ifelse((Significance.y %in% c("Not Significant", "NA") | is.na(Significance.y)) & Significance.x %in% c("Up", "Down"), "sc", "NA"))))
    print(paste0(j, "_resDE.GXE$Significance"))
    print(table(resDE.GXE$Significance))
    ##
    resDE.GXE <- resDE.GXE %>% mutate(Significance.25=ifelse(Significance.25.y %in% c("Up", "Down") &  Significance.25.x %in% c("Up", "Down"), "both", ifelse(Significance.25.y %in% c("Up", "Down") & (Significance.25.x %in% c("Not Significant", "NA") | is.na(Significance.25.x)), "GXE", ifelse((Significance.25.y %in% c("Not Significant", "NA") | is.na(Significance.25.y)) & Significance.25.x %in% c("Up", "Down"), "sc", "NA"))))
    print(paste0(j, "_resDE.GXE$Significance.25"))
    print(table(resDE.GXE$Significance.25))
    ##
    resDE.GXE <- resDE.GXE %>% mutate(Significance.0=ifelse(Significance.0.y %in% c("Up", "Down") &  Significance.0.x %in% c("Up", "Down"), "both", ifelse(Significance.0.y %in% c("Up", "Down") & (Significance.0.x %in% c("Not Significant", "NA") | is.na(Significance.0.x)), "GXE", ifelse((Significance.0.y %in% c("Not Significant", "NA") | is.na(Significance.0.y)) & Significance.0.x %in% c("Up", "Down"), "sc", "NA"))))
    print(paste0(j, "_resDE.GXE$Significance.0"))
    print(table(resDE.GXE$Significance.0))
    ##
    ##resDE.GXE <- resDE.GXE %>% dplyr::filter(!is.na(Significance))
    ##head(resDE.GXE)
    ##nrow(resDE.GXE)
    ##table(resDE.GXE$Significance)
    ##
    ##----logFC 0 and logFCvlogFC----
    logFC <- 0.5
    MCls <- unique(resDE.GXE$MCls)
    print(MCls)
    comb.2 <- MCls
    ##
    for(i in comb.2){
        print(i)
        res1 <- resDE.GXE %>% dplyr::filter(MCls==i)
        png(paste0("./3_treatEffect_3_GXE.outs/logFC_", logFC, "_all_scatter_", j, "_", i, "_control.png"))
        p <- ggplot(res1, aes(x=logFC, y=estimate, color=Significance)) +
            ##scale_colour_manual(values = c("#F8766D", "#00BA38", "darkgrey", "#619CFF")) +
            scale_colour_manual(values = sig.color) +
            geom_point(show.legend = TRUE)+
            xlab("GXE_data") + ylab("sc_data")+
            ##xlim(-4, 10) + ylim(-4, 5)+
            ggtitle(paste0("logFC_", logFC, "_all_", j, "_", i))+
            theme_bw()+
            theme(plot.title = element_text(color = "black", size = 18, face = "bold"), axis.text = element_text(size = 18), axis.title = element_text(size = 18))
        print(p)
        dev.off()
        ##
        png(paste0("./3_treatEffect_3_GXE.outs/logFC_", logFC, "_all_zvslogfc_scatter_", j, "_", i, "_control.png"))
        p <- ggplot(res1, aes(x=logFC, y=z_score, color=Significance)) +
            ##scale_colour_manual(values = c("#F8766D", "#00BA38", "darkgrey", "#619CFF")) +
            scale_colour_manual(values = sig.color) +
            geom_point(show.legend = TRUE)+
            xlab("GXE_data") + ylab("sc_data")+
            ##xlim(-4, 10) + ylim(-4, 5)+
            ggtitle(paste0("logFC_", logFC, "_all_zvslogFC_", j, "_", i))+
            theme_bw()+
            theme(plot.title = element_text(color = "black", size = 18, face = "bold"), axis.text = element_text(size = 18), axis.title = element_text(size = 18))
        print(p)
        dev.off()
    }
}

end
    
##         #----logFC 0 and zvslogFC----
##         logFC <- 0
##         MCls <- unique(resDE.GXE$MCls)
##         print(MCls)
##         comb.2 <- MCls
##         ## 
##         for(i in comb.2){
##                 print(i)
##                 res1 <- resDE.GXE %>% dplyr::filter(MCls==i)
##                 png(paste0("./3_treatEffect_3_GXE.outs/logFC_", logFC, "_all_zvslogfc_scatter_", j, "_", i, ".png"))
##                 p <- ggplot(res1, aes(x=logFC, y=z_score, color=Significance.0)) +
##                     scale_colour_manual(values = c("#F8766D", "#00BA38", "darkgrey", "#619CFF")) +
##                     geom_point(show.legend = TRUE)+
##                     xlab("GXE_data") + ylab("sc_data")+
## #                    xlim(-4, 10) + ylim(-4, 5)+
##                     ggtitle(paste0("logFC_", logFC, "_all_zvslogFC_", j, "_", i))+
##                     theme_bw()+
##                     theme(plot.title = element_text(color = "black", size = 18, face = "bold"), axis.text = element_text(size = 18), axis.title = element_text(size = 18))
##                 print(p)
##                 dev.off()
##         }}





## #################################
## # SCORR LOOP
## #################################
## #----scorr----
## #----spearman_corr----
## MCls <- unique(resDE$MCls)
## MCls
## comb.2 <- MCls

## treats <- unique(resGXE$contrast)
## treats

## logFC <- 0
## include <- "all"

## logFC <- 0
## conditions <- c("both", "sc", "GXE")
## include <- "union"

## logFC <- 0
## conditions <- c("both")
## include <- "both"



## #colnames(resDE)
## #as.numeric(cor.test(res1$z_score.y, res1$z_score.x, method = "spearman")$estimate)

## #----
## #change 4 things
## # significance
## # y-axis
## # ggtitle
## # png file name
## scorr_caff <- c()
## scorr_nic <- c()
## scorr_vitA <- c()
## scorr_vitD <- c()
## scorr_vitE <- c()
## scorr_zinc <- c()
## for (j in treats){
##         resDE2 <- resDE %>% dplyr::filter(contrast==j)
##         ##
##         resGXE2 <- resGXE %>% dplyr::filter(contrast==j)
##         ##
##         resGXE2.sum <- resGXE2 %>% arrange(padj)
##         ##
##         resGXE2.sum <- resGXE2.sum %>% distinct(g.id, .keep_all = TRUE)
##         ##
##         resDE.GXE <- left_join(resDE2, resGXE2.sum, by="contrast.gene")
##         resDE.GXE <- resDE.GXE %>% filter(!is.na(logFC))
##         ##
##         resDE.GXE <- resDE.GXE %>% mutate(Significance=ifelse(Significance.y %in% c("Up", "Down") &  Significance.x %in% c("Up", "Down"), "both", ifelse(Significance.y %in% c("Up", "Down") & (Significance.x %in% c("Not Significant", "NA") | is.na(Significance.x)), "GXE", ifelse((Significance.y %in% c("Not Significant", "NA") | is.na(Significance.y)) & Significance.x %in% c("Up", "Down"), "sc", "NA"))))
##         #print(paste0(j, "_resDE.GXE$Significance"))
##         #print(table(resDE.GXE$Significance))
##         ##
##         resDE.GXE <- resDE.GXE %>% mutate(Significance.25=ifelse(Significance.25.y %in% c("Up", "Down") &  Significance.25.x %in% c("Up", "Down"), "both", ifelse(Significance.25.y %in% c("Up", "Down") & (Significance.25.x %in% c("Not Significant", "NA") | is.na(Significance.25.x)), "GXE", ifelse((Significance.25.y %in% c("Not Significant", "NA") | is.na(Significance.25.y)) & Significance.25.x %in% c("Up", "Down"), "sc", "NA"))))
##         #print(paste0(j, "_resDE.GXE$Significance.25"))
##         #print(table(resDE.GXE$Significance.25))
##         ##
##         resDE.GXE <- resDE.GXE %>% mutate(Significance.0=ifelse(Significance.0.y %in% c("Up", "Down") &  Significance.0.x %in% c("Up", "Down"), "both", ifelse(Significance.0.y %in% c("Up", "Down") & (Significance.0.x %in% c("Not Significant", "NA") | is.na(Significance.0.x)), "GXE", ifelse((Significance.0.y %in% c("Not Significant", "NA") | is.na(Significance.0.y)) & Significance.0.x %in% c("Up", "Down"), "sc", "NA"))))
##         #print(paste0(j, "_resDE.GXE$Significance.0"))
##         #print(table(resDE.GXE$Significance.0))
##         ##
##         ##
##         ##CHANGE SIGNIFICANCE HERE
##         ##
##         ##
##         for (i in comb.2){
##         res1 <- resDE.GXE %>% dplyr::filter(contrast.x==j, MCls==i)               #                  %>% dplyr::filter(Significance %in% conditions)
##         ##
##         ##
##         ##CHANGE X-AXIS
##         ##
##         ##
##         #scorr <- as.numeric(cor.test(res1$logFC, res1$estimate, method = "spearman")$estimate)
##         scorr <- as.numeric(cor.test(res1$logFC, res1$z_score, method = "spearman")$estimate)
##         ##
##         ##
##         ##
##         ##
##         ##
##         print(scorr)
##         if (is.numeric(scorr) & j=="caffeine"){
##             scorr_caff <- c(scorr_caff, scorr)
##         } else if (is.numeric(scorr) & j=="nicotine"){
##             scorr_nic <- c(scorr_nic, scorr)
##         } else if (is.numeric(scorr) & j=="vitA"){
##             scorr_vitA <- c(scorr_vitA, scorr)
##         } else if (is.numeric(scorr) & j=="vitD"){
##             scorr_vitD <- c(scorr_vitD, scorr)
##         } else if (is.numeric(scorr) & j=="vitE"){
##             scorr_vitE <- c(scorr_vitE, scorr)
##         } else if (is.numeric(scorr) & j=="zinc"){
##             scorr_zinc <- c(scorr_zinc, scorr)
##         }
##     }
## }
## ##
## ##
## ##
## scorr_data <- data.frame(MCls, scorr_caff, scorr_nic,
##                          scorr_vitA, scorr_vitD, scorr_vitE,
##                          scorr_zinc)
## is.num <- sapply(scorr_data, is.numeric)
## scorr_data[is.num] <- lapply(scorr_data[is.num], round, 2)
## ##
## ##sig.matrix <- as.matrix(scorr_data)
## ##sig.matrix
## ##
## sig.melt <- melt(scorr_data)
## #sig.melt
## ##
## ##
## ##
## ##
## ##
## sig.p <- ggplot(sig.melt, aes(x = variable, y = MCls)) +
##     geom_tile(aes(fill=value), color="black")+
##     scale_fill_gradient(low = "white", high = "white")+
##     geom_text(aes(label = value), color = "black", size = 5) +
##     ##
##     ##
##     ##CHANGE GG TITLE
##     ##
##     ##
##     #ggtitle(paste0("logFC_", logFC, "_", include, "_logFC_corr_GXE_vs_sc")) +
##     ggtitle(paste0("logFC_", logFC, "_", include, "_logFCvszscore_corr_GXE_vs_sc")) +
##     ##
##     ##
##     theme(axis.text.x = element_text(angle = 90),
##                     axis.title.x = element_blank(),
##                     axis.title.y = element_blank(),
##                     legend.position="none")
## ##
## figure <- ggarrange(sig.p + font("x.text", size = 14))
## # ncol = 2, nrow = 2)
## ##
## ##
## ##CHANGE FILE NAME
## ##
## ##
## #png(paste0("./3_treatEffect_3_GXE.outs/logFC_", logFC, "_", include, "_logFC_corr.png"), width=500, height=500, pointsize=16, res=125)
## png(paste0("./3_treatEffect_3_GXE.outs/logFC_", logFC, "_", include, "_logFCvszscore_corr.png"), width=500, height=500, pointsize=16, res=125)
## ##
## ##
## print(figure)
## dev.off()







## ## #----heatmap----
## ## mat <- scorr_data
## ## rownames(mat) <- mat$MCls
## ## mat <- mat %>% dplyr::select(!MCls)
## ## mat

## ## mat2 <- as.matrix(mat)
## ## mat2 <- mat2[nrow(mat2):1, ]
## ## mat2


## ## fig1 <- pheatmap(mat2, border_color="NA",
## ##                                   cluster_rows=F, cluster_cols=F,
## ##                                   #                    annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend=T,
## ##                                   show_colnames=T,
## ##                                   show_rownames=T,
## ##                                   na_col="white",
## ##                                   fontsize_col=20,
## ##                                   fontsize_row=10)


## ## figfn <- "./1_GXE.outs/logFC_corr_heatmap.png"
## ## png(figfn, width=1800, height=2100,res=225)
## ## print(fig1)
## ## dev.off() 






## #################################
## # SCORR LOOP 2
## #################################
## #----scorr----
## #----spearman_corr----
## MCls <- unique(resDE$MCls)
## MCls
## comb.2 <- MCls

## treats <- unique(resGXE$contrast)
## treats

## #logFC <- 0
## #conditions <- c("both", "sc", "GXE")
## #include <- "union"


## logFC <- 0
## conditions <- c("both")
## include <- "both"



## logFC <- 0
## include <- "all"


## #colnames(resDE)
## #as.numeric(cor.test(res1$z_score.y, res1$z_score.x, method = "spearman")$estimate)

## #----
## #change 4 things
## # significance
## # y-axis
## # ggtitle
## # png file name
## scorr_caff <- c()
## scorr_nic <- c()
## scorr_vitA <- c()
## scorr_vitD <- c()
## scorr_vitE <- c()
## scorr_zinc <- c()
## for (j in treats){
##             resDE2 <- resDE %>% dplyr::filter(contrast==j)
##                     ##
##                     resGXE2 <- resGXE %>% dplyr::filter(contrast==j)
##                     ##
##                     # resGXE2.sum <- resGXE2 %>% arrange(padj)
##                     # resGXE2.sum <- resGXE2.sum %>% distinct(g.id, .keep_all = TRUE)
##                     ##
##                     resGXE2.sum <- resGXE2 %>% mutate(sig=ifelse(Significance.25 %in% c("Up", "Down"), "sig", "nonsig"))
##                     resGXE2.sum <- resGXE2.sum %>% arrange(desc(logFC))
##                     resGXE2.sum <- resGXE2.sum %>% arrange(padj)
##                     resGXE2.sum <- resGXE2.sum %>% arrange(desc(sig))
##                     resGXE2.sum <- resGXE2.sum %>% distinct(ensg, .keep_all = TRUE)
##                     ##
##                     resDE.GXE <- left_join(resDE2, resGXE2.sum, by="contrast.gene")
##                     resDE.GXE <- resDE.GXE %>% filter(!is.na(logFC))
##                     ##
##                     resDE.GXE <- resDE.GXE %>% mutate(Significance=ifelse(Significance.y %in% c("Up", "Down") &  Significance.x %in% c("Up", "Down"), "both", ifelse(Significance.y %in% c("Up", "Down") & (Significance.x %in% c("Not Significant", "NA") | is.na(Significance.x)), "GXE", ifelse((Significance.y %in% c("Not Significant", "NA") | is.na(Significance.y)) & Significance.x %in% c("Up", "Down"), "sc", "NA"))))
##                     #print(paste0(j, "_resDE.GXE$Significance"))
##                     #print(table(resDE.GXE$Significance))
##                     ##
##                     resDE.GXE <- resDE.GXE %>% mutate(Significance.25=ifelse(Significance.25.y %in% c("Up", "Down") &  Significance.25.x %in% c("Up", "Down"), "both", ifelse(Significance.25.y %in% c("Up", "Down") & (Significance.25.x %in% c("Not Significant", "NA") | is.na(Significance.25.x)), "GXE", ifelse((Significance.25.y %in% c("Not Significant", "NA") | is.na(Significance.25.y)) & Significance.25.x %in% c("Up", "Down"), "sc", "NA"))))
##                     #print(paste0(j, "_resDE.GXE$Significance.25"))
##                     #print(table(resDE.GXE$Significance.25))
##                     ##
##                     resDE.GXE <- resDE.GXE %>% mutate(Significance.0=ifelse(Significance.0.y %in% c("Up", "Down") &  Significance.0.x %in% c("Up", "Down"), "both", ifelse(Significance.0.y %in% c("Up", "Down") & (Significance.0.x %in% c("Not Significant", "NA") | is.na(Significance.0.x)), "GXE", ifelse((Significance.0.y %in% c("Not Significant", "NA") | is.na(Significance.0.y)) & Significance.0.x %in% c("Up", "Down"), "sc", "NA"))))
##                     #print(paste0(j, "_resDE.GXE$Significance.0"))
##                     #print(table(resDE.GXE$Significance.0))
##                     ##
##                     ##
##                     ##CHANGE SIGNIFICANCE HERE
##                     ##
##                     ##
##                     for (i in comb.2){
##                                 res1 <- resDE.GXE %>% dplyr::filter(contrast.x==j, MCls==i)               #                  %>% dplyr::filter(Significance %in% conditions)
##                                         ##
##                                         ##
##                                         ##CHANGE X-AXIS
##                                         ##
##                                         ##
##                                         #scorr <- as.numeric(cor.test(res1$logFC, res1$estimate, method = "spearman")$estimate)
##                                         scorr <- as.numeric(cor.test(res1$logFC, res1$z_score, method = "spearman")$estimate)
##                                         ##
##                                         ##
##                                         ##
##                                         ##
##                                         ##
##                                         print(scorr)
##                                         if (is.numeric(scorr) & j=="caffeine"){
##                                                         scorr_caff <- c(scorr_caff, scorr)
##                                                                 } else if (is.numeric(scorr) & j=="nicotine"){
##                                                                                 scorr_nic <- c(scorr_nic, scorr)
##                                                                                         } else if (is.numeric(scorr) & j=="vitA"){
##                                                                                                         scorr_vitA <- c(scorr_vitA, scorr)
##                                                                                                                 } else if (is.numeric(scorr) & j=="vitD"){
##                                                                                                                                 scorr_vitD <- c(scorr_vitD, scorr)
##                                                                                                                                         } else if (is.numeric(scorr) & j=="vitE"){
##                                                                                                                                                         scorr_vitE <- c(scorr_vitE, scorr)
##                                                                                                                                                                 } else if (is.numeric(scorr) & j=="zinc"){
##                                                                                                                                                                                 scorr_zinc <- c(scorr_zinc, scorr)
##                                                                                                                                                                                         }
##                                     }
##             }
## ##
## ##
## ##
## scorr_data <- data.frame(MCls, scorr_caff, scorr_nic,
##                                                   scorr_vitA, scorr_vitD, scorr_vitE,
##                                                   scorr_zinc)
## is.num <- sapply(scorr_data, is.numeric)
## scorr_data[is.num] <- lapply(scorr_data[is.num], round, 2)
## ##
## ##sig.matrix <- as.matrix(scorr_data)
## ##sig.matrix
## ##
## sig.melt <- melt(scorr_data)
## #sig.melt
## ##
## ##
## ##
## ##
## ##
## sig.p <- ggplot(sig.melt, aes(x = variable, y = MCls)) +
##         geom_tile(aes(fill=value), color="black")+
##         scale_fill_gradient(low = "white", high = "white")+
##         geom_text(aes(label = value), color = "black", size = 5) +
##         ##
##         ##
##         ##CHANGE GG TITLE
##         ##
##         ##
##         #ggtitle(paste0("logFC_", logFC, "_", include, "_logFC_corr_GXE_vs_sc")) +
##         ggtitle(paste0("logFC_", logFC, "_", include, "_logFCvszscore_corr_GXE_vs_sc")) +
##         ##
##         ##
##         theme(axis.text.x = element_text(angle = 90),
##                                   axis.title.x = element_blank(),
##                                   axis.title.y = element_blank(),
##                                   legend.position="none")
## ##
## figure <- ggarrange(sig.p + font("x.text", size = 14))
## # ncol = 2, nrow = 2)
## ##
## ##
## ##CHANGE FILE NAME
## ##
## ##
## #png(paste0("./3_treatEffect_3_GXE.outs/logFC_", logFC, "_", include, "_logFC_corr.png"), width=500, height=500, pointsize=16, res=125)
## png(paste0("./3_treatEffect_3_GXE.outs/logFC_", logFC, "_", include, "_logFCvszscore_corr.png"), width=500, height=500, pointsize=16, res=125)
## ##
## ##
## print(figure)
## dev.off()







## ## #----heatmap----
## ## mat <- scorr_data
## ## rownames(mat) <- mat$MCls
## ## mat <- mat %>% dplyr::select(!MCls)
## ## mat

## ## mat2 <- as.matrix(mat)
## ## mat2 <- mat2[nrow(mat2):1, ]
## ## mat2


## ## fig1 <- pheatmap(mat2, border_color="NA",
## ##                                   cluster_rows=F, cluster_cols=F,
## ##                                   #                    annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend=T,
## ##                                   show_colnames=T,
## ##                                   show_rownames=T,
## ##                                   na_col="white",
## ##                                   fontsize_col=20,
## ##                                   fontsize_row=10)


## ## figfn <- "./1_GXE.outs/logFC_corr_heatmap.png"
## ## png(figfn, width=1800, height=2100,res=225)
## ## print(fig1)
## ## dev.off()






#################################
# SCORR LOOP 3
#################################
#----scorr----
#----spearman_corr----
MCls <- unique(resDE$MCls)
MCls
comb.2 <- MCls

treats <- unique(resDE$contrast)
treats


#----all----
logFC <- "NA"
include <- "all"


#----Significance.0----
#logFC <- 0
#conditions <- c("both", "sc", "GXE")
#include <- "union"

#logFC <- 0
#conditions <- c("both")
#include <- "both"


#----Significance.25 ----
#logFC <- 0.25
#conditions <- c("both", "sc", "GXE")
#include <- "union"

#logFC <- 0.25
#conditions <- c("both")
#include <- "both"


#----Significance----
logFC <- 0.5
conditions <- c("both", "sc", "GXE")
include <- "union"

#logFC <- 0.5
#conditions <- c("both")
#include <- "both"


scorr_caff <- c()
scorr_nic <- c()
scorr_vitA <- c()
scorr_vitD <- c()
scorr_vitE <- c()
scorr_zinc <- c()
#scorr_water <- c()
scorr_caff_z <- c()
scorr_nic_z <- c()
scorr_vitA_z <- c()
scorr_vitD_z <- c()
scorr_vitE_z <- c()
scorr_zinc_z <- c()
#scorr_water_z <- c()


scorr_pval_caff <- c()
scorr_pval_nic <- c()
scorr_pval_vitA <- c()
scorr_pval_vitD <- c()
scorr_pval_vitE <- c()
scorr_pval_zinc <- c()
#scorr_pval_water <- c()
scorr_pval_caff_z <- c()
scorr_pval_nic_z <- c()
scorr_pval_vitA_z <- c()
scorr_pval_vitD_z <- c()
scorr_pval_vitE_z <- c()
scorr_pval_zinc_z <- c()
#scorr_pval_water_z <- c()



for (j in treats){
        resDE2 <- resDE %>% dplyr::filter(contrast==j)
        ##
        resGXE2 <- resGXE %>% dplyr::filter(contrast==j)
        ##
        # resGXE2.sum <- resGXE2 %>% arrange(padj)
        # resGXE2.sum <- resGXE2.sum %>% distinct(g.id, .keep_all = TRUE)
        ##
        resGXE2.sum <- resGXE2 %>% mutate(sig=ifelse(Significance.25 %in% c("Up", "Down"), "sig", "nonsig"))
        resGXE2.sum <- resGXE2.sum %>% arrange(desc(logFC))
        resGXE2.sum <- resGXE2.sum %>% arrange(padj)
        resGXE2.sum <- resGXE2.sum %>% arrange(desc(sig))
        resGXE2.sum <- resGXE2.sum %>% distinct(ensg, .keep_all = TRUE)
        ##
        resDE.GXE <- left_join(resDE2, resGXE2.sum, by="contrast.gene")
        resDE.GXE <- resDE.GXE %>% filter(!is.na(logFC))
        ##
        resDE.GXE <- resDE.GXE %>% mutate(Significance=ifelse(Significance.y %in% c("Up", "Down") &  Significance.x %in% c("Up", "Down"), "both", ifelse(Significance.y %in% c("Up", "Down") & (Significance.x %in% c("Not Significant", "NA") | is.na(Significance.x)), "GXE", ifelse((Significance.y %in% c("Not Significant", "NA") | is.na(Significance.y)) & Significance.x %in% c("Up", "Down"), "sc", "NA"))))
        ##
        resDE.GXE <- resDE.GXE %>% mutate(Significance.25=ifelse(Significance.25.y %in% c("Up", "Down") &  Significance.25.x %in% c("Up", "Down"), "both", ifelse(Significance.25.y %in% c("Up", "Down") & (Significance.25.x %in% c("Not Significant", "NA") | is.na(Significance.25.x)), "GXE", ifelse((Significance.25.y %in% c("Not Significant", "NA") | is.na(Significance.25.y)) & Significance.25.x %in% c("Up", "Down"), "sc", "NA"))))
        ##
        resDE.GXE <- resDE.GXE %>% mutate(Significance.0=ifelse(Significance.0.y %in% c("Up", "Down") &  Significance.0.x %in% c("Up", "Down"), "both", ifelse(Significance.0.y %in% c("Up", "Down") & (Significance.0.x %in% c("Not Significant", "NA") | is.na(Significance.0.x)), "GXE", ifelse((Significance.0.y %in% c("Not Significant", "NA") | is.na(Significance.0.y)) & Significance.0.x %in% c("Up", "Down"), "sc", "NA"))))
        ####################
        ## scorr_values ##
        ####################
        for (i in comb.2){
                res1 <- resDE.GXE %>% dplyr::filter(contrast.x==j, MCls==i)                                 %>% dplyr::filter(Significance %in% conditions)
                scorr <- as.numeric(cor.test(res1$logFC, res1$estimate, method = "spearman")$estimate)
                scorr_z <- as.numeric(cor.test(res1$logFC, res1$z_score, method = "spearman")$estimate)
                print(scorr)
                print(scorr_z)
                if (is.numeric(scorr) & j=="caffeine"){
                    scorr_caff <- c(scorr_caff, scorr)
                    scorr_caff_z <- c(scorr_caff_z, scorr_z)
                    } else if (is.numeric(scorr) & j=="nicotine"){
                    scorr_nic <- c(scorr_nic, scorr)
                    scorr_nic_z <- c(scorr_nic_z, scorr_z)
                    } else if (is.numeric(scorr) & j=="vitA"){
                    scorr_vitA <- c(scorr_vitA, scorr)
                    scorr_vitA_z <- c(scorr_vitA_z, scorr_z)
                    } else if (is.numeric(scorr) & j=="vitD"){
                    scorr_vitD <- c(scorr_vitD, scorr)
                    scorr_vitD_z <- c(scorr_vitD_z, scorr_z)
                    } else if (is.numeric(scorr) & j=="vitE"){
                    scorr_vitE <- c(scorr_vitE, scorr)
                    scorr_vitE_z <- c(scorr_vitE_z, scorr_z)
                    } else if (is.numeric(scorr) & j=="zinc"){
                    scorr_zinc <- c(scorr_zinc, scorr)
                    scorr_zinc_z <- c(scorr_zinc_z, scorr_z)
                    }
                    #else if (is.numeric(scorr) & j=="water"){
                    #scorr_water <- c(scorr_water, scorr)
                    #scorr_water_z <- c(scorr_water_z, scorr_z)
                    #}
                scorr_pval <- as.numeric(cor.test(res1$logFC, res1$estimate, method = "spearman")$p.value)
                scorr_pval_z <- as.numeric(cor.test(res1$logFC, res1$z_score, method = "spearman")$p.value)
                print(scorr_pval)
                print(scorr_pval_z)
                if (is.numeric(scorr_pval) & j=="caffeine"){
                    scorr_pval_caff <- c(scorr_pval_caff, scorr_pval)
                    scorr_pval_caff_z <- c(scorr_pval_caff_z, scorr_pval_z)
                    } else if (is.numeric(scorr_pval) & j=="nicotine"){
                    scorr_pval_nic <- c(scorr_pval_nic, scorr_pval)
                    scorr_pval_nic_z <- c(scorr_pval_nic_z, scorr_pval_z)
                    } else if (is.numeric(scorr_pval) & j=="vitA"){
                    scorr_pval_vitA <- c(scorr_pval_vitA, scorr_pval)
                    scorr_pval_vitA_z <- c(scorr_pval_vitA_z, scorr_pval_z)
                    } else if (is.numeric(scorr_pval) & j=="vitD"){
                    scorr_pval_vitD <- c(scorr_pval_vitD, scorr_pval)
                    scorr_pval_vitD_z <- c(scorr_pval_vitD_z, scorr_pval_z)
                    } else if (is.numeric(scorr_pval) & j=="vitE"){
                    scorr_pval_vitE <- c(scorr_pval_vitE, scorr_pval)
                    scorr_pval_vitE_z <- c(scorr_pval_vitE_z, scorr_pval_z)
                    } else if (is.numeric(scorr_pval) & j=="zinc"){
                    scorr_pval_zinc <- c(scorr_pval_zinc, scorr_pval)
                    scorr_pval_zinc_z <- c(scorr_pval_zinc_z, scorr_pval_z)
                    }
                    #else if (is.numeric(scorr_pval) & j=="water"){
                    #scorr_pval_water <- c(scorr_pval_water, scorr_pval)
                    #scorr_pval_water_z <- c(scorr_pval_water_z, scorr_pval_z)
                    #}       
        }
}


scorr_caff
scorr_nic
scorr_vitA 
scorr_vitD 
scorr_vitE 
scorr_zinc 
scorr_caff_z 
scorr_nic_z 
scorr_vitA_z 
scorr_vitD_z 
scorr_vitE_z 
scorr_zinc_z 
        

#################
## out of loop ##
#################  
#########################
## logFC vs logFC plot ##
#########################  
scorr_data <- data.frame(MCls, scorr_caff, scorr_nic,
                         scorr_vitA, scorr_vitD, scorr_vitE)
                         #scorr_water,
                         #scorr_zinc)
is.num <- sapply(scorr_data, is.numeric)
scorr_data[is.num] <- lapply(scorr_data[is.num], round, 2)
##
sig.melt <- melt(scorr_data)
##
sig.p <- ggplot(sig.melt, aes(x = variable, y = MCls)) +
    geom_tile(aes(fill=value), color="black")+
    scale_fill_gradient(low = "white", high = "white")+
    geom_text(aes(label = value), color = "black", size = 5) +
    ggtitle(paste0("logFC_", logFC, "_", include, "_logFC_corr_GXE_vs_sc")) +
    theme(axis.text.x = element_text(angle = 90),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position="none")

##
figure <- ggarrange(sig.p + font("x.text", size = 14))
##png(paste0("./3_treatEffect_3_GXE.outs/logFC_", logFC, "_", include, "_logFC_corr.estimate_contrastetOH.png"), width=700, height=700, pointsize=16, res=150)
png(paste0("./3_treatEffect_3_GXE.outs/logFC_", logFC, "_", include, "_logFC_corr.estimate_control.png"), width=700, height=700, pointsize=14, res=150)
print(figure)
dev.off()

#----pvalue----
scorr_pval_data <- data.frame(MCls, scorr_pval_caff, scorr_pval_nic,
                         scorr_pval_vitA, scorr_pval_vitD, scorr_pval_vitE)
                         #scorr_pval_water,
                         #scorr_pval_zinc)
is.num <- sapply(scorr_pval_data, is.numeric)
scorr_pval_data[is.num] <- lapply(scorr_pval_data[is.num], round, 2)
##
sig.melt <- melt(scorr_pval_data)
##
sig.p <- ggplot(sig.melt, aes(x = variable, y = MCls)) +
    geom_tile(aes(fill=value), color="black")+
    scale_fill_gradient(low = "white", high = "white")+
    geom_text(aes(label = value), color = "black", size = 5) +
    ggtitle(paste0("logFC_", logFC, "_", include, "_logFC_corr_GXE_vs_sc")) +
    theme(axis.text.x = element_text(angle = 90),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position="none")

##
figure <- ggarrange(sig.p + font("x.text", size = 14))
##png(paste0("./3_treatEffect_3_GXE.outs/logFC_", logFC, "_", include, "_logFC_corr.estimate_contrastetOH.png"), width=700, height=700, pointsize=16, res=150)
png(paste0("./3_treatEffect_3_GXE.outs/logFC_", logFC, "_", include, "_logFC_corr.pvalue_control.png"), width=700, height=700, pointsize=14, res=150)
print(figure)
dev.off()



#########################
## logFC vs zscore plot ##
#########################
scorr_data <- data.frame(MCls, scorr_caff_z, scorr_nic_z,
                         scorr_vitA_z, scorr_vitD_z,
                         scorr_vitE_z)
                         #scorr_water_z,
                         #scorr_zinc_z)
is.num <- sapply(scorr_data, is.numeric)
scorr_data[is.num] <- lapply(scorr_data[is.num], round, 2)
##
sig.melt <- melt(scorr_data)
##
sig.p <- ggplot(sig.melt, aes(x = variable, y = MCls)) +
    geom_tile(aes(fill=value), color="black")+
    scale_fill_gradient(low = "white", high = "white")+
    geom_text(aes(label = value), color = "black", size = 5) +
    ggtitle(paste0("logFC_", logFC, "_", include, "_logFCvszscore_corr_GXE_vs_sc")) +
    theme(axis.text.x = element_text(angle = 90),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position="none")

##
figure <- ggarrange(sig.p + font("x.text", size = 14))
png(paste0("./3_treatEffect_3_GXE.outs/logFC_", logFC, "_", include, "_logFCvszscore_corr.estimate_control.png"), width=700, height=700, pointsize=14, res=150)
print(figure)
dev.off()



#----pvalue----
scorr_pval_data <- data.frame(MCls, scorr_pval_caff_z, scorr_pval_nic_z,
                         scorr_pval_vitA_z, scorr_pval_vitD_z,
                         scorr_pval_vitE_z)
                         #scorr_pval_water_z,
                         #scorr_pval_zinc_z)
is.num <- sapply(scorr_pval_data, is.numeric)
scorr_pval_data[is.num] <- lapply(scorr_pval_data[is.num], round, 2)
##
sig.melt <- melt(scorr_pval_data)
##
sig.p <- ggplot(sig.melt, aes(x = variable, y = MCls)) +
    geom_tile(aes(fill=value), color="black")+
    scale_fill_gradient(low = "white", high = "white")+
    geom_text(aes(label = value), color = "black", size = 5) +
    ggtitle(paste0("logFC_", logFC, "_", include, "_logFCvszscore_corr_GXE_vs_sc")) +
    theme(axis.text.x = element_text(angle = 90),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position="none")

##
figure <- ggarrange(sig.p + font("x.text", size = 14))
png(paste0("./3_treatEffect_3_GXE.outs/logFC_", logFC, "_", include, "_logFCvszscore_corr.pvalue_control.png"), width=700, height=700, pointsize=14, res=150)
print(figure)
dev.off()














##############################################################
# SCORR LOOP 3 - for all treatments(sc) vs all treatments(GXE)
##############################################################
#----scorr----
#----spearman_corr----

head(resDE)

resGXE <- resGXE %>% dplyr::rename(gene = g.id)
#head(resGXE)


#MCls <- unique(resDE$MCls)
MCls <- c("0-CD4Naive", "6-Monocyte")
MCls
#comb.2 <- MCls


treats.DE <- unique(resDE$contrast)
treats.DE

treats <- unique(resGXE$contrast)
treats


#----all----
logFC <- "NA"
include <- "all"


#----Significance.0----
logFC <- 0
conditions <- c("both", "sc", "GXE")
include <- "union"

#logFC <- 0
#conditions <- c("both")
#include <- "both"


#----Significance.25 ----
#logFC <- 0.25
#conditions <- c("both", "sc", "GXE")
#include <- "union"

#logFC <- 0.25
#conditions <- c("both")
#include <- "both"


#----Significance----
#logFC <- 0.5
#conditions <- c("both", "sc", "GXE")
#include <- "union"

#logFC <- 0.5
#conditions <- c("both")
#include <- "both"


       

#oneMCl = "0-CD4Naive"
#oneMCl = "6-Monocyte"

for (i in MCls){
print(i)    
scorr_caff <- c()
scorr_nic <- c()
scorr_vitA <- c()
scorr_vitD <- c()
scorr_vitE <- c()
scorr_water <- c()
scorr_zinc <- c()
scorr_caff_z <- c()
scorr_nic_z <- c()
scorr_vitA_z <- c()
scorr_vitD_z <- c()
scorr_vitE_z <- c()
scorr_water_z <- c()
scorr_zinc_z <- c()
for (j in treats.DE){
    print(j)
    resDE2 <- resDE %>% dplyr::filter(MCls==i, contrast==j)
    for (k in treats){
        print(k)
        resGXE2 <- resGXE %>% dplyr::filter(contrast==k)
        resGXE2.sum <- resGXE2 %>% mutate(sig=ifelse(Significance.25 %in% c("Up", "Down"), "sig", "nonsig"))
        resGXE2.sum <- resGXE2.sum %>% arrange(desc(logFC))
        resGXE2.sum <- resGXE2.sum %>% arrange(padj)
        resGXE2.sum <- resGXE2.sum %>% arrange(desc(sig))
        resGXE2.sum <- resGXE2.sum %>% distinct(ensg, .keep_all = TRUE)
        ##
        resDE.GXE <- left_join(resDE2, resGXE2.sum, by="gene")
        resDE.GXE <- resDE.GXE %>% filter(!is.na(logFC))
        ##
        ##resDE.GXE <- resDE.GXE %>% mutate(Significance=ifelse(Significance.y %in% c("Up", "Down") &  Significance.x %in% c("Up", "Down"), "both", ifelse(Significance.y %in% c("Up", "Down") & (Significance.x %in% c("Not Significant", "NA") | is.na(Significance.x)), "GXE", ifelse((Significance.y %in% c("Not Significant", "NA") | is.na(Significance.y)) & Significance.x %in% c("Up", "Down"), "sc", "NA"))))
        ##resDE.GXE <- resDE.GXE %>% mutate(Significance.25=ifelse(Significance.25.y %in% c("Up", "Down") &  Significance.25.x %in% c("Up", "Down"), "both", ifelse(Significance.25.y %in% c("Up", "Down") & (Significance.25.x %in% c("Not Significant", "NA") | is.na(Significance.25.x)), "GXE", ifelse((Significance.25.y %in% c("Not Significant", "NA") | is.na(Significance.25.y)) & Significance.25.x %in% c("Up", "Down"), "sc", "NA"))))
        ##resDE.GXE <- resDE.GXE %>% mutate(Significance.0=ifelse(Significance.0.y %in% c("Up", "Down") &  Significance.0.x %in% c("Up", "Down"), "both", ifelse(Significance.0.y %in% c("Up", "Down") & (Significance.0.x %in% c("Not Significant", "NA") | is.na(Significance.0.x)), "GXE", ifelse((Significance.0.y %in% c("Not Significant", "NA") | is.na(Significance.0.y)) & Significance.0.x %in% c("Up", "Down"), "sc", "NA"))))
        ##
        ##resDE.GXE <- resDE.GXE %>% dplyr::filter(!is.na(Significance))
        ##head(resDE.GXE)
        ##nrow(resDE.GXE)
        ##table(resDE.GXE$Significance)
        ##
        ####################
        ## scorr_values ##
        ####################
        ##for (i in comb.2){
        tryCatch({
            logFC <- "NA"
            include <- "all"
            ##res1 <- resDE.GXE %>% dplyr::filter(contrast.x==j, MCls==i)
            res1 <- resDE.GXE
            scorr <- as.numeric(cor.test(res1$logFC, res1$estimate, method = "spearman")$estimate)
            print(scorr)
            if (is.numeric(scorr) & j=="caffeine"){
                scorr_caff <- c(scorr_caff, scorr)
            }  else if (is.numeric(scorr) & j=="nicotine"){
                scorr_nic <- c(scorr_nic, scorr)
            } else if (is.numeric(scorr) & j=="vitA"){
                scorr_vitA <- c(scorr_vitA, scorr)
            }  else if (is.numeric(scorr) & j=="vitD"){
                scorr_vitD <- c(scorr_vitD, scorr)
            } else if (is.numeric(scorr) & j=="vitE"){
                scorr_vitE <- c(scorr_vitE, scorr)
            } else if (is.numeric(scorr) & j=="water"){
                scorr_water <- c(scorr_water, scorr)
            } else if (is.numeric(scorr) & j=="zinc"){
                scorr_zinc <- c(scorr_zinc, scorr)
            }
            scorr_z <- as.numeric(cor.test(res1$logFC, res1$z_score, method = "spearman")$estimate)
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
scorr_caff
scorr_nic
scorr_vitA
scorr_vitD
scorr_vitE 
scorr_water
scorr_zinc
scorr_caff_z
scorr_vitA_z
scorr_zinc_z
scorr_nic_z
scorr_vitD_z
scorr_vitE_z
scorr_water_z
#################
## out of loop ##
#################  
#########################
## logFC vs logFC plot ##
#########################  
scorr_data <- data.frame(treats, scorr_caff, scorr_nic,
                         scorr_vitA, scorr_vitD, scorr_vitE,
                         scorr_water, scorr_zinc)
is.num <- sapply(scorr_data, is.numeric)
scorr_data[is.num] <- lapply(scorr_data[is.num], round, 2)
##
sig.melt <- melt(scorr_data)
##
##sig.p <- ggplot(sig.melt, aes(x = variable, y = MCls)) +
##head(sig.melt)
sig.p <- ggplot(sig.melt, aes(x = treats, y = variable)) +
    geom_tile(aes(fill=value), color="black")+
    scale_fill_gradient(low = "white", high = "white")+
    geom_text(aes(label = value), color = "black", size = 5) +
    ggtitle(paste0("logFC_", logFC, "_", include, "_logFC_corr_GXE_vs_sc")) +
    theme(axis.text.x = element_text(angle = 90),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position="none")
##
figure <- ggarrange(sig.p + font("x.text", size = 14))
png(paste0("./3_treatEffect_3_GXE.outs/logFC_", logFC, "_", include, "_logFC_corr.", i, ".alltreats.includewater.png"), width=600, height=600, pointsize=16, res=125)
#png(paste0("./3_treatEffect_3_GXE.outs/logFC_", logFC, "_", include, "_logFC_corr.", i, ".alltreats.includewater.pvalue.png"), width=500, height=500, pointsize=16, res=125)
print(figure)
dev.off()
#########################
## logFC vs zscore plot ##
#########################
scorr_data <- data.frame(treats, scorr_caff_z, scorr_nic_z,
                         scorr_vitA_z, scorr_vitD_z,
                         scorr_vitE_z,
                         scorr_water_z, scorr_zinc_z)
is.num <- sapply(scorr_data, is.numeric)
scorr_data[is.num] <- lapply(scorr_data[is.num], round, 2)
##
sig.melt <- melt(scorr_data)
##
##sig.p <- ggplot(sig.melt, aes(x = variable, y = MCls)) +
sig.p <- ggplot(sig.melt, aes(x = treats, y = variable)) +
    geom_tile(aes(fill=value), color="black")+
    scale_fill_gradient(low = "white", high = "white")+
    geom_text(aes(label = value), color = "black", size = 5) +
    ggtitle(paste0("logFC_", logFC, "_", include, "_logFCvszscore_corr_GXE_vs_sc")) +
    theme(axis.text.x = element_text(angle = 90),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position="none")


##
figure <- ggarrange(sig.p + font("x.text", size = 14))
png(paste0("./3_treatEffect_3_GXE.outs/logFC_", logFC, "_", include, "_logFCvszscore_corr.", i, ".alltreats.includewater.png"), width=600, height=600, pointsize=16, res=125)
#png(paste0("./3_treatEffect_3_GXE.outs/logFC_", logFC, "_", include, "_logFCvszscore_corr.", i, ".alltreats.includewater.pvalue.png"), width=500, height=500, pointsize=16, res=125)
print(figure)
dev.off()
}      
