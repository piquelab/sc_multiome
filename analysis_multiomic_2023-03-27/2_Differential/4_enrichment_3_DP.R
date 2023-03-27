##
library(tidyverse)
library(annotables)
library(clusterProfiler)
library(org.Hs.eg.db)
##
library(EnsDb.Hsapiens.v86)
#library(EnsDb.Hsapiens.v75)
library(ChIPseeker,lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
##
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(ggExtra)
library(RColorBrewer)
library(gtable)
library(ggsignif)
library(pheatmap)
library(corrplot)
library(viridis)
theme_set(theme_grey())

rm(list=ls())

###
outdir <- "./4_enrichment_2.outs/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)



##########################
## resDP                ##
##########################
opfn <- paste0("./3_treatEffect_2_2.outs/resDP_control.rds")
#write_rds(resDP, opfn)
resDP <- read_rds(opfn)



## #----peakAnno----
peakAnno <- read_rds("./2.2_compareRNAandATAC.outs/1_annot.signac_macs2_0.1_cn.rds")
##head(peakAnno)

peakAnno <- peakAnno %>% mutate(query_region=gsub(":", "-", query_region))
peakAnno <- peakAnno %>% mutate(dis=abs(distance))
colnames(peakAnno)[which(names(peakAnno) == "query_region")] <- "peak"
##peakAnno <- peakAnno %>% dplyr::select(peak, gene_id, distance, dis)# flank_geneIds) 
head(peakAnno)
nrow(peakAnno)




##----background gene----
#gene <- unique(unlist(str_split(peakAnno$flank_geneIds, ";")))
gene <- unique(peakAnno$gene_id)

head(grch38)
x <- grch38%>%dplyr::filter(ensgene%in%gene, grepl("protein_coding", biotype))
head(x)
nrow(x)


geneBG <- bitr(x$ensgene, fromType="ENSEMBL", toType=c("ENTREZID","SYMBOL"),
                              OrgDb=org.Hs.eg.db)
head(geneBG)
nrow(geneBG)




##----gene cluster for enrichment data analysis----
head(resDP)
nrow(resDP)

res2 <- resDP%>%as.data.frame()%>%
    ##dplyr::rename("peak_region"="gene")%>%
    mutate(comb=paste(MCls, contrast, sep="_"), direction=ifelse(estimate>0, "Up", "Down"))

res2 <- res2 %>% dplyr::filter(abs(estimate)>0.5, p.adjusted<0.1)
##head(res2)
##nrow(res2)

res2 <- res2%>%left_join(peakAnno, by="peak")
head(res2)
nrow(res2)

##
## comb <- sort(unique(res2$comb))
## df0 <- map_dfr(comb, function(ii){
##    tmp <- res2%>%dplyr::filter(comb==ii)%>%mutate(dis=abs(distanceToTSS))%>%
##       group_by(geneId)%>%top_n(-1, dis)
##       dplyr::filter(dis==min(dis,na.rm=T))%>%ungroup()
##    tmp
## })

##df_trial <- res2%>%group_by(contrast, MCls, gene_id) %>%top_n(-1, dis)
##head(df_trial)
##nrow(df_trial)

df0 <- res2%>%group_by(contrast, MCls, gene_id)%>%top_n(-1, dis)%>%ungroup()%>%as.data.frame()
head(df0)
nrow(df0)

df2 <- df0%>%dplyr::filter(dis<100000)
head(df2)
nrow(df2)

x <- bitr(df2$gene_id, fromType="ENSEMBL",
                    toType=c("ENTREZID","SYMBOL"), OrgDb=org.Hs.eg.db)
head(x)
nrow(x)

geneCluster <- df2%>%inner_join(x, by=c("gene_id"="ENSEMBL"))
head(geneCluster)
nrow(geneCluster)








#######################
##   GO enrichment   ##
#######################
# "GO analysis for direction (up and down) directly", "\n")
cgDP <- compareCluster(ENTREZID~contrast+MCls+direction,
                         data=geneCluster,
                         universe=geneBG$ENTREZID,
                         fun="enrichGO",
                         OrgDb="org.Hs.eg.db",
                         pvalueCutoff=1,
                         qvalueCutoff=1,
                         ont="ALL",
                         minGSSize=0,
                         maxGSSize=1000)

cgDP <- cg.DP

head(cgDP)
nrow(cgDP)

write_rds(cgDP, "./4_enrichment_2.outs/1_enrichGO_universe.rds")

#cgDP <- readRDS("./4_enrichment_2.outs/1_enrichGO.rds")
#cgDP <- readRDS("./4_enrichment_2.outs/1_enrichGO_universe.rds")




####################
### show figures ###
####################
cgDP.df <- as.data.frame(cgDP)
##cgDP <- cgDP%>%mutate(Cluster1=gsub("\\.((Down)|(Up))", "", Cluster))

colnames(cgDP.df)

head(cgDP.df$Cluster)

MCls <- unique(cgDP.df$MCls)
#MCls
contrast <- unique(cgDP.df$contrast)
#contrast



##----celltype----
cl <- paste(rep(contrast, times=length(unique(cgDP.df$MCls))),
                        rep(MCls, each=length(unique(cgDP.df$contrast))), sep=".")
cl <- paste(rep(cl,times=2), rep(c("Up","Down"), each=48), sep=".")
cl


ii <- paste(rep(c("A","B","C","D", "E", "F", "G", "H"),each=6),rep(1:6,times=8), sep="")
cl2 <- paste(rep(c("X","Y"),each=48), rep(ii,times=2), sep=".")
cl2

cluster2 <- setNames(cl2, cl)
head(cluster2)

##lab2 <- setNames(gsub("-","+",cl),cl2)
lab2 <- setNames(cl,cl2)
head(lab2)

cg0 <- as.data.frame(cgDP)
head(cg0)
head(cg0$Cluster)

x <- cluster2[as.character(cg0$Cluster)]
head(x)

cgDP2<- cgDP%>%
        mutate(ClusterNew=x,
                          maxGSSize=as.numeric(gsub("/.*", "", BgRatio)),
                          ngene=as.numeric(gsub(".*/", "", GeneRatio)))
head(cgDP2 %>% dplyr::select(!geneID))




##----treatment----
cl <- paste(rep(contrast, times=length(unique(cgDP.df$MCls))),
                        rep(MCls, each=length(unique(cgDP.df$contrast))), sep=".")
cl <- paste(rep(cl,times=2), rep(c("Up","Down"), each=48), sep=".")
cl


ii <- paste(rep(c("A","B","C","D","E","F"),times=8),rep(1:8,each=6), sep="")
cl2 <- paste(rep(c("X","Y"),each=48), rep(ii,times=2), sep=".")
cl2

cluster2 <- setNames(cl2, cl)
head(cluster2)

##lab2 <- setNames(gsub("-","+",cl),cl2)
lab2 <- setNames(cl,cl2)
head(lab2)

cg0 <- as.data.frame(cgDP)
#head(cg0)
head(cg0$Cluster)

x <- cluster2[as.character(cg0$Cluster)]
head(x)

cgDP2<- cgDP%>%
        mutate(ClusterNew=x,
               maxGSSize=as.numeric(gsub("/.*", "", BgRatio)),
               ngene=as.numeric(gsub(".*/", "", GeneRatio)))
head(cgDP2 %>% dplyr::select(!geneID))





####################
### plot figures ###
####################
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
require(DOSE)


#sig_cgDP2 <- cgDP2 %>% dplyr::filter(p.adjust<0.1)
#nrow(sig_cgDP2)


sig_cgDP2 <- cgDP2 %>% dplyr::filter(maxGSSize>10, maxGSSize<500, ngene>10, p.adjust<0.1)
head(sig_cgDP2)
nrow(sig_cgDP2)


fig1 <- enrichplot::dotplot(sig_cgDP2, x="ClusterNew", showCategory=5)+
        scale_x_discrete("", labels=lab2)+
        theme(axis.text.x=element_text(angle=90, hjust=1,size=10),
                        axis.text.y=element_text(size=10))
#fig1 <- dotplot(cg2, showcategory=5, split=".sign") + facet_grid(.~.sign)
#fig1 <- dotplot(cg2, x="ClusterNew", showCategory=5)


###
figfn <- "./4_enrichment_2.outs/Figure1.1_GO_DP_ClusterNew_treat_universe_highres.png"
png(figfn, width=4000, height=6000, res=175)
print(fig1)
dev.off()






##
fig1 <- cnetplot(sig_cgDP2)



end








#####################
### KEGG analysis ###
#####################
ck <- compareCluster(ENTREZID~contrast+MCls+direction,
                         data=geneCluster,
                         universe=geneBG$ENTREZID,
                         fun="enrichKEGG",
                         ## OrgDb="org.Hs.eg.db",
                         pvalueCutoff=1,
                         qvalueCutoff=1,
                         ## ont="ALL",
                         minGSSize=0,
                         maxGSSize=1000)

ck <- ck%>%mutate(Cluster1=gsub("\\.((Down)|(Up))", "", Cluster))
write_rds(ck,"./4_enrichment.outs/2_enrichKEGG.rds")


###
### figures

ck <- read_rds("./4_enrichment.outs/2_enrichKEGG.rds")

cl <- paste( rep(c("LPS", "PHA", "LPS-DEX", "PHA-DEX"),each=4),
               rep(c("Bcell", "Monocyte", "NKcell", "Tcell"), times=4), sep=".")
cl <- paste(rep(cl,times=2), rep(c("Up","Down"), each=16), sep=".")

ii <- paste(rep(c("A","B","C","D"),each=4),rep(1:4,times=4), sep="")
cl2 <- paste(rep(c("X","Y"),each=16), rep(ii,times=2), sep=".")
cluster2 <- setNames(cl2, cl)
lab2 <- setNames(gsub("-","+",cl),cl2)

##
ck0 <- as.data.frame(ck)
x <- cluster2[as.character(ck0$Cluster)]
ck2 <- ck%>%
       mutate(ClusterNew=x,
                        maxGSSize=as.numeric(gsub("/.*", "", BgRatio)),
                        ngene=as.numeric(gsub(".*/", "", GeneRatio)))%>%
       dplyr::filter(maxGSSize>10, maxGSSize<500, ngene>10, p.adjust<0.1)

fig2 <- enrichplot::dotplot(ck2, x="ClusterNew", showCategory=5)+
       scale_x_discrete("", labels=lab2)+
       theme(axis.text.x=element_text(angle=60, hjust=1,size=12),
                      axis.text.y=element_text(size=10))

###
figfn <- "./4_enrichment.outs/Figure2.1_KEGG.png"
png(figfn, width=2000, height=1500, res=150)
print(fig2)
dev.off()

##
figfn <- "./4_enrichment.outs/Figure2.1_KEGG.pdf"
pdf(figfn, width=20, height=15)
print(fig2)
dev.off()













## #################
## ### examples  ###
## #################

## ###
## odds.fun <-  function(df){
##     ###
##     res <- map_dfr(1:nrow(df), function(i){
##         Diff <- as.numeric(df[i, c("Diff.in", "Diff.not")])
##         Bg <- as.numeric(df[i, c("Bg.in", "Bg.not")])
##         dat <- data.frame(Diff=Diff, Bg=Bg)
##         rownames(dat) <- c("in.category", "not.category")
##         fish <- fisher.test(dat)
##         res0 <- data.frame(odds=as.numeric(fish$estimate),
##                            down=fish$conf.int[1],
##                            up=fish$conf.int[2])
##         res0
##     })
##     ###
##     df$odds <- res$odds
##     df$down <- res$down
##     df$up <- res$up
##     df
## }


## ExampleGOplot <- function(cg){

##     ### prepare data
##        ## x <- str_split(cg$GeneRatio, "/", simplify=T)
##        ## GeneRatio <- as.numeric(x[,1])/as.numeric(x[,2])
##        Drt2 <- c("Up"=1, "Down"=2)
##        cg <- cg%>%mutate(Direction2=Drt2[direction.x],
##                                contrast2=paste(Direction2, contrast.x, sep="."))%>%
##                  mutate(contrast2=gsub("-", "+", contrast2))
##        ## cg$size <- rep(1,nrow(cg))
##        ## cg$size[GeneRatio>=0.05&GeneRatio<0.15] <- 2
##        ## cg$size[GeneRatio>=0.15] <- 3
##        #
##        cg <- cg%>%drop_na(odds)
##        fig0 <- ggplot(cg, aes(x=contrast2, y=MCls.x))+
##                  geom_point(aes(size=odds, colour=p2))+
##                  scale_x_discrete(labels=c("1.LPS"="LPS.Up", "2.LPS"="LPS.Down",
##                                            "1.LPS+DEX"="LPS+DEX.Up", "2.LPS+DEX"="LPS+DEX.Down",
##                                            "1.PHA"="PHA.Up", "2.PHA"="PHA.Down",
##                                            "1.PHA+DEX"="PHA+DEX.Up", "2.PHA+DEX"="PHA+DEX.Down"))+
##            scale_colour_gradient(name="p.adjust",
##                                  low="blue", high="red", na.value=NA, trans="reverse", n.breaks=5,
##                                  guide=guide_colourbar(order=1))+    #"#ffa500"
##            scale_size_binned("odds ratio",
##                              guide=guide_bins(show.limits=TRUE, axis=TRUE,
##                                               axis.show=arrow(length=unit(1.5,"mm"), ends="both"), order=2),
##                              n.breaks=4)+
##            theme_bw()+
##            theme(axis.title=element_blank(),
##                  ## axis.text.y=element_text(size=12),
##                  legend.background=element_blank(),
##                  ## legend.title=element_text(size=8),
##                  ## legend.text=element_text(size=6),
##                  legend.key.size=grid::unit(0.6, "lines"))
##     fig0
## }



## ###

## contrast <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")
## MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
## rn <- paste(rep(contrast, each=8), rep(rep(MCls, each=2), times=4),
##                         rep(rep(c("Down", "Up"),times=4), times=4), sep=".")
## tmp <- data.frame(contrast=rep(contrast, each=8),
##                      MCls=rep(rep(MCls, each=2), times=4),
##                      direction=rep(rep(c("Down", "Up"),times=4), times=4))%>%
##        mutate(rn=paste(contrast, MCls, direction, sep="."))


## ###
## ###
## cg <- read_rds("./4_enrichment.outs/1_enrichGO.rds")

## ###
## ### response to LPS
## cg2 <- cg%>%dplyr::filter(Description=="response to lipopolysaccharide")

## cg2 <- cg2%>%as.data.frame()%>%
##         mutate(Diff.in=as.numeric(gsub("/.*","",GeneRatio)),
##                           Diff.total=as.numeric(gsub(".*/","",GeneRatio)),
##                           Diff.not=Diff.total-Diff.in)
## cg2 <- cg2%>%
##         mutate(Bg.in=as.numeric(gsub("/.*","", BgRatio)),
##                           Bg.total=as.numeric(gsub(".*/","", BgRatio)),
##                           Bg.not=Bg.total-Bg.in)

## cg2 <- odds.fun(cg2)
## cg2 <- cg2%>%full_join(tmp, by=c("Cluster"="rn"))
## cg2 <- cg2%>%mutate(p2=ifelse(p.adjust>0.2, NA, p.adjust))

## fig <- ExampleGOplot(cg2)+
##        ggtitle("response to lipopolysaccharide")+
##        theme(axis.text.x=element_text(angle=-90, size=8, hjust=0, vjust=0.5),
##                       plot.title=element_text(hjust=0.5, size=12))

## figfn <-"./4_enrichment.outs/Figure3.1_response_to_LPS.png"
## png(figfn, width=550, height=400, res=120)
## print(fig)
## dev.off()


## ###
## ### type I interferon
## cg2 <- cg%>%dplyr::filter(Description=="type I interferon signaling pathway")

## cg2 <- cg2%>%as.data.frame()%>%
##         mutate(Diff.in=as.numeric(gsub("/.*","",GeneRatio)),
##                           Diff.total=as.numeric(gsub(".*/","",GeneRatio)),
##                           Diff.not=Diff.total-Diff.in)
## cg2 <- cg2%>%
##         mutate(Bg.in=as.numeric(gsub("/.*","", BgRatio)),
##                           Bg.total=as.numeric(gsub(".*/","", BgRatio)),
##                           Bg.not=Bg.total-Bg.in)

## cg2 <- odds.fun(cg2)
## cg2 <- cg2%>%full_join(tmp, by=c("Cluster"="rn"))
## cg2 <- cg2%>%mutate(p2=ifelse(p.adjust>0.2, NA, p.adjust))

## fig <- ExampleGOplot(cg2)+
##        ggtitle("Type I IFN signaling")+
##        theme(axis.text.x=element_text(angle=-90, size=8, hjust=0, vjust=0.5),
##                       plot.title=element_text(hjust=0.5, size=12))

## figfn <-"./4_enrichment.outs/Figure3.2_IFN.png"
## png(figfn, width=550, height=400, res=120)
## print(fig)
## dev.off()


## ###
## ###


