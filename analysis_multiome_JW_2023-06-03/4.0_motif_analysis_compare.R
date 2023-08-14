##
library(Matrix)
library(tidyverse)
library(data.table)
library(irlba)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(Signac)
library(SeuratWrappers)
library(SeuratObject)
library(qvalue)


##
library(cowplot)
library(RColorBrewer)
library(ComplexHeatmap)
library(viridis)
library(circlize)
library(ggrepel)
library(ggrastr)
library(openxlsx)

rm(list=ls())

outdir <- "./4.0_motif_compare.outs/"
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)

###
### 3.2_motif_diff_fitInd_JW.rds I generated using Mhmd's script considering individual
### in the differetial analysis and the final used for genetic analysis
### 3.3_motif_diff_JW.rds I generated using Mhmd's script not considering individual
### 2_motif.ave_macs2_0.1_cn_control.rds average motif activity across cells per combinations




###
### motif name
fn <- "./sc_multiome_data/3_motif/2_motif.activities.outs/1_scATAC.motifActivities.rds"
atac <- read_rds(fn)
x <- Motifs(atac)
motifs <- unlist(x@motif.names)
## motifs_name <- unlist(unname(motifs))
 
########################################################
### compare the two results files for response motif ###
########################################################

###
### response motifs recent file
fn <- "./sc_multiome_data/3_motif/2_motif.activities.outs/3_motif.diff.results_macs2_0.1_cn_indfilt_incsampleID_all_cols_control.rds"
res <- read_rds(fn)%>%mutate(comb=paste(MCls, contrast, sep="_"), zscore=beta/stderr)
res <- res%>%mutate(rn_id=paste(comb, gene, sep="_"))%>%
   select(rn_id, comb, MCls, contrast, zscore)

comb_set <- sort(unique(res$comb))


###
res2 <- map_dfr(comb_set, function(ii){
   ###
   tmp <- res%>%dplyr::filter(comb==ii)
   th0 <- quantile(abs(tmp$beta), probs=0.9)
   tmp2 <- tmp%>%filter(qval<0.1, abs(beta)>th0)%>%dplyr::select(motif_ID=gene, motif_name=motif, comb)
   tmp2
})
summ <- res2%>%group_by(comb)%>%summarize(ny=n(), .groups="drop")%>%as.data.frame()
    




###
### file I used for genetic analysis
fn2 <- "./sc_multiome_data/3_motif/2_motif.activities.outs/3_motif.diff.results_macs2_0.1_cn.rds"
res2 <- read_rds(fn2)%>%mutate(comb=paste(MCls, contrast, sep="_"), zscore=beta/stderr)
res2 <- res2%>%filter(comb%in%comb_set)%>%
    mutate(rn_id=paste(comb, gene, sep="_"))%>%
    select(rn_id, zscore_old=zscore)


###
### plot data
res_comb <- res%>%inner_join(res2, by="rn_id")


### colors
col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
  "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
   "6_Monocyte"="#984ea3", "7_dnT"="black")
col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
       "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3")


### plots
p <- ggplot(res_comb, aes(x=zscore_old, y=zscore, color=contrast))+
   geom_point(size=0.2)+
   facet_wrap(~MCls, ncol=4, scales="free")+
   scale_color_manual(values=col2, guide=guide_legend(override.aes=list(size=2)))+
   geom_abline(color="grey")+ 
   xlab(bquote(italic(Z)~"-score from old file for genetic analysis"))+
   ylab(bquote(italic(Z)~"-score from new file"))+
   theme_bw()+
   theme(legend.title=element_blank(),
         legend.text=element_text(size=9),
         legend.key.size=unit(0.5, "cm"),
         axis.title=element_text(size=10),
         axis.text=element_text(size=10))

figfn <- paste(outdir, "Figure1_compare_zscore.png", sep="")
ggsave(figfn, p, device="png", width=920, height=420, units="px", dpi=120)



################################################################
### use the Mohammed script to perform differential analysis ###
################################################################

fn <- "./4.0_motif_compare.outs/2_motif.ave_macs2_0.1_cn_control.rds"
mat <- read_rds(fn) 

bti <-colnames(mat)
cvt0 <- str_split(bti, "_", simplify=T)%>%as.data.frame()

cvt <- data.frame(bti=bti, MCls=paste0(cvt0[,1], "_", cvt0[,2]), treat=cvt0[,3], sampleID=cvt0[,4])

MCls <- unique(cvt$MCls)
MCls





# filtered
#-------------------------------------------------------------------------
cvt <- cvt%>%mutate(MCls_IDs=paste(cvt$MCls, ".", cvt$sampleID, sep=""))


MCls.IDs <- unique(cvt$MCls_IDs)
#MCls.IDs

cvt.filtered <- map_dfr(MCls.IDs, function(oneX){
   ###
   cat(oneX, "\n") 
   cvt.2 <- cvt%>%
       dplyr::filter(MCls==str_split(oneX, "\\.", simplify=T)[1],
                    sampleID==str_split(oneX, "\\.", simplify=T)[2])
   ##
   cvt.2 <- cvt.2%>%mutate(count=ifelse(nrow(cvt.2)>=7,1,0))
   cvt.2
 })
#cvt.filtered

cvt.2 <- cvt.filtered %>% dplyr::filter(count==1)



MCls <- unique(cvt.2$MCls)

 
### myDE #adding Z
myDE <- function(y, X, Z, gene){
    ##
    con.ls <- list("caffeine"=c("control", "caffeine"),
                   "nicotine"=c("control", "nicotine"),
                   "vitA"=c("control", "vitA"),
                   "vitD"=c("control", "vitD"),
                   "vitE"=c("control", "vitE"),
                   "zinc"=c("control", "zinc"))
    ##
    ## con.ls <- list("caffeine"=c("water", "caffeine"),
    ##               "nicotine"=c("water", "nicotine"),
    ##               "vitA"=c("etOH", "vitA"),
    ##               "vitD"=c("etOH", "vitD"),
    ##               "vitE"=c("etOH", "vitE"),
    ##               "water"=c("etOH", "water"),
    ##               "zinc"=c("water", "zinc"))
    ##
    ## con.ls <- list("LPS"=c("CTRL","LPS"),
    ##    "LPS-DEX"=c("LPS","LPS-DEX"),
    ##    "PHA"=c("CTRL", "PHA"),
    ##    "PHA-DEX"=c("PHA","PHA-DEX"))
    ## Contrast <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")
    ##Contrast <- c("caffeine", "nicotine", "vitA", "vitD",
    ##             "vitE", "water", "zinc") 
    Contrast <- c("caffeine", "nicotine", "vitA", "vitD",
                 "vitE", "zinc") 
    ##
    ## linear regression
    x1 <- X[,1]
    z1 <- Z[,1]
    lm0 <- try(lm(y~0+x1+z1), silent=T)
    ##
    if ( (class(lm0)!="try-error")){
        b <- coef(lm0)
        nn <- gsub("x[12]", "", names(b))
        names(b) <- nn
        vb <- diag(vcov(lm0))
        names(vb) <- nn
        ##
        ## Contrast
        dd <- lapply(Contrast,function(one){
            con0 <- con.ls[[one]]
            if ( all(con0%in%nn)){
                bhat <- b[con0[2]]-b[con0[1]]
                sdhat <- sqrt(vb[con0[2]]+vb[con0[1]])
                z <- bhat/sdhat
                p <- 2*pnorm(-abs(z))
                dd1 <- data.frame(gene=gene, beta=bhat, stderr=sdhat,
                                  pval=p, contrast=one,
                                  coeff_treat=b[con0[2]], coeff_con=b[con0[1]],
                                  vb_treat=vb[con0[2]], vb_con=vb[con0[1]])
            }else{
                dd1 <- NA                                                                       }
            dd1
        })
        dd <- dd[!is.na(dd)]
        dd <- do.call(rbind,dd)
        ##
    }else{
        dd <- NA
    } ## End if try-error
    dd
}
    




#----function----
#change cvt either cvt or cvt.2
res <- map_dfr(MCls, function(oneMCl){
    ##
    cvti <- cvt.2%>%dplyr::filter(MCls==oneMCl)
    ##print(head(cvti))
    mati <- mat[,cvti$bti]
    ##print(head(mati))
    X <- data.frame(x1=cvti$treat)
    ##print(head(X))
    ##adding this line
    Z <- data.frame(x1=cvti$sampleID)
    ##print(head(Z))
    ##
    rn <- rownames(mati)
    ##print(head(rn))
    TMP <- lapply(rn, function(ii){
        y <- mati[ii,]
        dd <- myDE(y, X, Z, ii) #adding Z in this
        dd
    }) ##,mc.cores=1)
    TMP <- TMP[!is.na(TMP)]
    TMP <- do.call(rbind, TMP)%>%as.data.frame()%>%mutate(MCls=oneMCl)
    TMP 
})    



 
### add qvalue
res2 <- res%>%group_by(MCls, contrast)%>%
   mutate(qval=p.adjust(pval, "BH"))%>%
   ungroup()%>%
   as.data.frame()
res2 <- res2%>%mutate(motif_name=motifs[gene])

opfn <- paste(outdir, "3.2_motif_diff_fitInd_JW.rds", sep="")
write_rds(res2, opfn)



###############################################################
### compare z-score between my script and Mohammed' script  ###
###############################################################

###
### response motifs recent file
## fn <- "./sc_multiome_data/3_motif/2_motif.activities.outs/3_motif.diff.results_macs2_0.1_cn_indfilt_incsampleID_all_cols_control.rds"
fn <- "./sc_multiome_data/3_motif/2_motif.activities.outs/3_motif.diff.results_macs2_0.1_cn.rds"
res <- read_rds(fn)%>%mutate(comb=paste(MCls, contrast, sep="_"), zscore=beta/stderr)
res <- res%>%mutate(rn_id=paste(comb, gene, sep="_"))%>%
   select(rn_id, comb, MCls, contrast, zscore)

comb_set <- sort(unique(res$comb))



###
### file I generated using Mohammed's script
## fn2 <- "./4.0_motif_compare.outs/3.2_motif_diff_fitInd_JW.rds"
fn2 <- "./4.0_motif_compare.outs/3.3_motif_diff_JW.rds"
res2 <- read_rds(fn2)%>%mutate(comb=paste(MCls, contrast, sep="_"), zscore=beta/stderr)
res2 <- res2%>%filter(comb%in%comb_set)%>%
    mutate(rn_id=paste(comb, gene, sep="_"))%>%
    select(rn_id, zscore_JW=zscore)


###
### plot data
res_comb <- res%>%inner_join(res2, by="rn_id")


### colors
col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
  "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
   "6_Monocyte"="#984ea3", "7_dnT"="black")
col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
       "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3")


### plots
p <- ggplot(res_comb, aes(x=zscore_JW, y=zscore, color=contrast))+
   geom_point(size=0.2)+
   facet_wrap(~MCls, ncol=4, scales="free")+
   scale_color_manual(values=col2, guide=guide_legend(override.aes=list(size=2)))+
   geom_abline(color="grey")+ 
   xlab(bquote(italic(Z)~"-score I generated using Mohammed's script(not ind)"))+
   ylab(bquote(italic(Z)~"-score from new file"))+
   theme_bw()+
   theme(legend.title=element_blank(),
         legend.text=element_text(size=9),
         legend.key.size=unit(0.5, "cm"),
         axis.title=element_text(size=10),
         axis.text=element_text(size=10))

figfn <- paste(outdir, "Figure1.2_compare_zscore_JWvsMohammed.png", sep="")
ggsave(figfn, p, device="png", width=920, height=420, units="px", dpi=120)




#################################################################
### Compare the differential file currently used and old used ###
#################################################################

###
### response motifs recent file
## fn <- "./sc_multiome_data/3_motif/2_motif.activities.outs/3_motif.diff.results_macs2_0.1_cn_indfilt_incsampleID_all_cols_control.rds"
fn <- "./sc_multiome_data/3_motif/2_motif.activities.outs/3_motif.diff.results_macs2_0.1_cn.rds"
res <- read_rds(fn)%>%mutate(comb=paste(MCls, contrast, sep="_"), zscore=beta/stderr)
res <- res%>%mutate(rn_id=paste(comb, gene, sep="_"))%>%
   select(rn_id, comb, MCls, contrast, zscore)

comb_set <- sort(unique(res$comb))


###
### file I generated using Mohammed's script
fn2 <- "./4.0_motif_compare.outs/3.2_motif_diff_fitInd_JW.rds"
res2 <- read_rds(fn2)%>%mutate(comb=paste(MCls, contrast, sep="_"), zscore=beta/stderr)
res2 <- res2%>%filter(comb%in%comb_set)%>%
    mutate(rn_id=paste(comb, gene, sep="_"))%>%
    select(rn_id, zscore_JW=zscore)


###
### plot data
res_comb <- res%>%inner_join(res2, by="rn_id")


### colors
col1 <- c("0_CD4Naive"="#ffaa00", "1_TCM"="pink", "2_NKcell"="#aa4b56",
  "3_TEM"="blue", "4_Bcell"="#4daf4a", "5_CD8Naive"="green",
   "6_Monocyte"="#984ea3", "7_dnT"="black")
col2 <- c("caffeine"="red", "nicotine"="tan", "vitA"="tan4",
       "vitD"="seagreen4", "vitE"="salmon3", "zinc"="maroon3")


### plots
p <- ggplot(res_comb, aes(x=zscore_JW, y=zscore, color=contrast))+
   geom_point(size=0.2)+
   facet_wrap(~MCls, ncol=4, scales="free")+
   scale_color_manual(values=col2, guide=guide_legend(override.aes=list(size=2)))+
   geom_abline(color="grey")+ 
   xlab(bquote(italic(Z)~"-score I generated using Mohammed's script (fit ind)"))+
   ylab(bquote(italic(Z)~"-score from old file"))+
   theme_bw()+
   theme(legend.title=element_blank(),
         legend.text=element_text(size=9),
         legend.key.size=unit(0.5, "cm"),
         axis.title=element_text(size=10),
         axis.text=element_text(size=10))

figfn <- paste(outdir, "Figure2_compare_zscore_newJW_OldMhmd.png", sep="")
ggsave(figfn, p, device="png", width=950, height=480, units="px", dpi=120)



### monocytes
plotDF <- res_comb%>%filter(MCls=="6_Monocyte")
p2 <- ggplot(plotDF, aes(x=zscore_JW, y=zscore, color=contrast))+
   geom_point(size=0.2)+
   facet_wrap(~contrast, ncol=3, scales="free")+
   scale_color_manual(values=col2, guide=guide_legend(override.aes=list(size=2)))+
   geom_abline(color="grey")+ 
   xlab(bquote(italic(Z)~"-score I generated using Mohammed's script (fit ind)"))+
   ylab(bquote(italic(Z)~"-score from old file"))+
   theme_bw()+
   theme(legend.title=element_blank(),
         legend.text=element_text(size=9),
         legend.key.size=unit(0.5, "cm"),
         axis.title=element_text(size=10),
         axis.text=element_text(size=10))

figfn <- paste(outdir, "Figure2.2_compare_zscore_Monocyte.png", sep="")
ggsave(figfn, p2, device="png", width=820, height=480, units="px", dpi=120)



####
#### number response motifs at different threshold 
fn <- "./4.0_motif_compare.outs/3.2_motif_diff_fitInd_JW.rds"
res <- read_rds(fn2)%>%mutate(comb=paste(MCls, contrast, sep="_"), zscore=beta/stderr)

###
comb <- sort(unique(res$comb))
summ <- map_dfr(comb, function(ii){
   ###
   tmp <- res%>%dplyr::filter(comb==ii)

   ###
   ths <- quantile(abs(tmp$beta), probs=c(0.5, 0.75, 0.8, 0.9)) 
   nums <- sapply(ths, function(th0){
   ###     
      motifs <- tmp%>%filter(qval<0.1, abs(beta)>th0)%>%pull(gene)%>%unique()
      length(motifs)
   })
   df2 <- data.frame(comb=ii, t(nums))
   names(df2) <- c("comb", "motifs_th0.5", "motifs_th0.75", "motifs_th0.8", "motifs_th0.9")
   df2     
})

###
fn <- "./4.0_motif_compare.outs/3.2_motif_diff_fitInd_JW.rds"
res <- read_rds(fn2)%>%mutate(comb=paste(MCls, contrast, sep="_"), zscore=beta/stderr)
ii <- comb[37]
res2 <- res%>%filter(comb==ii)%>%arrange(desc(zscore))
res3 <- res2%>%mutate(qval2=qvalue(pval)$qvalue)
 
###
###
fn <- "./sc_multiome_data/3_motif/2_motif.activities.outs/3_motif.diff.results_macs2_0.1_cn.rds"
res <- read_rds(fn)%>%mutate(comb=paste(MCls, contrast, sep="_"), zscore=beta/stderr)
res2 <- res%>%filter(comb==ii)%>%arrange(desc(zscore))




#########################################
### Summary number of response motifs ###
#########################################


fn <- "./sc_multiome_data/3_motif/2_motif.activities.outs/3_motif.diff.results_macs2_0.1_cn_indfilt_incsampleID_all_cols_control.rds"
res <- read_rds(fn)%>%mutate(comb=paste(MCls, contrast, sep="_"))

fn2 <- "./4.0_motif_compare.outs/3.2_motif_diff_fitInd_JW.rds"
res2 <- read_rds(fn2)%>%mutate(comb=paste(MCls, contrast, sep="_"), zscore=beta/stderr)

comb <- sort(unique(res$comb))

## compare
df_olap <- map_dfr(comb, function(ii){
   ###
   cat(ii, "\n")
    
   x <- res%>%filter(comb==ii)
   th0 <- quantile(abs(x$beta), probs=0.9)
   motif <- x%>%filter(qval<0.1, abs(beta)>th0)%>%pull(gene)%>%unique()

   ##
   x2 <- res2%>%filter(comb==ii)
   th0 <- quantile(abs(x2$beta), probs=0.9)
   motif2 <- x2%>%filter(qval<0.1, abs(beta)>th0)%>%pull(gene)%>%unique()    

   summ <- data.frame(comb=ii, nmotif_Mhmd=length(motif), nmotif_JW=length(motif2),
                      nolap=length(intersect(motif, motif2)))
   summ
})

opfn <- paste(outdir, "0_compare_Mhmd_JW.xlsx", sep="")
write.xlsx(df_olap, file=opfn, overwrite=T)

### motif activities
## fn <- "./sc_multiome_data/3_motif/2_motif.activities.outs/2_motif.ave_macs2_0.1_cn_control.rds"
## mat_motif <- read_rds(fn)
## rn <- rownames(mat_motif)
## rn2 <- motifs[rn]

### motif gene expression


## fn <- "./sc_multiome_data/2_Differential/1_DiffRNA_2.outs/1_YtX.sel_clusters_0.1_cn.rds"
## mat_RNA <- read_rds(fn)
## sum(rn2%in%rownames(mat_RNA))
