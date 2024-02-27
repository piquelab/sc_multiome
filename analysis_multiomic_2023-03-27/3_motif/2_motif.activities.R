##
library(tidyverse)
library(Seurat)
library(parallel)
library(SeuratDisk)
library(SeuratData)
library(SeuratObject)
library(Signac)
library(SeuratWrappers)
library(SeuratData)
library(chromVAR)
library(JASPAR2020)
library(JASPAR2022)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38)
#library(BSgenome.Hsapiens.1000genomes.hs37d5, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
#library(BSgenome.Hsapiens.1000genomes.hs37d5)
#library(ChIPseeker, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(ChIPseeker)
###
library(ggplot2)
#library(cowplot, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(cowplot)
library(grid)
library(gridExtra)
library(ggExtra)
library(RColorBrewer)
library(ggsci)
library(viridis)
library(ComplexHeatmap)
library(circlize)
theme_set(theme_grey())


###
outdir <- "./2_motif.activities.outs/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


##################################################
#---- calculate motif----
##################################################
fn <- "./1_motif.outs/1_scATAC.motif.rds"
atac <- read_rds(fn)

## atac <- RunChromVAR(object=atac,
##    genome=BSgenome.Hsapiens.1000genomes.hs37d5)


atac <- RunChromVAR(object=atac,
   genome=BSgenome.Hsapiens.UCSC.hg38)

#opfn <- "./2_motif.activities.outs/1_scATAC.motifActivities.rds"
#write_rds(atac, opfn)

atac <- read_rds("./2_motif.activities.outs/1_scATAC.motifActivities.rds")
#atac
#head(atac@meta.data)


## #----tsne plots----
## atac <- read_rds("./2_motif.activities.outs/1_scATAC.motifActivities.rds")
## X <- atac@assays$chromvar@data

## library(tsne)
## xtsne <- tsne(t(X))
## write_rds(xtsne, "./2_motif.activities.outs/1.1_tsne.rds")








#####################################################
#----differential motif activities----
#####################################################
atac <- read_rds("./2_motif.activities.outs/1_scATAC.motifActivities.rds")
atac

#atac2 <- subset(atac, subset=MCls!="DC")
atac2 <-atac

Y <- atac2@assays$chromvar@data
Y

head(Y)
ncol(Y)

meta <- atac2@meta.data
head(meta)
nrow(meta)

## meta <- meta%>%
##    mutate(treat=gsub(".*-ATAC-|_.*", "", NEW_BARCODE),
##           bti=paste(MCls, treat, SNG.BEST.GUESS, sep="_"))%>%
##    dplyr::select(NEW_BARCODE, bti) 


meta$treat <- meta$treats

colnames(meta)

table(meta$treat)

#----cell labeling----
#----cell labeling----

meta <- meta %>% mutate(treat=ifelse(treat %in% c("water", "etOH"), "control", treat))
table(meta$treat)

meta <- meta %>% mutate(bti=paste(wsnn_res.0.1, MCls, treat, SNG.BEST.GUESS, sep="_"))%>% dplyr::select(NEW_BARCODE, bti) 

head(meta)
ncol(meta)
unique(meta$bti)

###
dd <- meta%>%
   group_by(bti)%>%
   summarise(ncell=n(), .groups="drop")%>%
   ungroup()%>%as.data.frame()
head(dd)

###
bti <- factor(meta$bti)
head(bti)
length(bti)

X <- model.matrix(~0+bti)
head(X)
ncol(X)
nrow(X)

YtX <- Y %*% X
YtX <- as.matrix(YtX)
colnames(YtX) <- gsub("^bti", "", colnames(YtX))

head(YtX)
ncol(YtX)


## Julong script
YtX_ave <- sweep(YtX, 2, dd$ncell, "/")
##YtX_ave <- sweep(YtX_sel, 2, dd$ncell, "/")

head(YtX_ave)
ncol(YtX_ave)

dd0 <- dd%>%dplyr::filter(ncell>20)
YtX_sel <- YtX_ave[,dd0$bti]
head(YtX_sel)
##write_rds(YtX_sel, "./2_motif.activities.outs/2_motif.ave_macs2_0.1_cn_control.rds")


## ## new way
## dd <- dd %>%filter(ncell>20)
## head(dd)
## ##we dont use >100 filter since there are no counts, but chromVAR activities 
## YtX_sel <- YtX[,dd$bti]
## head(YtX_sel)
## ncol(YtX_sel)
## ##YtX_ave <- sweep(YtX, 2, dd$ncell, "/")
## YtX_ave <- sweep(YtX_sel, 2, dd$ncell, "/")
## head(YtX_ave)
## ncol(YtX_ave)
## write_rds(YtX_ave, "./2_motif.activities.outs/2_motif.ave_macs2_0.1_cn_control.rds")




#----heatmap----
#mat <- read_rds("./2_motif.activities.outs/2_motif.ave_macs2_0.1_cn_control.rds")
mat <- YtX_sel
head(mat)
ncol(mat)

meta <- str_split(colnames(mat), "_", simplify=T)%>%as.data.frame()
names(meta) <- c("cluster", "MCls", "treat", "sampleID")
head(meta)
unique(meta$treat)


###
## b <- as.vector(mat)
## b
## breaks <- quantile(b, probs=seq(0, 1, length.out=100))
## col_fun <-  colorRamp2(breaks,
##    colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100))

y <- as.numeric(mat)
y
max(y)
min(y)

scale <- quantile(abs(y),probs=0.99)
scale

mybreaks <- seq(-scale, scale, length.out=100)
mybreaks
names(mybreaks) <- NULL

mycol <- colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100)



## column_ha <- HeatmapAnnotation(
##    celltype=meta$MCls,
##    treatment=meta$treat,
##    col=list(
##       celltype=c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
##                  "NKcell"="#aa4b56", "Tcell"="#ffaa00"),
##       treatment=c("CTRL"="#828282", "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
##                   "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")))

column_ha <- HeatmapAnnotation(
    celltype=paste0(meta$cluster, "_", meta$MCls),
    treatment=meta$treat,
    col=list(celltype=c("4_Bcell"="#4daf4a", "6_Monocyte"="#984ea3",
                        "2_NKcell"="#aa4b56", "0_CD4Naive"="#ffaa00",
                        "1_TCM"="pink","3_TEM"="blue",
                        "5_CD8Naive"="green", "7_dnT"="black", "8_MAIT"="grey",
                        "9_Platelet"="purple4", "10_DC"="slategray4"),
             treatment=c("caffeine"="red", "nicotine"="tan",
                         "vitA"="tan4", "vitD"="seagreen4",
                         "vitE"="salmon3", "zinc"="maroon3",
                         "control"="grey")))
                         ##"water"="grey", "etOH"="darkgrey")))

## fig <- Heatmap(mat, col=col_fun,
##    cluster_rows=T, cluster_columns=T,
##    show_row_dend=T, show_column_dend=T,
##    top_annotation=column_ha,
##    heatmap_legend_param=list(title="motif activities",
##       title_gp=gpar(fontsize=10),
##       labels_gp=gpar(fontsize=10)),
##    show_row_names=F, show_column_names=F,
##    use_raster=F, raster_device="png")

## figfn <- "./2_motif.activities.outs/Figure1.1_heatmap.motif.activities_macs2_0.1_cn.png"
## png(figfn, height=800, width=900, res=120)
## set.seed(0)
## fig <- draw(fig)
## dev.off()


fig <- Heatmap(mat,
               cluster_rows=T, cluster_columns=F,
               show_row_dend=T, show_column_dend=F,
               top_annotation=column_ha,
               heatmap_legend_param=list(title="motif activities",
                                         title_gp=gpar(fontsize=18),
                                         labels_gp=gpar(fontsize=18)),
               show_row_names=F, show_column_names=F,
               use_raster=F, raster_device="png")

figfn <- "./2_motif.activities.outs/Figure1.1_heatmap.motif.activities_macs2_0.1_cn_2_nocluster_control.png"
png(figfn, height=1800, width=2100, res=250)
#set.seed(0)
fig <- draw(fig)
dev.off()








###############################
### Differential activities ###
###############################
mat <- read_rds("./2_motif.activities.outs/2_motif.ave_macs2_0.1_cn_control.rds")
head(mat)
ncol(mat)

bti <-colnames(mat)
head(bti)

cvt0 <- str_split(bti, "_", simplify=T)%>%as.data.frame()
head(cvt0)

cvt <- data.frame(bti=bti, MCls=paste0(cvt0[,1], "_", cvt0[,2]), treat=cvt0[,3], sampleID=cvt0[,4])
head(cvt)
#cvt
#unique(cvt$treat)
table(cvt$MCls, cvt$treat)

### Differential procedure
#MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
MCls <- unique(cvt$MCls)
MCls





# filtered
#-------------------------------------------------------------------------
cvt <- cvt%>%mutate(MCls_IDs=paste(cvt$MCls, ".", cvt$sampleID, sep=""))
head(cvt)

MCls.IDs <- unique(cvt$MCls_IDs)
#MCls.IDs

cvt.filtered <- map_dfr(MCls.IDs, function(oneX){
          time0 <- Sys.time()
                  ###
                  cvt.2 <- cvt%>%dplyr::filter(MCls==str_split(oneX, "\\.", simplify=T)[1], sampleID==str_split(oneX, "\\.", simplify=T)[2])
                  ##
                  cvt.2 <- cvt.2%>%mutate(count=ifelse(nrow(cvt.2)>=7,1,0))
                  ##
                  time1 <- Sys.time()
                  elapsed <- difftime(time1, time0, units="mins")
                  cat(oneX, elapsed, "Done\n")
                  cvt.2
                })
#cvt.filtered

cvt.2 <- cvt.filtered %>% dplyr::filter(count==1)
head(cvt.2)
#unique(cvt.2$treat)
#table(cvt.2$MCls, cvt$treat)

comb.2 <- table(cvt.2$MCls, cvt.2$treat)
comb.2
#write.csv(comb.2, paste0("./1_DiffRNA_2.outs/combinations_2_", reso, "_cn.csv"), quote = FALSE, row.names = TRUE)
#comb.2 <- read.csv("./1_DiffRNA_2.outs/combinations_2.csv")
#comb.2

## bti3 <- cvt.2$bti
## #bti3

## #ncol(YtX)
## YtX <- YtX[,bti3]
## #ncol(YtX)

## dd2 <- dd %>% dplyr::filter(bti %in% bti3)
## #dd2
## #sum(dd2$ncell)


MCls <- unique(cvt.2$MCls)
MCls


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
    TMP <- mclapply(rn, function(ii){
        y <- mati[ii,]
        dd <- myDE(y, X, Z, ii) #adding Z in this
        dd
    },mc.cores=1)
    TMP <- TMP[!is.na(TMP)]
    TMP <- do.call(rbind, TMP)%>%as.data.frame()%>%mutate(MCls=oneMCl)
    TMP 
})    

head(res)



### add qvalue
res2 <- res%>%group_by(MCls, contrast)%>%
   mutate(qval=p.adjust(pval, "BH"))%>%
   ungroup()%>%
   as.data.frame()
head(res2)

atac <- read_rds("./2_motif.activities.outs/1_scATAC.motifActivities.rds")
motif <- Motifs(atac)

motif

x <- unlist(motif@motif.names)
head(x)

res2 <- res2%>%mutate(motif=x[gene])
head(res2)

###
#opfn <- "./2_motif.activities.outs/3_motif.diff.results_macs2_0.1_cn.rds"
opfn <- "./2_motif.activities.outs/3_motif.diff.results_macs2_0.1_cn_indfilt_incsampleID_all_cols_control.rds"
write_rds(res2, opfn)


#res <- read_rds("./2_motif.activities.outs/3_motif.diff.results_macs2_0.1_cn.rds")
#res <- read_rds("./2_motif.activities.outs/3_motif.diff.results_macs2_0.1_cn_all_cols.rds")
#res <- read_rds("./2_motif.activities.outs/3_motif.diff.results_macs2_0.1_cn_indfilt_all_cols.rds")

#head(res)

res <- res2

end









## #----individual step----
MCls = "6_Monocyte"

cvti <- cvt%>%dplyr::filter(MCls=="6_Monocyte")
head(cvti)
unique(cvti$treat)

mati <- mat[,cvti$bti]
head(mati)

X <- data.frame(x1=cvti$treat)
head(X)

Z <- data.frame(x1=cvti$sampleID)
head(Z)

rn <- rownames(mati)
head(rn)

y <- mati["MA0030.1",]
head(y)
length(y)

#dd <- myDE(y, X, "MA0030.1")
#dd


#myDE <- function(y, X, gene){
    ###


con.ls <- list("caffeine"=c("control", "caffeine"),
                   "nicotine"=c("control", "nicotine"),
                   "vitA"=c("control", "vitA"),
                   "vitD"=c("control", "vitD"),
                   "vitE"=c("control", "vitE"),
                   "zinc"=c("control", "zinc"))
    
con.ls <- list("caffeine"=c("water", "caffeine"),
                  "nicotine"=c("water", "nicotine"),
                  "vitA"=c("etOH", "vitA"),
                  "vitD"=c("etOH", "vitD"),
                  "vitE"=c("etOH", "vitE"),
                  "water"=c("etOH", "water"),
                  "zinc"=c("water", "zinc"))

   ## con.ls <- list("LPS"=c("CTRL","LPS"),
   ##    "LPS-DEX"=c("LPS","LPS-DEX"),
   ##    "PHA"=c("CTRL", "PHA"),
   ##    "PHA-DEX"=c("PHA","PHA-DEX"))
   ## Contrast <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")


Contrast <- c("caffeine", "nicotine", "vitA", "vitD",
                 "vitE", "zinc")

Contrast <- c("caffeine", "nicotine", "vitA", "vitD",
                 "vitE", "water", "zinc")


###
## linear regression
x1 <- X[,1]
x1

z1 <- Z[,1]
z1

lm0 <- try(lm(y~0+x1+z1), silent=T)
lm0

### 
#if ( (class(lm0)!="try-error")){
    
b <- coef(lm0)
b

nn <- gsub("x[12]", "", names(b))
nn

names(b) <- nn
b

vb <- diag(vcov(lm0))
vb

names(vb) <- nn
vb

###
## Contrast

#dd <- lapply(Contrast,function(one){

con0 <- con.ls[["caffeine"]]
con0


#if (

all(con0%in%nn)

#){

bhat <- b[con0[2]]-b[con0[1]]
bhat

sdhat <- sqrt(vb[con0[2]]+vb[con0[1]])
sdhat

z <- bhat/sdhat
z

#pnorm #this is a function
p <- 2*pnorm(-abs(z))
p

dd1 <- data.frame(gene="MA0030.1", beta=bhat, stderr=sdhat, pval=p, contrast=con0[1])
dd1


#}else{
#dd1 <- NA                                                                       }
#dd1
#})

#dd <- dd[!is.na(dd)]
#dd <- do.call(rbind,dd)

###
#   }else{
#      dd <- NA
#   } ## End if try-error
#   dd 
#} ###


TMP <- mclapply(rn, function(ii){
    y <- mati[ii,]
    dd <- myDE(y, X, ii)
     dd
   },mc.cores=1)

TMP <- TMP[!is.na(TMP)]
head()

TMP <- do.call(rbind, TMP)%>%as.data.frame()%>%mutate(MCls=oneMCl)
head()

TMP 
head()
## #----individual step----

## myDE(y, X, "MA0030.1")



cvti <- cvt.2%>%dplyr::filter(MCls=="6_Monocyte")

head(cvti)

X <- data.frame(x1=cvti$treat)
X

Z <- data.frame(x1=cvti$sampleID)
Z

x1 <- X[,1]
x1

z1 <- Z[,1]
z1




    





## #----check for monocytes----
head(res)


## table(res$MCls, res$contrast)

unique(res$MCls)

res.mo.water <- res %>% dplyr::filter(contrast == "water", MCls=="6_Monocyte")
head(res.mo.water)
nrow(res.mo.water)


res.mo.caff <- res %>% dplyr::filter(contrast == "caffeine", MCls=="6_Monocyte")
head(res.mo.caff)
nrow(res.mo.caff)


res.mo.vitA <- res %>% dplyr::filter(contrast == "vitA", MCls=="6_Monocyte")
head(res.mo.vitA)
nrow(res.mo.vitA)


res.cd4.water <- res %>% dplyr::filter(contrast == "water", MCls=="0_CD4Naive")
head(res.cd4.water)
nrow(res.cd4.water)


## max(res.mo.water$beta)

## sig.res.mo.water <- res.mo.water%>% dplyr::filter(qval<0.1,abs(beta)>0)

## nrow(sig.res.mo.water)

## write.csv(res.mo.water,"./2_motif.activities.outs/res_mo_water.csv", row.names = FALSE)

end




## #----check for monocytes function loop----
## mat <- read_rds("./2_motif.activities.outs/2_motif.ave_macs2_0.1_cn.rds")
## head(mat)

## bti <-colnames(mat)
## head(bti)

## cvt0 <- str_split(bti, "_", simplify=T)%>%as.data.frame()
## head(cvt0)

## cvt <- data.frame(bti=bti, MCls=paste0(cvt0[,1], "_", cvt0[,2]), treat=cvt0[,3], sampleID=cvt0[,4])

## head(cvt)
## #cvt
## #unique(cvt$treat)
## #table(cvt$MCls, cvt$treat)

## MCls <- unique(cvt$MCls)
## MCls


## rn <- rownames(mat)
## rn

## dat <- map_dfr(MCls, function(oneMCl){
##        ###
##        cvti <- cvt%>%dplyr::filter(MCls==oneMCl)
##        mati <- mat[,cvti$bti]
##        X <- data.frame(x1=cvti$treat)
##        rn <- rownames(mati)
##        TMP <- mclapply(rn, function(ii){
##                 y <- mati[ii,]
##                      dd <- myDE(y, X, ii)
##                      dd
##                    },mc.cores=1)
##        TMP <- TMP[!is.na(TMP)]
##        TMP <- do.call(rbind, TMP)%>%as.data.frame()%>%mutate(MCls=oneMCl)
##        TMP
##     })    

## dat <- data.frame()

## ##dat <- rbind.data.frame(dat,

## for(i in MCls){
##     tryCatch({
##     cvti <- cvt%>%dplyr::filter(MCls==i)
##     #head(cvti)
##     #unique(cvti$treat)
##     mati <- mat[,cvti$bti]
##     #head(mati)
##     X <- data.frame(x1=cvti$treat)
##     #head(X)
##     rn <- rownames(mati)
##     #head(rn)
##     for (j in rn){
##         y <- mati[j,]
##         #head(y)
##         #length(y)
##         #dd <- myDE(y, X, "MA0030.1")
##         #dd
##         #myDE <- function(y, X, gene){
##         ###
##         con.ls <- list("caffeine"=c("water", "caffeine"),
##                        "nicotine"=c("water", "nicotine"),
##                        "vitA"=c("etOH", "vitA"),
##                        "vitD"=c("etOH", "vitD"),
##                        "vitE"=c("etOH", "vitE"),
##                        "water"=c("etOH", "water"),
##                        "zinc"=c("water", "zinc"))
##         Contrast <- c("caffeine", "nicotine", "vitA", "vitD",
##                       "vitE", "water", "zinc")
##         ## linear regression
##         x1 <- X[,1]
##         #x1
##         lm0 <- try(lm(y~0+x1), silent=T)
##         #lm0
##         ### 
##         #if ( (class(lm0)!="try-error")){
##         b <- coef(lm0)
##         #b
##         nn <- gsub("x[12]", "", names(b))
##         #nn
##         names(b) <- nn
##         #b
##         vb <- diag(vcov(lm0))
##         #vb
##         names(vb) <- nn
##         #vb
##         ###
##         ## Contrast
##         #dd <- lapply(Contrast,function(one){
##         for(k in Contrast){ 
##             con0 <- con.ls[[k]]
##             #con0
##             ##if (
##             ##all(con0%in%nn)
##             ##){
##             bhat <- b[con0[2]]-b[con0[1]]
##             #bhat
##             sdhat <- sqrt(vb[con0[2]]+vb[con0[1]])
##             #sdhat
##             z <- bhat/sdhat
##             #z
##             ##pnorm #this is a function
##             p <- 2*pnorm(-abs(z))
##             #p
##             dd1 <- data.frame(gene="MA0030.1", beta=bhat, stderr=sdhat, pval=p, contrast="water")
##             ##dd1
##             return(c(j, i, k, b[con0[2]], b[con0[1]], bhat, sdhat, z, p))
##                   ##}else{
##                   ##dd1 <- NA                                                                       ##}
##                   ##dd1
##                   ##})
##         }
##     }
##     },   error=function(e){cat("ERROR", "\n")})
## })

## head(dat)

## #dd <- dd[!is.na(dd)]
## #dd <- do.call(rbind,dd)
## ###
## #   }else{
## #      dd <- NA
## #   } ## End if try-error
## #   dd 
## #} ###
## TMP <- mclapply(rn, function(ii){
##     y <- mati[ii,]
##     dd <- myDE(y, X, ii)
##      dd
##    },mc.cores=1)

## TMP <- TMP[!is.na(TMP)]
## head()

## TMP <- do.call(rbind, TMP)%>%as.data.frame()%>%mutate(MCls=oneMCl)
## head()

## TMP 
## head()
















#----significant res----
###
## fn <- "./1.2_motif.outs/3_motif.enrich.direction.rds"
## enrich <- read_rds(fn)
## drt2 <- c("0"="Down", "1"="Up")
## enrich <- enrich%>%mutate(direction2=drt2[as.character(direction)],
##    cluster=paste(contrast, MCls, direction2, sep="."),
##    newCluster=clst2[cluster])
###
## topmotif <- enrich%>%
##    group_by(cluster)%>%
##    top_n(n=6, wt=fold.enrichment)%>%ungroup()%>%
##     dplyr::pull(motif)%>%unique()

## topmotif2 <- enrich%>%
##     filter(motif%in%topmotif, fold.enrichment>1.41, qvalue.fisher<0.1)%>%
##     dplyr::pull(motif.name)%>%unique()

sigs <- res%>%dplyr::filter(qval<0.1,abs(beta)>1.41)%>%
   group_by(MCls, contrast)%>%
   summarise(ny=n(), .groups="drop")
head(sigs)

# log Fc of 0.25 use and 0
# https://www.omnicalculator.com/math/log-2
# 0.5 <- 1.41
# 0.25 <- 1.19
# 0 <- 1
topmotif <- res%>%
   dplyr::filter(qval<0.1, abs(beta)>0)%>%
   dplyr::pull(motif)
#head(topmotif)
length(topmotif)

topmotif <- sort(unique(topmotif))
#head(topmotif)
length(topmotif)





#------------------------------------------------
#----data for heatmap and data with topmotifs----
#------------------------------------------------
res <- res%>%mutate(condition=paste(MCls, contrast, sep="_"))
head(res)
nrow(res)
#table(res$condition)
table(res$MCls, res$contrast)
length(unique(res$motif))


res.filt <-  res%>%dplyr::filter(qval<0.1, abs(beta)>0) 
head(res.filt)
nrow(res.filt)

sig_motifs_table <- table(res.filt$MCls, res.filt$contrast)
sig_motifs_table

length(unique(res.filt$motif))

library(reshape2)
sig.matrix <-as.matrix(sig_motifs_table) # create a numeric matrix object
sig.melt <- melt(sig.matrix)

sig.p <- ggplot(sig.melt, aes(x = Var2, y = Var1)) +
        geom_tile(aes(fill=value), color="black")+
        scale_fill_gradient(low = "white", high = "white")+
        geom_text(aes(label = value), color = "black", size = 5) +
        ggtitle(paste0("Significant_motifs_coef_0")) +
        theme(axis.text.x = element_text(angle = 90),
                        axis.title.x = element_blank(),
                        axis.title.y = element_blank(),
                        legend.position="none")



library("ggpubr")
figure <- ggarrange(sig.p + font("x.text", size = 14))



png(paste0("./2_motif.activities.outs/sig_motifs_table_coef_0.png"), width=500, height=500, pointsize=16, res=125)
print(figure)
dev.off()



## res.topmotif <- res%>%dplyr::filter(motif %in% topmotif)
## head(res.topmotif)
## nrow(res.topmotif)
## table(res.topmotif$MCls, res.topmotif$contrast)
## length(unique(res.topmotif$motif))


## condition <- sort(unique(res$condition))
## condition




#---------------------------
#----mat.5topmotifs.qval----
#---------------------------
head(res)

mat.5topmotifs.qval <- map_dfr(condition, function(ii){
  # res2 <- res.topmotif %>% dplyr::filter(condition==ii)
   res2 <- res.filt %>% dplyr::filter(condition==ii)
   res2 <- res2 %>% arrange(qval)
   res2 <- res2 [1:5, ]
   #b <- res2$beta
   #names(b) <- res2$motif
   #b[topmotif]   
})
head(mat.5topmotifs.qval)
nrow(mat.5topmotifs.qval)
table(mat.5topmotifs.qval$MCls, mat.5topmotifs.qval$contrast)
#length(unique(mat.5topmotifs.qval$motif))

topmotif.qval <- unique(mat.5topmotifs.qval$motif)
head(topmotif.qval)
length(topmotif.qval)


mat.5topmotifs.qval.2 <- map_dfc(condition, function(ii){
    res2 <- mat.5topmotifs.qval %>% dplyr::filter(condition==ii)
    #res2 <- res2 %>% arrange(desc(qval))
    #res2 <- res2[1:5, ]
    #res2
    b <- res2$beta
    names(b) <- res2$motif
    b[topmotif.qval]
    })
head(mat.5topmotifs.qval.2)
nrow(mat.5topmotifs.qval.2)


mat <- as.matrix(mat.5topmotifs.qval.2)
head(mat)

mat[is.na(mat)] = 0
head(mat)



#----mat all values but selecting only top5 motifs, NA not equal to zero----
condition

mat <- map_dfc(condition, function(ii){
   res2 <- res%>%dplyr::filter(condition==ii)
   b <- res2$beta
   names(b) <- res2$motif
   b[topmotif.qval]
})
head(mat)

mat <- as.matrix(mat)
head(mat)
nrow(mat)

any(is.na(mat))

is.na(mat) %>% table()


mat[is.na(mat)] = 0 #either this or step below
head(mat)
#----mat all values but selecting only top5 motifs, NA not equal to zero----


#----RESUME----
colnames(mat) <- condition
rownames(mat) <- topmotif.qval
head(mat)
nrow(mat)

#----OPTIONAL----
mat<- na.omit(mat)


x <- as.data.frame(str_split(condition, "_", simplify=T))
x <- x %>% dplyr::arrange(x$V3) %>% mutate(cond=paste0(V1, "_", V2, "_", V3))
head(x)

condition.2 <- x$cond
condition.2

mat3 <- mat[, condition.2]
head(mat3)


#----heatmap----
b <- as.vector(mat)

breaks <- quantile(b, probs=seq(0, 1, length.out=100),na.rm=T)

col_fun <-  colorRamp2(breaks,
   colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100))

condition.2

celltype <- str_split(condition.2, "_", simplify=T)
celltype
celltype.2 <- paste(celltype[,1], "_", celltype[,2], sep="")
celltype.2
treatment.2 <- celltype[,3]
treatment.2


celltype <- str_split(condition, "_", simplify=T)
celltype
celltype.2 <- paste(celltype[,1], "_", celltype[,2], sep="")
celltype.2
treatment.2 <- celltype[,3]
treatment.2


column_ha <- HeatmapAnnotation(
        celltype=celltype.2,
        treatment=treatment.2,
        col=list(celltype=c("4_Bcell"="#4daf4a", "6_Monocyte"="#984ea3",
                            "2_NKcell"="#aa4b56", "0_CD4Naive"="#ffaa00",
                            "1_TCM"="pink", "3_TEM"="blue",
                            "5_CD8Naive"="green", "7_dnT"="black",
                            "8_MAIT"="grey"),
                 treatment=c("caffeine"="red", "nicotine"="tan",
                             "vitA"="tan4", "vitD"="seagreen4",
                             "vitE"="salmon3", "water"="grey", "zinc"="maroon3"))
)

fig <- Heatmap(mat,
   cluster_rows=T,
   cluster_columns=F,
   show_row_dend=T, show_column_dend=T,
   top_annotation=column_ha,
   heatmap_legend_param=list(title="Diff motif",
      title_gp=gpar(fontsize=10),
      labels_gp=gpar(fontsize=10)),
   show_row_names=T, show_column_names=T,
   row_names_side="left",
   column_names_gp=gpar(fontsize=9),
   row_names_gp=gpar(fontsize=6),
   use_raster=F, raster_device="png")

figfn <- "./2_motif.activities.outs/Figure1.2_heatmap.motif.activities_macs2_0.1_cn_noclust_celltypes_FC_0_top5_qval_filt_data_allvalues.png"
png(figfn, height=3200, width=2800, res=225)
#set.seed(0)
fig <- draw(fig)
dev.off()




#---------------------------
#----mat.5topmotifs.beta----
#---------------------------
mat.5topmotifs.beta <- map_dfr(condition, function(ii){
    res2 <- res %>% dplyr::filter(condition==ii)
    res2 <- res2 %>% arrange(desc(beta))
    res2 <- res2[1:5, ]
    res2
    #b <- res2$beta
    #names(b) <- res2$motif
    #b[topmotif]
    })
head(mat.5topmotifs.beta)
nrow(mat.5topmotifs.beta)
table(mat.5topmotifs.beta$motif)


topmotif.beta <- sort(unique(mat.5topmotifs.beta$motif))
topmotif.beta


mat.5topmotifs.beta.2 <- map_dfc(condition, function(ii){
    res2 <- mat.5topmotifs.beta %>% dplyr::filter(condition==ii)
    #res2 <- res2 %>% arrange(desc(beta))
    #res2 <- res2[1:5, ]
    #res2
    b <- res2$beta
    names(b) <- res2$motif
    b[topmotif.beta]
    })
head(mat.5topmotifs.beta.2)
nrow(mat.5topmotifs.beta.2)


mat <- as.matrix(mat.5topmotifs.beta.2)
head(mat)

mat[is.na(mat)] = 0
head(mat)

colnames(mat) <- condition
rownames(mat) <- topmotif.beta
head(mat)


x <- as.data.frame(str_split(condition, "_", simplify=T))
x <- x %>% dplyr::arrange(x$V3) %>% mutate(cond=paste0(V1, "_", V2, "_", V3))
head(x)

condition.2 <- x$cond
condition.2

mat3 <- mat[, condition.2]
head(mat3)


#----heatmap----
b <- as.vector(mat)

breaks <- quantile(b, probs=seq(0, 1, length.out=100),na.rm=T)

col_fun <-  colorRamp2(breaks,
   colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100))

condition.2

celltype <- str_split(condition.2, "_", simplify=T)
celltype

celltype.2 <- paste(celltype[,1], "_", celltype[,2], sep="")
celltype.2

treatment.2 <- celltype[,3]
treatment.2

column_ha <- HeatmapAnnotation(
        celltype=celltype.2,
        treatment=treatment.2,
        col=list(celltype=c("4_Bcell"="#4daf4a", "6_Monocyte"="#984ea3",
                            "2_NKcell"="#aa4b56", "0_CD4Naive"="#ffaa00",
                            "1_TCM"="pink", "3_TEM"="blue",
                            "5_CD8Naive"="green", "7_dnT"="black",
                            "8_MAIT"="grey"),
                 treatment=c("caffeine"="red", "nicotine"="tan",
                             "vitA"="tan4", "vitD"="seagreen4",
                             "vitE"="salmon3", "water"="grey", "zinc"="maroon3"))
)

fig <- Heatmap(mat3,
   cluster_rows=T,
   cluster_columns=F,
   show_row_dend=T, show_column_dend=T,
   top_annotation=column_ha,
   heatmap_legend_param=list(title="Diff motif",
      title_gp=gpar(fontsize=10),
      labels_gp=gpar(fontsize=10)),
   show_row_names=T, show_column_names=T,
   row_names_side="left",
   column_names_gp=gpar(fontsize=9),
   row_names_gp=gpar(fontsize=6),
   use_raster=F, raster_device="png")

figfn <- "./2_motif.activities.outs/Figure1.2_heatmap.motif.activities_macs2_0.1_cn_noclust_treat_FC_0_top5_beta.png"
png(figfn, height=3200, width=2800, res=225)
#set.seed(0)
fig <- draw(fig)
dev.off()




#---------------------------
#----topmotif----
#---------------------------
mat <- map_dfc(condition, function(ii){
   res2 <- res%>%dplyr::filter(condition==ii)
   b <- res2$beta
   names(b) <- res2$motif
   b[topmotif]
})
head(mat)

mat <- as.matrix(mat)
head(mat)

colnames(mat) <- condition
rownames(mat) <- topmotif
head(mat)


x <- as.data.frame(str_split(condition, "_", simplify=T))
x <- x %>% dplyr::arrange(x$V3) %>% mutate(cond=paste0(V1, "_", V2, "_", V3))
head(x)

condition.2 <- x$cond
condition.2

mat3 <- mat[, condition.2]
head(mat3)


#----heatmap----
###
###
## condition2 <- c("Bcell_"

b <- as.vector(mat)

breaks <- quantile(b, probs=seq(0, 1, length.out=100),na.rm=T)

col_fun <-  colorRamp2(breaks,
   colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100))

## column_ha <- HeatmapAnnotation(
##    celltype=gsub("_.*", "", condition),
##    treatment=gsub(".*_", "", condition),
##    col=list(
##       celltype=c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
##                  "NKcell"="#aa4b56", "Tcell"="#ffaa00"),
##       treatment=c("LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
##                   "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")))

## fig <- Heatmap(mat, col=col_fun,
##    cluster_rows=T, cluster_columns=T,
##    show_row_dend=F, show_column_dend=T,
##    top_annotation=column_ha,
##    heatmap_legend_param=list(title="Diff motif",
##       title_gp=gpar(fontsize=10),
##       labels_gp=gpar(fontsize=10)),
##    show_row_names=T, show_column_names=T,
##    row_names_side="left",
##    column_names_gp=gpar(fontsize=9),
##    row_names_gp=gpar(fontsize=6),
##    use_raster=F, raster_device="png")


condition

celltype <- str_split(condition.2, "_", simplify=T)
celltype

celltype.2 <- paste(celltype[,1], "_", celltype[,2], sep="")
celltype.2

treatment.2 <- celltype[,3]
treatment.2

column_ha <- HeatmapAnnotation(
        celltype=celltype.2,
        treatment=treatment.2,
        col=list(celltype=c("4_Bcell"="#4daf4a", "6_Monocyte"="#984ea3",
                            "2_NKcell"="#aa4b56", "0_CD4Naive"="#ffaa00",
                            "1_TCM"="pink", "3_TEM"="blue",
                            "5_CD8Naive"="green", "7_dnT"="black",
                            "8_MAIT"="grey"),
                 treatment=c("caffeine"="red", "nicotine"="tan",
                             "vitA"="tan4", "vitD"="seagreen4",
                             "vitE"="salmon3", "water"="grey", "zinc"="maroon3"))
)

fig <- Heatmap(mat3,
   cluster_rows=T,
   cluster_columns=F,
   show_row_dend=T, show_column_dend=T,
   top_annotation=column_ha,
   heatmap_legend_param=list(title="Diff motif",
      title_gp=gpar(fontsize=10),
      labels_gp=gpar(fontsize=10)),
   show_row_names=T, show_column_names=T,
   row_names_side="left",
   column_names_gp=gpar(fontsize=9),
   row_names_gp=gpar(fontsize=6),
   use_raster=F, raster_device="png")

figfn <- "./2_motif.activities.outs/Figure1.2_heatmap.motif.activities_macs2_0.1_cn_noclust_treat_FC_0.png"
png(figfn, height=3200, width=2800, res=225)
#set.seed(0)
fig <- draw(fig)
dev.off()



### swap
## row_ha <- HeatmapAnnotation(
##    celltype=gsub("_.*", "", condition),
##    treatment=gsub(".*_", "", condition),
##    col=list(
##       celltype=c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
##                  "NKcell"="#aa4b56", "Tcell"="#ffaa00"),
##       treatment=c("LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
##                   "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")))

## fig <- Heatmap(t(mat), col=col_fun,
##    cluster_rows=T, cluster_columns=T,
##    show_row_dend=T, show_column_dend=F,
##    top_annotation=row_ha,
##    heatmap_legend_param=list(title="Diff motif",
##       title_gp=gpar(fontsize=10),
##       labels_gp=gpar(fontsize=10), legend_direction="horizontal"),
##    show_row_names=T, show_column_names=T,
##    row_names_side="left",
##    column_names_gp=gpar(fontsize=6),
##    row_names_gp=gpar(fontsize=9),
##    use_raster=F, raster_device="png")

## figfn <- "./2_motif.activities.outs/Figure1.3_heatmap.motif.activities.png"
## png(figfn, height=550, width=700, res=120)
## set.seed(0)
## fig <- draw(fig)
## dev.off()












#################################
### write differential motif ###
################################
res <- read_rds("./2_motif.activities.outs/3_motif.diff.results_macs2_0.1_cn.rds")
##
## motif.name <- sort(unique(res$motif))
## motif.id <- ConvertMotifID(object=motif, name=motif.name)
## motif.annot <- data.frame(motif.name=motif.name, motif.id=motif.id)

## res <- res%>%left_join(motif.annot, by=c("motif"="motif.name"))

sigs <- res%>%mutate(condition=paste(MCls, contrast, sep="_"))%>%
   dplyr::filter(qval<0.1, beta>0)%>%  ### abs(beta)>0.32)%>% ##90%
   group_by(condition)%>%
   summarise(ny=n(), .groups="drop")
head(sigs)


###
###
res2 <- res%>%mutate(condition=paste(MCls, contrast, sep="_"))%>%
   dplyr::filter(qval<0.1, beta>0)##abs(beta)>0.32)
head(res2)

condition <- sort(unique(res2$condition))
motif_response <- map(condition, function(one){
    res2%>%dplyr::filter(condition==one)%>%dplyr::pull(gene)
})
names(motif_response) <- condition
head(motif_response)

opfn <- paste(outdir, "4.2_response.motif.positive_macs2_0.1_cn.rds", sep="")
write_rds(motif_response, opfn)


##################################
### define  response motif (1) ###
##################################

## rm(list=ls())

## outdir <- "./2_motif.activities.outs/"

## res <- read_rds("./2_motif.activities.outs/3_motif.diff.results.rds")

## ##
## ii <- 4
## oneMCl <- "Tcell"
## res2 <- res%>%filter(MCls==oneMCl)

## th0 <- quantile(abs(res2$beta),probs=0.9)
## ## th0 <- 0.21
## ### CTRL
## res3 <- res2%>%filter(contrast=="LPS"|contrast=="PHA", qval<0.1, abs(beta)>th0)%>%
##     arrange(desc(abs(beta)))
## opfn <- paste(outdir, "5_", ii, "_",  oneMCl, "_CTRL.motif.txt", sep="")
## write.table(unique(res3$gene), opfn, row.names=F, quote=F, col.names=F)
## length(unique(res3$gene))

## ### LPS-EtOH
## res3 <- res2%>%filter(contrast=="LPS", qval<0.1, abs(beta)>th0)%>%
##     arrange(desc(abs(beta)))
## opfn <- paste(outdir, "5_",  ii, "_", oneMCl, "_LPS-EtOH.motif.txt", sep="")
## write.table(unique(res3$gene), opfn, row.names=F, quote=F, col.names=F)
## length(unique(res3$gene))


## ### LPS-DEX
## res3 <- res2%>%filter(contrast=="LPS-DEX", qval<0.1, abs(beta)>th0)%>%
##     arrange(desc(abs(beta)))
## opfn <- paste(outdir, "5_",  ii, "_", oneMCl, "_LPS-DEX.motif.txt", sep="")
## write.table(unique(res3$gene), opfn, row.names=F, quote=F, col.names=F)
## length(unique(res3$gene)) 


## ### PHA-EtOH
## res3 <- res2%>%filter(contrast=="PHA", qval<0.1, abs(beta)>th0)%>%
##     arrange(desc(abs(beta)))
## opfn <- paste(outdir, "5_", ii, "_", oneMCl, "_PHA-EtOH.motif.txt", sep="")
## write.table(unique(res3$gene), opfn, row.names=F, quote=F, col.names=F)
## length(unique(res3$gene))


## ### PHA-DEX
## res3 <- res2%>%filter(contrast=="PHA-DEX", qval<0.1, abs(beta)>th0)%>%
##     arrange(desc(abs(beta)))
## opfn <- paste(outdir, "5_", ii, "_", oneMCl, "_PHA-DEX.motif.txt", sep="")
## write.table(unique(res3$gene), opfn, row.names=F, quote=F, col.names=F)
## length(unique(res3$gene))




##################################
### define  response motif (2) ###
##################################

### consider direction, positive and negative

## rm(list=ls())


## outdir <- "./2_motif.activities.outs/"

## res <- read_rds("./2_motif.activities.outs/3_motif.diff.results.rds")

## ##
## ii <- 4
## oneMCl <- "Tcell"
## res2 <- res%>%filter(MCls==oneMCl)

## th0 <- quantile(abs(res2$beta),probs=0.9)

## ### CTRL
## res3 <- res2%>%filter(contrast=="LPS"|contrast=="PHA", qval<0.1, abs(beta)>th0, beta<0)
## ###
## opfn <- paste(outdir, "6_", ii, "_",  oneMCl, "_CTRL.motif.txt", sep="")
## write.table(unique(res3$gene), opfn, row.names=F, quote=F, col.names=F)
## length(unique(res3$gene))

## ### LPS-EtOH
## res3.1 <- res2%>%filter(contrast=="LPS", qval<0.1, abs(beta)>th0, beta>0)
## ##
## res3.2 <- res2%>%filter(contrast=="LPS-DEX", qval<0.1, abs(beta)>th0, beta<0)
## res3 <- rbind(res3.1, res3.2)
## opfn <- paste(outdir, "6_",  ii, "_", oneMCl, "_LPS-EtOH.motif.txt", sep="")
## write.table(unique(res3$gene), opfn, row.names=F, quote=F, col.names=F)
## length(unique(res3$gene))


## ### LPS-DEX
## res3 <- res2%>%filter(contrast=="LPS-DEX", qval<0.1, abs(beta)>th0, beta>0)
## ##
## opfn <- paste(outdir, "6_",  ii, "_", oneMCl, "_LPS-DEX.motif.txt", sep="")
## write.table(unique(res3$gene), opfn, row.names=F, quote=F, col.names=F)
## length(unique(res3$gene)) 


## ### PHA-EtOH
## res3.1 <- res2%>%filter(contrast=="PHA", qval<0.1, abs(beta)>th0, beta>0)
## res3.2 <- res2%>%filter(contrast=="PHA-DEX", qval<0.1, abs(beta)>th0, beta<0) 
## res3 <- rbind(res3.1, res3.2)
## ##
## opfn <- paste(outdir, "6_", ii, "_", oneMCl, "_PHA-EtOH.motif.txt", sep="")
## write.table(unique(res3$gene), opfn, row.names=F, quote=F, col.names=F)
## length(unique(res3$gene))


## ### PHA-DEX
## res3 <- res2%>%filter(contrast=="PHA-DEX", qval<0.1, abs(beta)>th0, beta>0)
## opfn <- paste(outdir, "6_", ii, "_", oneMCl, "_PHA-DEX.motif.txt", sep="")
## write.table(unique(res3$gene), opfn, row.names=F, quote=F, col.names=F)
## length(unique(res3$gene))



#################################
### define response motif (3) ###
#################################

### union of treatment motifs

res <- read_rds("./2_motif.activities.outs/3_motif.diff.results_macs2_0.1_cn.rds")

##
ii <- 4

oneMCl <- "1_TCM"

res2 <- res%>%filter(MCls==oneMCl)
head(res2)

th0 <- quantile(abs(res2$beta),probs=0.9)
th0

res3 <- res2%>%filter(qval<0.1, abs(beta)>th0)
head(res3)


##
opfn <- paste(outdir, "7_", ii, "_", oneMCl, "_union.motif_macs2_0.1_cn.txt", sep="")
write.table(unique(res3$gene), opfn, row.names=F, quote=F, col.names=F)
length(unique(res3$gene))

outdir <- "2_motif.activities.outs/"
#fn <- paste(outdir, "7_4_Tcell_union.motif.txt", sep="")
fn <- paste(outdir, "7_4_1_TCM_union.motif_macs2_0.1_cn.txt", sep="")
x <- read.table(fn)
x
##############################################################
### identify motif with high activities in some celll type ###
##############################################################
## atac2 <- subset(atac2, subset=MCls!="DC")

## Y <- atac2@assays$chromvar@data
## meta <- atac2@meta.data%>%
##    mutate(treat=gsub(".*-ATAC-|_.*", "", NEW_BARCODE),
##           bti=paste(MCls, treat, SNG.BEST.GUESS, sep="_"))%>%
##    dplyr::select(NEW_BARCODE, bti) 
