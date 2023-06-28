###
library(tidyverse)
library(data.table)
library(openxlsx)
library(ggrastr)
library(cowplot)

###
### traits list
fn <- "../traits/traits_of_interest.xlsx"
x <- read.xlsx(fn)

trait2 <- x%>%filter(Inclusion==1)%>%pull(traits)

opfn <- "traits_of_interest.txt"
write.table(trait2, opfn, row.names=F, col.names=F, quote=F)

opfn <- "traits_of_interest.xlsx"
x2 <- x%>%filter(Inclusion==1)%>%dplyr::select(traits, ngene_ptwas)
write.xlsx(x2, file=opfn, overwrite=T)



###
###
outdir <- "./gwas_0/"
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)


traits <- sort(read.table("traits_of_interest.txt")$V1)


traitDF <- data.frame(trait_new=traits[c(9, 1, 8, 12)],
     trait_old=c("UKB_asthma", "CAD", "Hypertension", "Height"))

###
### check the traits
fig_ls <- lapply(1:nrow(traitDF), function(i){
    ####
    inew <- traitDF$trait_new[i]
    fn <- paste("./gwas_source/", inew, ".gambit.vcf.gz", sep="")
    summ <- fread(fn, header=F, data.table=F)
    summ <- summ%>% mutate(id_b38=paste(V1, V2, V4, V3, "b38", sep="_"))

    ## ,
    ##         pval=pnorm(abs(V7), lower.tail=F)*2)
    ###
    summ2 <- summ%>%dplyr::select(id_b38, z_new=V7)

 
## id_new <- as.character(summ$id_b38)
## snp <- fread("WBL_snplist.txt", header=F)$V1
    iold <- traitDF$trait_old[i]
    fn <- paste("/wsu/home/groups/piquelab/gtex_v8/enloc_analysis/gwas/", iold, ".torus.zval.gz", sep="")
    old <- fread(fn, header=F, data.table=F)
    old2 <- old%>%dplyr::select(id_b38=V1, z_old=V3)

    cat(iold, "\n")
    
    DF <- summ2%>%inner_join(old2, by="id_b38")
    ##
    p <- ggplot(DF, aes(x=z_old, y=z_new))+
       rasterize(geom_point(size=0.1, color="grey"), dpi=100)+
       theme_bw()+
       ggtitle(iold)+
       xlab(bquote(~italic(Z)~"score from enloc_analysis"))+
       ylab(bquote(~italic(Z)~"score from gtex_gwas"))+         
       theme(axis.text=element_text(size=10),
             axis.title=element_text(size=10),
             plot.title=element_text(size=12, hjust=0.5))
    p
})

##
figfn <- "./torm/Fig1_comb.png" 
png(figfn, width=680, height=680, res=100)
plot_grid(plotlist=fig_ls, ncol=2)
dev.off()



###
### compare two asthma datasets

## william sent to me
fn <- "imputed_UKB_20002_1111_self_reported_asthma.txt.gz"
data1 <- fread(fn, header=T, data.table=F)
data1 <- data1%>%dplyr::select(id_b38=panel_variant_id, z_imput=zscore) 

### from gtex_gwas folder
fn <- "./gwas_source/UKB_20002_1111_self_reported_asthma.gambit.vcf.gz"
data2 <- fread(fn, header=F, data.table=F)
data2 <- data2%>% mutate(id_b38=paste(V1, V2, V4, V3, "b38", sep="_"))
data2 <- data2%>%dplyr::select(id_b38, z_new=V7)
 
DF <- data1%>%inner_join(data2, by="id_b38")
##
p <- ggplot(DF, aes(x=z_imput, y=z_new))+
   rasterize(geom_point(size=0.1, color="grey"), dpi=100)+
   theme_bw()+
   xlab(bquote(~italic(Z)~"score from imputed"))+
   ylab(bquote(~italic(Z)~"score from gtex_gwas"))+ 
   ggtitle("Asthma")+
       theme(axis.text=element_text(size=10),
             axis.title=element_text(size=10),
             plot.title=element_text(size=12, hjust=0.5))
###
figfn <- "./torm/Fig1.2_asthma.png"
png(figfn, width=420, height=420, res=100)
p
dev.off()
