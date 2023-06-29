###
library(tidyverse)
library(data.table)
library(openxlsx)
library(ggrastr)
library(cowplot)
options(scipen=13)

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
### gwas data 
outdir <- "./gwas_0/"
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)


###
###
traits <- sort(read.table("traits_of_interest.txt")$V1)
for ( ii in traits){
    ###
    fn <- paste("./gwas_source/", ii, ".gambit.vcf.gz", sep="")
    summ <- fread(fn, header=F, data.table=F)
    summ <- summ%>%
       mutate(id_b38=paste(V1, V2, V4, V3, "b38", sep="_"),
           pval=pnorm(abs(V7), lower.tail=F)*2)
    
    ### output
    summ2 <- summ%>%dplyr::select(id_b38, chr=V1, pos=V2, zscore=V7, pval)
    opfn <- gzfile(paste(outdir, ii, "_gwas.txt.gz", sep=""))
    write.table(summ2, file=opfn, row.names=F, col.names=T, quote=F)

    cat(ii, nrow(summ2), "\n")
}    

### End

## id_new <- as.character(summ$id_b38)
## snp <- fread("WBL_snplist.txt", header=F)$V1
