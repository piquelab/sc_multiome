###
###
library(tidyverse)
library(data.table)
options(scipen=15)

rm(list=ls())

### parsing argument
args=commandArgs(trailingOnly=T)
if ( length(args)>0){
   ##
   trait <- args[1] 
   snpFile <- args[2]
}else{
   trait <- "CARDIoGRAM_C4D_CAD"
   snpFile <- "zzz_splitSNP0000"
}


###
### directory for outputs
outdir2 <- paste("./gwas_imputefile/", trait, "/", sep="")


###
### gwas data
fn <- paste("./gwas_0/", trait, "_gwas.txt.gz", sep="")
summ <- fread(fn, header=T, data.table=F)


###
### missing SNP files
fn <- paste(outdir2, snpFile, sep="")
snp <- read.table(fn)$V1
cvt <- str_split(snp, "_", simplify=T)
bed <- data.frame(id_b38=snp, chr=cvt[,1], pos=cvt[,2])

### 1859556, millions SNPs missing
## fn <- paste(outdir2, "zzz_allmissing.txt", sep="")
## x <- read.table(fn)$V1 


###
### imputation 
nsnp <- nrow(bed)
time0 <- Sys.time()
tmp <- lapply(1:nsnp, function(i){
   ##
   if (i%%10==0) cat(i,"\n") 
   chr_i <- bed$chr[i]
   pos_i <- as.integer(bed$pos[i])

   ## +/- 10 kb 
   s0 <- pos_i-1e+04
   s1 <- pos_i+1e+04
   
   tmp2 <- summ%>%filter(chr==chr_i, as.numeric(pos)>=s0, as.numeric(pos)<=s1)
   if ( nrow(tmp2)!=0){
      dtss <- abs(as.numeric(tmp2$pos)-pos_i)
      imin <- which.min(dtss)
      tmp3 <- data.frame(id_b38=bed$id_b38[i], chr=chr_i, pos=pos_i,
                         zscore=tmp2$zscore[imin], pval=tmp2$pval[imin],
                         id_b38_close=tmp2$id_b38[imin])
   }else{
      tmp3 <- NULL
   }   
   tmp3
})

time1 <- Sys.time()
diff0 <- difftime(time1, time0, units="secs")
cat(diff0, "\n")

##
tmp <- do.call(rbind, tmp)


###
### output
ii <- gsub("zzz_", "",  snpFile)
gfn <- gzfile(paste(outdir2, "gwas_", ii, ".txt.gz", sep=""))
write.table(tmp, file=gfn, row.names=F, col.names=F, quote=F)


###END
