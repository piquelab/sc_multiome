
library(JASPAR2022)
library(TFBSTools)
library(gdata)


###

outdir <- "./JASPAR2022_core/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)

### get motif data
pfm <- getMatrixSet(x=JASPAR2022, opts=list(species=9606, all_version=F))

### motif list
motif <- sort(as.character(names(pfm)))
write.table(motif, "motifList2022.txt", row.names=F, quote=F, col.names=F)

### write motif
for (i in motif){
    pfm0 <- pfm[[i]]
    opfn <- paste(outdir, i, ".jaspar", sep="")
    write.fwf(as.matrix(pfm0), opfn, sep="\t", colnames=F)
}    


## pwm
motif <- read.table("motifList2022.txt")$V1

pwm <- getMatrixSet(x=JASPAR2022, opts=list(specifies=9606, matrixtype="PWM"))

