##
library(tidyverse)
library(DESeq2)
library(Matrix)
library(qvalue)

rm(list=ls())

outdir <- "./4_ATAC_summary.outs/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)




###
fn <- "./2_inter_ATAC.outs/contrast_list2.txt"
con_df <- read.table(fn, header=T)


fn <- "./2_inter_ATAC.outs/2.2_each_pair.results.rds"
results <- read_rds(fn)
## names(results)[8] <- "condition"
## opfn2 <- fn
## write_rds(results, opfn2)

###
### Effects
Bmat <- results%>%
    pivot_wider(id_cols=gene, names_from=condition, values_from=estimate, values_fill=NA)%>%
    column_to_rownames(var="gene")

###
### SE
SEmat <- results%>%
    pivot_wider(id_cols=gene, names_from=condition, values_from=stderror, values_fill=NA)%>%
    column_to_rownames(var="gene")

identical(rownames(Bmat), rownames(SEmat))
identical(colnames(Bmat), colnames(SEmat))



###
### compare effects in cell-typeX_treatment vs 0_CD4Naive_treatment
res0 <- results%>%filter(condition=="0_CD4Naive_caffeine")%>%dplyr::select(gene, baseMean)
identical(res0$gene, rownames(Bmat))

res_compare <- map_dfr(1:nrow(con_df), function(i){
    ##
    con1 <- con_df$con1[i]
    con0 <- con_df$con0[i]

    cat (con1, "\n")

    ### beta and se
    b1 <- Bmat[,con1]
    se1 <- SEmat[,con1]

    b0 <- Bmat[,con0]
    se0 <- SEmat[,con0]

    diff <- b1-b0
    se_diff <- sqrt(se1^2+se0^2)
    zscore <- diff/se_diff
    pval <- pnorm(abs(zscore), lower.tail=FALSE)*2

    ### fdr
    fdr <- p.adjust(pval, "BH")
    ## ii_nna <- !is.na(pval)
    ## fdr_ii <- p.adjust(pval[ii_nna]) ##qvalue(pval[ii_nna])$qvalues
    ## fdr[ii_nna] <- fdr_ii
    ###
    res0 <- res0%>%
        mutate(estimate=diff, stderror=se_diff, statistic=zscore, p.value=pval, p.adjusted=fdr,
               conditionX=con1, conditionY=con0)
    res0
})



res_compare <- res_compare%>%mutate(padj2=p.adjust(p.value, "BH"))

###
###
opfn <- paste(outdir, "1_MCls.diff.rds", sep="")
write_rds(res_compare, file=opfn)



###
### END

