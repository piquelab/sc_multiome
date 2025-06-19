###
###
library(tidyverse)
library(data.table)
library(openxlsx)
library(annotables)
##

rm(list=ls())

## outdir <- ""
## if ( ! file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)

###
### Criteria of filtering gene
### 1-FDR_INTACT<0.05
### 2-annotated in Response motif
### 3-PIP>0.2 and FDR_twas<1e-03


fn <- "traits_name_df.xlsx"
traits_df <- read.xlsx(fn)

 
res_all <- NULL
for (i in 1:11){
   ###    
   trait <- traits_df$traits_long[i]
   trait2 <- traits_df$traits_short[i]
    
   ### annotation of sigs 
   fn <- paste("../3_examples_output/", trait, "/", trait, "_summ_annot_fdr0.05.xlsx", sep="")
   res_anno <- read.xlsx(fn)
   res_anno <- res_anno[,c(1,5:33)]
   ### all results for genes
   fn <- paste("../3_examples_output/", trait, "/", trait, "_comb_infor.txt", sep="")
   res <- read.table(fn, header=T)%>%filter(FDR_intact<0.05, nannot>1)
   ### combine
   res2 <- res%>%inner_join(res_anno, by="Gene")%>%arrange(FDR_intact)

   ### filtering, PIP and FDR_twas  
   res3 <- res2%>%filter(PIP_conditions>0.2, FDR_twas<1e-03)
   cat(i, trait2, nrow(res3), "\n")
 
   if ( nrow(res3)>0){
      ###    
      res3$traits <- trait2
      res_all <- rbind(res_all, res3)
   }    
}    

opfn <- "1_candidate_genes.xlsx"
write.xlsx(res_all, file=opfn)


#######################
### rename symbol
#########################

res <- read.xlsx("1_candidate_genes.xlsx")


anno <- grch38%>%filter(ensgene%in%res$Gene, chr%in%as.character(1:22))%>%dplyr::select(ensgene, symbol)
gene_name <- anno$symbol
names(gene_name) <- anno$ensgene

res <- res%>%mutate(symbol=gene_name[Gene])

write.xlsx(res, "1_candidate_genes.xlsx")



##########################
### get snp rs id
#########################


library(biomaRt)

res <- read.xlsx("1_candidate_genes.xlsx")

id_b38 <- sort(unique(res$id_b38))
x  <- str_split(id_b38, "_", simplify=T)
x2 <- data.frame("id_b38"=id_b38,
   chr_pos=paste(gsub("chr", "", x[,1]), x[,2], x[,2], sep=":"),
   allele0=paste(x[,3], "/", x[,4], sep=""))


pos <- unique(x2$chr_pos)
 
ensembl <- useEnsembl(biomart="snps", dataset="hsapiens_snp")

###
### get snp rs
pos_df <- map_dfr(1:length(pos), function(i){  
   ##
   ii <- pos[i] 
   cat(i, ii, "\n") 
   df0 <- try(getBM(attributes=c("refsnp_id", "allele", "chrom_start"),
      filters="chromosomal_region",
      values=ii,
      mart=ensembl), silent=T)
   ##
   if ( class(df0)!="try-error" & nrow(df0)>0){
       df0$chr_pos <- ii
       df0
    }
})    
opfn <- "0_snpinfor.txt"
write.table(pos_df, file=opfn, quote=F, row.names=F)

pos_df2 <- pos_df%>%mutate(pos2=as.integer(gsub(".*:", "", chr_pos)))
pos_df2 <- pos_df2%>%filter(chrom_start==pos2)


###
### combine rs id with 
x3 <- x2%>%full_join(pos_df2, by="chr_pos")%>%
    dplyr::select(id_b38, refsnp_id, allele0, allele)%>%
    drop_na(id_b38, refsnp_id, allele0, allele)

###
### get rs id for each snp after remove replicates

rs_df <- map_dfr(unique(x3$id_b38), function(ii){
   ##
   df0 <- x3%>%filter(id_b38==ii)
   if (nrow(df0)>1){
       cat(ii, nrow(df0), "\n")
       ### select column
       isel <- sapply(1:nrow(df0), function(i){
           ##
           s0 <- unlist(str_split(df0$allele0[i], "/"))
           s1 <- unlist(str_split(df0$allele[i],"/"))
           all(s0%in%s1)
       })
       df2 <- df0[isel,1:2]       
   }else{
       df2 <- df0[,1:2]
   }
   df2 
})


###
### add rs id into candidate gene files and output
res2 <- res%>%left_join(rs_df, by="id_b38")
write.xlsx(res2, file="1.2_candidate_genes.xlsx")
