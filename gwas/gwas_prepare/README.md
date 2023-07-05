# prepare gwas files for the following analysis

We used the gwas data in the path `/wsu/home/groups/piquelab/gtex_v8/gtex_gwas/`. The target traits for analysis is `traits_of_interest.txt`. The SNPs list  in WBL_GTEx is `WBL_snplist`.
In this folder, we prepare three specific files (1) for torus enrich analysis, (2) PIP files for colocalization analysis and (3) for TWAS analysis.  
We firstprepare torus format file using the script `2_torus_submit.sh` and `2_torus_format.R`.    
Next we impute the gwas summary data in the following steps:
- First get missing SNPs that appeared in GTEx but don't in the gwas summary data using the script `3.0_missing_SNPs.R`. 
- Then break down the missing SNP files into bundle file with 60000 SNPs for each trait using the script `3.0_split_SNPs.sh`. 
- `3_impute.R` and `3_impute_submit.sh` impute the missing SNPs using the closest SNP around $\pm$ 10 kb region.
-  The script `4_concate.sh` and `4_combine.R` combined all the bundle files and the orginal gwas data into the imputed files, which are the final files used for the downstream analysis. 
  
  