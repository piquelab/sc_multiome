# Colocalization analysis and INTACT analysis 

We use fastenloc software to do colocalization analysis. The three steps are as follows,
- Prepare eqtl files from `DAP-G` outputs in the directory `../../gtex_v8/enloc_prepare/` using `perl summarize_dap2enloc.pl -dir dap_rst_dir -vcf ./phASER_GTEx_v8_merged.vcf.gz |gzip - >fastenloc.eqtl.vcf.gz`.  
- Prepare gwas files using `bash torus_PIP.sh` in the directory `./gwas_analysis` and followed by `gzip xxx.pip`. 
- Run `./fastenloc.static -eqtl ./eQTL_results/WBL_enloc.eqtl.vcf.gz -go ./gwas/Asthma.pip.gz -total_variants 11946126 &`

We also run `2_summmary.R` R script to get the INTACT results, the file named `ALOFT_intact.txt`. Run INTACT under the `R version test_4.2.2`

We replicated colocalization analysis using new response motif annotation in the folder `analysis_correct`. 


