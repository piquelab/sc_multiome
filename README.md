# Project GxE single multiomic analysis

We performed a series of  analyses in the following folders,
- The folder `analysis_multiome_JW_2023-06-03` contains script for displaying results and plots for manuscript.   
- The folder `analysis_multiomic_2023-03-27` contains the script for Baseline analyses By Mohammed. 
- The folder `genetic_variant_annot` contains the script for annotation genetic variants derived from sc_multiome data. This annotation files would be  used for `torus` and `dap-g` in eQTL analysis
- The folder `gtex_v8` contains the scripts for Whole blood eQTL analysis, including the folder `torus` for enrichment analysis and `dap-g` for fine-mapped eQTLs. 
- The folder `gwas` contains the scripts for identification risk genes, including the folder `gwas_prepare` for preparing correct format of gwas data, the folder `enloc_analysis` for   
  colocalization analysis and the folder `twas_analysis_SMR` for TWAS analysis.    

We used "BH" from p.adjust to calculate FDR as default approach.    
In default, we used all the genes passing filtering conditions to calculate FDR and after that we extracte protein coding genes. 
