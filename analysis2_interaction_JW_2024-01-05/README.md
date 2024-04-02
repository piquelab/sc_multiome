# Interaction analysis of single multiome data

In the path `/rs/rs_grp_scgxemultiome/` we perform the new analysis for single multiome data,
- The folder `1_interaction` for interaction analysis of gene expression and chromatin accessibility;
  - Interaction analysis of gene expression
    - `1_interaction_RNA.R` includes the two kinds of interaction analysis, (1) lrt overall test between full model w/ interaction term and reduced model w/o, output `./1_inter_RNA.outs/1.2_lrt_results.rds`; (2) consider all the combination of cell types and treatment, output `2_full_RNA.dds.rds`. Also combine (2) analysis results from 48 contrasts.  
    - `1.1_submit.sh` and `1.1_extract_results.R`, extract the results of (2) interaction analysis for 48 contrasts, output `2.2_each_pair.results.rds`. 
    - `3.1_RNA_compare.R` compare other cell types to the contrast cell type (CD4 Naive cells) for each treatment, 42 contrasts (6 treatments $\times$ 7 cell types). 
    - `3.2_RNA_summary.R` summarize the results 
  - Interaction analysis of chromatin accessibility
    - `2_interaction_ATAC.R` includes two kinds of interaction analysis, (1) lrt overall test between full model w/ interaction term and reduced model w/o, output `./2_inter_ATAC.outs/1.2_lrt_results.rds`; (2) consider all the combination of cell types and treatment, output `./2_inter_ATAC.outs/2_full_ATAC.dds.rds`. Also combine (2) analysis results from 48 contrasts, output `./2_inter_ATAC.outs/2.2_each_pair.results.rds`
    - `2.1_submit.sh` and `2.1_extract_results.R`, extract the results of (2) interaction analysis for 48 contrasts, output `./2_inter_ATAC.outs/dds_results/`.
    - `4.1_ATAC_compare.R` compare other cell types to the contrast cell type (CD4 Naive cells) for each treatment, 42 contrasts (6 treatments $\times$ 7 cell types). 
    - `4.2_ATAC_summary.R` summarize the results            
- The folder `2_estimate_var` for estimation of variance by cell types, treatment and individual for gene expression and chromatin accessibility.
 
