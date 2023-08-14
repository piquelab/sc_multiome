# Using SMR approach to conduct TWAS

We used three approaches to select the gene representative SNP:
- Pick the minimum pvalue SNP using FastQTL result. The FastQTL result is from `/wsu/home/groups/piquelab/gtex_v8/Data/fastqtl/`, which is by William.
- Pick the maximum PIP SNP using DAP-G result with integration of multiple condition annotations. 
- Pick the maximum PIP SNP using DAP-G with baseline annotation. This result is from `sc-atac` project in the `./genetic_analysis_GTExV8/DAP-G/5_summary.outs/`. 

The folder `analysis_correct` is the final analysis using the newly generated response motifs
      
