# Analyis of gene expression, chromatin accessibility and TFs

The primary analysis script are completed by Mohammed and in the path `/nfs/rprdata/julong/sc_multiome/analysis_multiomic_2023-03-07`
The scripts in this directory are mainly used for plots in the pulication and deep analyses, including the following parts. Set the working directory(`/nfs/rprdata/julong/sc_multiome/analysis_multiome_JW_2023-06-03/`)  
1. Preprocessing analysis includes  main and supplementary Figures 
   - [1_processing_plots.R](1_processing_plots.R) for Figure1
   - The folder `1.2_ArchR_process` contains the script for calculating gene score using [ArchR](https://www.archrproject.com/)
2. Differential analysis (gene expression+chromatin accessibility)
   - [2.3_summary_diff.R](2.3_summary_diff.R) for main Figure2 and supplementary Figures in differential analysis (gene expression+chromatin accessbility)     
   - [2.3_MA_diff.R](2.3_MA_diff.R) for MA and qq plots in the differential analysis
   - [2.3_clusterProfiler.R](2.3_clusterProfiler.R) for enrichment analysis for DEGs  
3. Compare response changes between gene expression and chromatin accessibility
   - [3.2_compareRNAandATAC.R](3.2_compareRNAandATAC.R) for correlation analysis between RNA and ATAC
   - [3.2_track_plots.R](3.2_track_plots.R) for chromatin track plots
4. Compare the response changes betwee TF genes and TF activity
   - [4.2_motif_analysis.R](4.2_motif_analysis.R) for geting  the data of TF-genes and the pairwise heatmap
   - [4.2_TF-genes_cca.R](4.2_TF-genes_cca.R) CCA analysis for TF-genes and TF activity
5. The folder `5_LDA_analysis` contains the scripts for DLDA analysis
   - [1_prepare_data.R](./5_LDA_analysis/1_prepare_data.R) generate seurat object for each cell type
   - [2_RNA.DLDA.R](./5_LDA_analysis/2_RNA.DLDA.R) caculate 6 treatments' DLDA for each cell type
   - [3_summary_results.R](./5_LDA_analysis/3_summary_results.R) making scatter plots and Heatmap plots. 
     