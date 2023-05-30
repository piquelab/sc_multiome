#!/bin/bash

cd $PWD


sbatch -q primary --mem=80G --time=1-23:00:00 -n 1 -N 1-1 --job-name=prepare_WBL --output=slurm_WBL.out --wrap "
   perl summarize_dap2enloc.pl -dir dap_rst_dir -vcf ./phASER_GTEx_v8_merged.vcf.gz |gzip - > WBL_enloc.eqtl.vcf.gz "
