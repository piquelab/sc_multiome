#!/bin/bash

cd $PWD

###loop by tissue
cat traits_ls.txt | \
while read trait;
do 
   echo ${trait}
   sbatch -q primary --mem=20G --time=1-10:00:00 -n 1 -N 1-1 --job-name=enloc_${trait} --output=slurm_${trait}.out --wrap "
   /wsu/home/ha/ha21/ha2164/Bin/fastenloc.static -eqtl ./eQTL_results/WBL_enloc.eqtl.vcf.gz -go ./gwas/${trait}.pip.gz -total_variants 11946126 -prefix ./enloc_output/${trait}"
   ###
done
      
### Aorta, 1,078,801
### Coronary, 1,057,064 
### Tibial, 1,034,996
