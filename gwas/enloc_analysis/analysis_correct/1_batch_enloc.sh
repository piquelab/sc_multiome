#!/bin/bash

cd $PWD

###loop by tissue
cat ../traits_of_interest.txt | \
while read trait;
do 
   echo ${trait}
   # outdir=./enloc_output/${trait}/
   # if [ ! -d ${outdir} ]; then
   #    ##
   #    mkdir -p ${outdir}
   # fi
   ###
   sbatch -q express -p erprp --mem=40G --time=2-10:00:00 -n 1 -N 1-1 --job-name=enloc_${trait} --output=slurm_enloc_${trait}.out --wrap "
   /wsu/home/ha/ha21/ha2164/Bin/fastenloc.static -eqtl ./eQTL_results/WBL_enloc.eqtl.vcf.gz -go ../gwas_PIP/${trait}.pip.gz -total_variants 9428357 -prefix ./enloc_output/${trait}"
   ###
done
      
### Aorta, 1,078,801
### Coronary, 1,057,064 
### Tibial, 1,034,996
