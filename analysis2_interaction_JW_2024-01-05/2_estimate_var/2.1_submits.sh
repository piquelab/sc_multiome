#!/bin/bash/

cd $PWD

cat ./2_vars_ATAC.outs/peakList_files.txt | \
while read geneFile; do 
###
   sbatch -q primary --mem=10G --time=12:00:00 -N 1-1 -n 1 --output=slurm_vars_${geneFile}.outs --wrap "
   module load R;
   R CMD BATCH --no-save --no-retore '--args ${geneFile}' 2.1_var_ATAC.R vars_${geneFile}.Rout"
   echo ${geneFile}
   sleep 0.5;
###
done    
