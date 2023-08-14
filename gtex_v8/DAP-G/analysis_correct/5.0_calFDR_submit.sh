#!/bin/bash/


cd $PWD
 
cat ../geneList_files.txt |
while read geneFile; do
   ##
   sbatch -q express -p erprp --mem=4G --time=1-00:00:00 -N 1-1 -n 1 --job-name=calFDR_${geneFile} --output=slurm_calFDR_${geneFile}.output --wrap "
   module load R;
   R CMD BATCH --no-save --no-retore '--args ${geneFile}' 5.0_calFDR.R calFDR_${geneFile}.Rout"
   echo ${geneFile}
   sleep 0.5;
done 
