#!/bin/bash/

cat motifList_files.txt | \
while read motifFile;
do 

echo ${motifFile}
sbatch --export=motifFile=${motifFile} --output=slurm_${motifFile}_out 2_run_scanPwmVar_one.sh

sleep 1;

done

##End
