#!/bin/bash/

for i in {1..48}; do

sbatch -q express -p erprp --mem=20G --time=12:00:00 -n 1 -N 1-1 --job-name=results_${i} --output=slurm_results_${i}.out --wrap "

module load R;

R CMD BATCH --no-save --no-restore '--args ${i}' 2.1_extract_results.R 2.1_extract_${i}.Rout"

echo ${i}

sleep 1;

done
