#!/bin/bash/


cd $PWD

cat ../traits_of_interest.txt | \
while read trait; do
   ###
   sbatch -q express -p erprp --mem=120G --time=3-00:00:00 -N 1-1 -n 1 --job-name=SMR_${trait} --output=slurm_SMR_${trait}.output --wrap "
   module load R;
   R CMD BATCH --no-save --no-retore '--args ${trait}' 1_SMR_WBL.R zzz_SMR_${trait}.Rout"
   echo ${trait}
   sleep 1;
done
