#!/bin/bash/

cd $PWD

cat traits_not_finished.txt | \
while read trait; do

###trait="Hypertension"
# for trait in Hypertension UKB_BMI;

cat ./missing_snpFile/${trait}_missing_snpFile.txt | \
while read snpFile; do
   ##
   sbatch -q express -p erprp --mem=20G --time=3-00:00:00 -N 1-1 -n 1 --job-name=impute_${trait}_${snpFile} --output=slurm_impute_${trait}_${snpFile}.output --wrap "
   module load R; 
   R CMD BATCH --no-save --no-restore '--args ${trait} ${snpFile}' 3_impute.R impute_${trait}_${snpFile}.Rout"
   echo ${trait} ${snpFile}
   sleep 1;
done 

done
