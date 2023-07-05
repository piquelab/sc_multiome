#!/bin/bash


###
### 
cd $PWD
cat traits_of_interest.txt | \
while read trait;
do 
echo ${trait}
###
sbatch -q express -p erprp --mem=30G --time=2-01:00:00 -n 1 -N 1-1 --job-name=${trait}_torus_PIP --output=slurm_${trait}_torus_PIP.out --wrap "
   cd $PWD;
   module load misc;
   torus --load_zval -d ./gwas/${trait}_torus.txt.gz -dump_pip  ./gwas_PIP/${trait}.pip "
###
done
