#!/bin/bash


###
### 
cd $PWD
cat traits_ls.txt | \
while read trait;
do 
echo ${trait}
###
sbatch -q express -p erprp --mem=20G --time=2-01:00:00 -n 1 -N 1-1 --job-name=${trait}_torus_PIP --output=${trait}_torus_PIP.out --wrap "
   cd $PWD;
   module load misc;
   torus --load_zval -d ./gwas/${trait}.torus.zval.gz -dump_pip  ./gwas/${trait}.pip "
###
done
