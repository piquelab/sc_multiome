#!/bin/bash/

cd $PWD

###################################################################
### torus enrichment analysis for combine2 and Union annotation ###
###################################################################

outdir=./torus_output/
if [ ! -d ${outdir} ]; then
   mkdir -p ${outdir}
fi


## loop for each cell-type
cat traits_ls.txt | \
while read ii;
do

echo ${ii}

sbatch -q primary --mem=80G --time=5-23:00:00 -N 1-1 -n 1 --job-name=torus0.1_${ii} --output=slurm0.1_${ii}_torus.out --wrap "
  module load misc;               
  torus --load_zval -d ./gwas_input/${ii}.torus.zval.gz \
     -annot ./torus_input/zzz_torus.annot.gz \
     -est > ${outdir}${ii}.est "
##  
sleep 1;
done
### End 


 

