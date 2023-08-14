#!/bin/bash


###
### create output directory

outdir=./torus_output/
if [ ! -d ${outdir} ]; then
   ##
   mkdir -p ${outdir}
fi


###
### submit jobs
cat annot_ls.txt | \
while read ii; 
do

echo ${ii} 

sbatch -q primary --mem=50G --time=2-23:00:00 -n 1 -N 1-1 --job-name=torus0.1_${ii} --output=slurm_torus_${ii}.output --wrap "
  module load misc;
  torus -d ../../torus_input/Whole_Blood.eQTL.txt.gz \
    -smap ../../torus_input/zzz_snp.map.gz \
    -gmap ../../torus_input/zzz_gene.map.gz \
    -annot ./torus_input/${ii}.torus.gz \
    -est > ${outdir}${ii}_WBL.est \
    -dump_prior ${outdir}${ii}_WBL.dump.prior"
sleep 1;

done


### End

