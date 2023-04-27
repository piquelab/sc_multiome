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

sbatch -q primary --mem=100G --time=2-23:00:00 -n 1 -N 1-1 --job-name=torus0.1_WBL --output=slurm_torus.output --wrap "
  module load misc;
  torus -d ./torus_input/Whole_Blood.eQTL.txt.gz \
    -smap ./torus_input/zzz_snp.map.gz \
    -gmap ./torus_input/zzz_gene.map.gz \
    -annot ./torus_input/zzz_torus.annot.gz \
    -est > ${outdir}Whole_Blood.est \
    -dump_prior ${outdir}Whole_Blood.dump.prior"

### End

