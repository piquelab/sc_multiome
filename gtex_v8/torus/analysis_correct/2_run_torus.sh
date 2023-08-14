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

sbatch -q primary --mem=100G --time=2-23:00:00 -n 1 -N 1-1 --job-name=torus0.1_WBL_union --output=slurm_torus_union.output --wrap "
  module load misc;
  torus -d ../torus_input/Whole_Blood.eQTL.txt.gz \
    -smap ../torus_input/zzz_snp.map.gz \
    -gmap ../torus_input/zzz_gene.map.gz \
    -annot ./torus_input/zzz2_union_torus.annot.gz \
    -est > ${outdir}Union_WBL.est \
    -dump_prior ${outdir}Union_WBL.dump.prior"

### End

