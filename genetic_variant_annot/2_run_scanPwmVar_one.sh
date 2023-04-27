#!/bin/bash
#SBATCH -q primary
#SBATCH --mem=20G
#SBATCH --time=2-01:00:00
#SBATCH -N 1-1
#SBATCH -n 1


module load ucscbrowser/2019-05-25

## create output directory
outdir=./annot_jaspar2022/
if [ ! -d ${outdir} ]; then
  mkdir -p ${outdir}
fi


cat ./Motif_file/${motifFile} | \
while read motif;
do
echo ${motif}

if  [ -f ./JASPAR2022_core/${motif}.jaspar ]; then
   scanPwmVar ./JASPAR2022_core/${motif}.jaspar -j -base=2.0 /wsu/home/groups/piquelab/data/RefGenome/hg38.2bit snp.bed.gz|\
     awk -v OFS='\t' '$5>10||$8>10 {print $1, $11+1, $11, $1"_"$11+1"_"$12"_"$13}' | bgzip > ./annot_jaspar2022/allsnp_${motif}.bed.gz

fi

done
