#!/bin/bash
#SBATCH -q primary
#SBATCH --mem=120G
#SBATCH --time=24:00:00
#SBATCH -N 1-1
#SBATCH -n 1
#SBATCH --job-name=vcf_dosage
#SBATCH --output=slurm_dosage.output


module load bcftools

outdir=./Example_genes_vcf/
if [ ! -d ${outdir} ]; then 
   mkdir -p ${outdir}
fi


cat gene_infor.txt | \
while read ens symbol chr pos1 pos2;
do

chr=chr${chr}
genfn=${ens}_${symbol}_${chr}
region=${chr}:${pos1}-${pos2}

echo ${ens} ${symbol} ${chr} ${pos1} ${pos2}

bcftools plugin dosage phASER_GTEx_v8_merged.vcf.gz -r ${region} | bgzip > ${outdir}${genfn}.vcf.gz

# sbatch -q express -p erprp --mem=120G --time=24:00:00 -N 1-1 -n 1 --job-name=${symbol}_dosage --output=dosage_${symbol}.output --wrap "
# module load bcftools;
# bcftools plugin dosage ../../phASER_GTEx_v8_merged.vcf.gz -r ${region} | bgzip > ${genfn}.vcf.gz"

done


 
