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

bcftools query phASER_GTEx_v8_merged.vcf.gz -r ${region} -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\n' | bgzip > ${outdir}${genfn}.SNPid.txt.gz


done


 
