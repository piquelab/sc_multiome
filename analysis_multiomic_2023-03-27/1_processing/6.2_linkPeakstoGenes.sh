#!/bin/bash
#04-29-2021
#Mohammed Husain Bharmal
# macs2 peak calling for each sample

mkdir -p 6.2_linkPeakstoGenes_2.2_slurm

##mkdir -p slurm


##cat comb_list_try.txt  |  while read sample; do
##  sbatch -q primary  -N1-1 -n 2 --mem=20G -t 1000 --job-name=${sample} -o "6.2_linkPeakstoGenes_slurm/slurm_${sample}_-%A_%a.out" --wrap "Rscript 6.2_linkPeakstoGenes.R $sample"
##sleep 0.5
##done


cat comb_list_test.txt  |  while read sample; do
echo ${sample}
sbatch -q primary  -N1-1 -n 2 --mem=240G -t 2000 --job-name=${sample} -o "6.2_linkPeakstoGenes_2.2_slurm/slurm_${sample}" --wrap "Rscript 6.2_linkPeakstoGenes_2.R $sample"
sleep 0.5
done
