#!/bin/bash/

## cp files from GTEx project to the working directory
### JW, 6-15-2023

cat traits_of_interest.txt | \
while read trait; do
   ##
   cp /wsu/home/groups/piquelab/gtex_v8/gtex_gwas/${trait}.*vcf.gz* ./gwas_source/ &
   echo ${trait}
done 
 
