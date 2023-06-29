#!/bin/bash/

cd $PWD

cat traits_of_interest.txt | \
while read trait; do

###
outdir2=./gwas_imputefile/${trait}/
cat ${outdir2}/gwas_*.txt.gz > gwas_concate.txt.gz 

echo ${trait}

done 
