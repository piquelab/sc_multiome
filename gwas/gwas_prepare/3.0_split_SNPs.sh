#!/bin/bash/

cd $PWD

cat traits_of_interest.txt | \
while read trait; do
###
  dir=./gwas_imputefile/${trait}
  # rm ${dir}/zzz_splitSNP*

  split -l 60000 ${dir}/zzz_allmissing.txt -d -a 4 ${dir}/zzz_splitSNP
  sleep 2;
  ls ${dir} | grep zzz_splitSNP > ./missing_snpFile/${trait}_missing_snpFile.txt 
  echo ${trait} 
##
done
