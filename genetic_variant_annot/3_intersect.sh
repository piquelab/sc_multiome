#!/bin/bash

module load bedtools 

ls ./Peak_bed/ |sed 's/_peak.bed//'| uniq | \
while read oneMCl;
do

echo ${oneMCl}
bedtools intersect -a snp2.bed.gz -b ./Peak_bed/${oneMCl}_peak.bed -c|awk '{print $1, $2, $3, $4}' OFS='\t'|bgzip > ./Peak_bed/SNP_${oneMCl}_peak.bed.gz &
###
done
