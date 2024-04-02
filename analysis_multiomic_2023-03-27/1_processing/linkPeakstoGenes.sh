#!/bin/bash


#FILE=3_macs2.outs/${sample}_peaks.narrowPeak
#if [[ -f "$FILE" ]]; then
#    echo "$FILE exists."
#    exit 0
#fi


#############
#############
sample=$1

##bedtools intersect -a ./3_macs2.outs/merged.filtered.bed.gz  -b ./2_clean_bam/${sample}_clean.bam -c -sorted | bgzip >./4_merging_and_intersecting.outs/${sample}_intersect.bed.gz

Rscript 
