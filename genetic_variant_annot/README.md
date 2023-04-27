# The annotation of genetic variant for torus 

## First step, get motif list 
We get motif 692 core motifs from JASPAR2022 by running the script `1_run_getJaspar2022.R`. And split 692 motif into 18 files boundle
```ruby
mkdir Motif_file
cd Motif_file
split ../motifList2022.txt -l 40 -d -a 3 splitMotif
cd ..
ls ./Motif_file/ |grep splitMotif > motifList_files.txt
```

## Second step, annotate genetic variants in TF binding sites  
- we first generated bed file for genetic variants from GTEx v8 project using the 
```ruby
## extract snp information file from GTEx v8 project
bcftools query /wsu/home/groups/piquelab/gtex_v8/GTEx_Analysis_v8_phASER/phASER_GTEx_v8_merged.vcf.gz -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' | bgzip > gtex_v8_snpinfor.txt.gz &

## generate bed file
zcat gtex_v8_snpinfor.txt.gz |awk -v OFS='\t' '{print $1,$2,$4,$5}'|bgzip > snp.bed.gz &

## the script for generating bed file for bedtools 
zcat snp.bed.gz |awk -v OFS='\t' '{print $1, $2, $2+1}'|bgzip > snp2.bed.gz & 

```
