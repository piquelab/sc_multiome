## How to get eQTL source files
The summary data are from `cp /wsu/home/groups/piquelab/gtex_v8/Data/summary_data/Whole_Blood.allpairs.txt.gz .`       
We get the SNP list by `zcat Whole_Blood.allpairs.txt.gz |grep -v gene |cut -f 2 |sort|uniq > snpList.txt &`       
We get  the gene list by `zcat Whole_Blood.allpairs.txt.gz` |grep -v gene |awk '{print $1}' |sort |uniq > geneList.txt &`.   
We describe the structure of folders,
- Currect directory is default setting, define response motif for each condition at the top10% of response change and jointly estimate the enrichment of eQTL signals. 
- The folder `analysis2_th0.2`, define response motif for each condition at the top 20% of response change and jointly estimate the enrichment of eQTL signals. We don't use it in the paper.
- The folder `analysis3_each, define response motif for each condition at the top10% of response change and separately estimate the enrichment of eQTL signals.     
- The folder `analysis_correct`, I generated new response motifs and annotation used for the final analysis. This results will be used in the paper. 