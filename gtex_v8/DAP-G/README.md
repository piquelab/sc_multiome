### Run DAP-G for eQTL mapping in GTEx Whole blood  

whole blood sbams files, `link -s /wsu/home/groups/piquelab/gtex_v8/Tissues/Whole_Blood/sbams sbams`. 
split geneList files, 
```
cd geneList
split -l 500 geneList.txt -d -a 3 splitGene
cd ..
ls ./geneList |grep split > geneList_files.txt
```
run dap-g script `bash 3_run_dap-g.submit.sh`

We used joint annotation of condition response motifs for torus and the prior file from torus for fine-mapping dap-g. `cp ../torus/torus_output/Whole_Blood.dump.prior . `. 
 