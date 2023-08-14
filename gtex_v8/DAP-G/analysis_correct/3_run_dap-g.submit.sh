#!/bin/bash/

cat ../geneList_files.txt | \
while read geneFile;
do
echo ${geneFile}
sbatch --export=geneFile=${geneFile} --job-name=dap-g_${geneFile} --output=slurm.WBL_${geneFile}_out 3_run_dap-g.one.sh
sleep 0.5;
done
## End 
