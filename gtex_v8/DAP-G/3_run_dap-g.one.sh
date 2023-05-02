#!/bin/bash
#SBATCH -q primary
##SBATCH -p erprp
#SBATCH --mem=40G
#SBATCH --time=2-01:00:00
#SBATCH -N 1-1
#SBATCH -n 4


outdir=./dap-g_outs/Whole_Blood
if [ ! -d ${outdir} ]; then
   mkdir -p ${outdir}
fi


cat ./geneList/${geneFile} | \
while read ENSG;
do
   ##
   if [ -f ./Whole_Blood/sbams/${ENSG}.sbams.dat ] && [ -f ./Whole_Blood.dump.prior/${ENSG}.prior ]; then 
      echo ${ENSG}
      /wsu/home/ha/ha21/ha2164/Bin/dap-g.static -d ./Whole_Blood/sbams/${ENSG}.sbams.dat \
         -p ./Whole_Blood.dump.prior/${ENSG}.prior \
         -t 4 \
         -o  ${outdir}/${ENSG}.rst -l ${outdir}/${ENSG}.log 
    fi
##
done
