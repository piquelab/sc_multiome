#!/bin/bash
#SBATCH -q primary
###SBATCH -p erprp
#SBATCH --mem=4G
#SBATCH --time=2-01:00:00
#SBATCH -N 1-1
#SBATCH -n 1


### output directory
outdir=./dap-g_outs/Whole_Blood

##
cat ../geneList/${geneFile} | \
while read ENSG;
do
##
   if [ -f ${outdir}/${ENSG}.rst ]; then
   ##
     echo ${ENSG}
      grep '\[' ${outdir}/${ENSG}.rst > ${outdir}/${ENSG}.model.out
      grep '((' ${outdir}/${ENSG}.rst > ${outdir}/${ENSG}.SNP.out
      grep '{'  ${outdir}/${ENSG}.rst > ${outdir}/${ENSG}.cluster.out
   fi
##
done





