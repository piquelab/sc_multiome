#!/bin/bash

cd $PWD

## create directory
outdir=./dap_rst_dir/
if [ ! -d ${outdir} ]; then
   mkdir -p ${outdir}
fi

cat geneList.txt | \
while read gene; do
   ##
   if [ -f ../DAP-G/dap-g_outs/Whole_Blood/${gene}.rst ]; then
       cp ../DAP-G/dap-g_outs//Whole_Blood/${gene}.rst ${outdir}
   fi
done

