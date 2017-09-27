#!/bin/bash

# author: Emily Tsang

set -o nounset -o errexit -o pipefail

# little script to run the commands for getting multiple primers based on the initial site list
# to get this to happen faster, run the python design script with sub batches of the sites.
# will combine them afterwards
# would be great to parallelize the python script instead, but don't have time to implement that.

## The variable in caps on the RHS should be set globally, e.g. in bashrc ## 
outdir=${EGTEXDIR}/primerDesign 
#################################################################

scriptdir=`dirname \$(readlink -f "\$0")`

# take one argument as input which is the pass name
if [ $# -ne 4 ]; then
    echo "usage: primer_design_pass.sh pass_name nparallel start_count master_file"
    exit
fi
pass=$1
nsplit=$2
start=$3
masterfile=$4

# make output directories
for (( i=1; i<=$nsplit; i++ ))
do
    mkdir ${outdir}/${pass}${i}
    inputfile=${outdir}/${pass}${i}/selected.variants.${pass}${i}.txt
    if [ $i -lt $nsplit ]; then
	cat $masterfile | head -n $((i*100)) | tail -n 100 > $inputfile
    else
	cat $masterfile | tail -n +$((i*100-99)) > $inputfile
    fi

    python ${scriptdir}/run_yamRTPCR_cDNA.py --poolSize 11 --numberPools 9 --counterStart $((start + (i-1)*10 - (i-1)*1)) --geneAnnotation UCSC_TableBrowser_gencode19_selected_exons UCSC_TableBrowser_gencode19_selected_transcripts --directory ${outdir}/${pass}${i} --inputSiteFile $inputfile &> ${outdir}/${pass}${i}/log.txt &
done

