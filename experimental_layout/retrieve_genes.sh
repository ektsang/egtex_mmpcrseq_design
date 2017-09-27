#!/bin/bash

# author: Emily Tsang

## # Wrapper script to run retrieve.genes.py

set -o nounset -o errexit -o pipefail

## The variables in caps on the RHS should be set globally, e.g. in bashrc ##
gtf=$GENCODE19
dir=${EGTEXDIR}/geneSelection
#################################################################

scriptdir=`dirname \$(readlink -f "\$0")`

python ${scriptdir}/retrieve_genes.py $gtf ${dir}/selected.genes.txt | sort -k1,1 -k2,2n > ${dir}/selected.genes.exons.bed
python ${scriptdir}/retrieve_genes.py $gtf ${dir}/selected.genes.v2.txt | sort -k1,1 -k2,2n > ${dir}/selected.genes.exons.v2.bed
