#!/bin/bash

# author: Emily Tsang

## Wrapper to run retrieve.sites.py

set -o nounset -o errexit -o pipefail

## The variables in caps on the RHS should be set globally, e.g. in bashrc ##
# 
dir=${EGTEXDIR}/sampleInfo/genotypes
#################################################################

scriptdir=`dirname \$(readlink -f "\$0")`

python ${scriptdir}/retrieve_sites.py -o ${dir}/OMNI_V6_86ind_SNP_PASS_ONLY.bed ${dir}/OMNI_V6_86ind_SNP_PASS_ONLY.vcf.gz &

python ${scriptdir}/retrieve_sites.py -o ${dir}/WGS_V7_79ind_SNP_PASS_ONLY.bed ${dir}/WGS_V7_79ind_SNP_PASS_ONLY.vcf.gz &

python ${scriptdir}/retrieve_sites.py -o ${dir}/WES_V7_80ind_SNP_PASS_ONLY.bed ${dir}/WES_V7_80ind_SNP_PASS_ONLY.vcf.gz

wait

echo "done!"
