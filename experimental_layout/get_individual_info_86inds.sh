#!/bin/bash

# author: Emily Tsang 

set -o nounset -o errexit -o pipefail

# get subject information for the 86 individuals included in the eGTEx project

## The variables in caps on the RHS should be set globally, e.g. in bashrc ##
dir=${EGTEXDIR}/sampleInfo
anno=$GTEX_SUBJECTSv7
#################################################################

output=${dir}/subject.phenotypes.v7.86inds.txt

# get header of annotation file
head -n 1 $anno > $output

# then append relevant individuals
grep -f <(cat ${dir}/participant.info.txt | tail -n +2 | awk '{print $1}') $anno >> $output
