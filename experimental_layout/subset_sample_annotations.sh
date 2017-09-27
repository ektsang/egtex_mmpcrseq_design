#!/bin/bash

# author: Emily Tsang

set -o nounset -o errexit -o pipefail

# get subject, sample, tissue, and whether or not to use the sample from the sample annotation file
# retrieve this for all RNAseq samples available for the given set of individuals

## The variables in caps on the RHS should be set globally, e.g. in bashrc ##
dir=${EGTEXDIR}/sampleInfo
anno=$GTEX_SAMPLESv7
#################################################################

subjects=${dir}/participant.info.txt
out=${dir}/subject.sample.mapping.v7.txt

echo -e "subject\tsample\ttissue\tstatus" > $out

tail -n +2 $subjects | awk '{print $1}' | grep -f - $anno | cut -f1,14,26,27,28 | \
    awk 'BEGIN{FS=OFS="\t"}
         {split($1,sample,"-");
          if($3!="") {
            print sample[1]"-"sample[2],$1,$2,$3":"$4
          } else if ($5!=""){
            print sample[1]"-"sample[2],$1,$2,$5
          } else {
            print sample[1]"-"sample[2],$1,$2,"NOT_USED"
          }}' | \
    sed 's/ - /_/g' | sed 's/ /_/g' | sed 's/(//g' | sed 's/)//g' \
    >> $out
