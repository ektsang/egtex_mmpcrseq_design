#!/bin/bash

set -o nounset

outliers=/srv/scratch/restricted/goats/data/outliers_medz_picked.txt
inds=${EGTEXDIR}/sampleInfo/participant.info.txt
out=${EGTEXDIR}/geneSelection/other/medz.outliers.txt

cat $inds | tail -n +2 | awk '{print $1}' | grep -f - $outliers > $out
