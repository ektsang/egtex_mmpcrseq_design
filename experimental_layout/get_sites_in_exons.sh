#!/bin/bash

# author: Emily Tsang

set -o nounset -o errexit -o pipefail

## The variables in caps on the RHS should be set globally, e.g. in bashrc ##
dir=${EGTEXDIR}
vep=${GTEX_WGSv7vep}
bias=${BIAS}
mappability=${MAPPABILITY}
################################################################# 

outdir=${dir}/selectedSites
gtdir=${dir}/sampleInfo/genotypes
exons=${dir}/geneSelection/selected.genes.exons.bed
omni=${gtdir}/OMNI_V6_86ind_SNP_PASS_ONLY.bed
wgs=${gtdir}/WGS_V7_79ind_SNP_PASS_ONLY.bed
transcripts=${outdir}/transcripts.txt

passed=${outdir}/mappable.sites.txt
omnifilt=${gtdir}/OMNI_V6_86ind_SNP_PASS_ONLY_MAPPABLE.bed
wgsfilt=${gtdir}/WGS_V7_79ind_SNP_PASS_ONLY_MAPPABLE.bed

scriptdir=`dirname \$(readlink -f "\$0")`

## first remove potentially biased sites from the bed files
## this means keep only sites with a mappability of 1 and without demonstrated bias in simulations
cat $bias | awk 'NR>1{print "chr"$2"\t"$3-1"\t"$3}' | \
    bedtools subtract -a $mappability -b stdin | \
    awk '$4==1{print $1"\t"$2"\t"$3}' > $passed

bedtools intersect -wa -a $omni -b $passed > ${omnifilt}
bedtools intersect -wa -a $wgs -b $passed > ${wgsfilt}

## then run the intersection with the genes
python ${scriptdir}/get_sites_in_exons.py --extravar $omnifilt --vep $vep -o ${outdir}/genes.intersect.WGS79.OMNI86 $wgsfilt $exons $transcripts


