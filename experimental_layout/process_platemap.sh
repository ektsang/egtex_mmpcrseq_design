#!/bin/bash

# author: Emily Tsang

set -o nounset -o errexit -o pipefail

# break the master plate map into useful reference files

## The variables in caps on the RHS should be set globally, e.g. in bashrc ##
dir=${EGTEXDIR}/sampleInfo
wgs=$GTEX_WGSv7
wes=$GTEX_EXOMEv7
omni=$GTEX_OMNIv6
rnaseq=${GTEX_RNAv7}/GTEx_Analysis_2016-01-15_v7_RSEMv1.2.22_transcript_tpm.txt.gz
#################################################################

map=${dir}/master.platemap.txt

cd $dir

# get list of individuals with WGS and WES in v7 and array data in v6
zcat $wgs | head -n 300 | grep "#CHROM" | sed 's/\t/\n/g' | grep "GTEX" | \
    tr '-' ' ' | awk '{print $1"-"$2}' > wgs.samples

zcat $wes | head -n 300 | grep "#CHROM" | sed 's/\t/\n/g' | grep "GTEX" | \
    tr '-' ' ' | awk '{print $1"-"$2}' > exome.samples

zcat $omni | head -n 300 | grep "#CHROM" | sed 's/\t/\n/g' | grep "GTEX" | \
    tr '-' ' ' | awk '{print $1"-"$2}' > omni.samples

# get participant IDs and the number of tissues each
cat $map | tail -n +2 | cut -f12 | sort | uniq -c | \
    awk 'BEGIN{OFS="\t"; 
           while((getline < "wgs.samples") > 0){
             wgs[$1];
           };
           while((getline < "exome.samples") > 0){
             exome[$1];
           };
           while((getline < "omni.samples") > 0){ 
             omni[$1];
           };
           print "ParticipantID\tTissueCount\thas.WGS.V7\thas.WES.V7\thas.OMNI.V6"
         }{print $2,$1,$2 in wgs,$2 in exome,$2 in omni}' | \
    column -t > participant.info.txt

rm wgs.samples exome.samples omni.samples

# get tissue names and the number of samples each
cat $map | tail -n +2 | cut -f16 | \
    sed 's/ - /_/g' | sed 's/ /_/g' | sed 's/(//g' | sed 's/)//g' | \
    sort | uniq -c | \
    awk 'BEGIN{OFS="\t"; print "Tissue\tnSamples"}{print $2,$1}' | \
    column -t > tissue.sample.counts.txt

# get list of samples with RNA-seq data in v7
zcat $rnaseq | head -n 1 | sed 's/\t/\n/g' | grep "GTEX" | \
    tr '-' ' ' | awk '{if(NF==5){print $1"-"$2"-"$3} else {print $1"-"$2"-"$3"-"$4}}' > rnaseq.samples 

# get list of samples we received that already have rna-seq data in v7
cat $map | tail -n +2 | cut -f13 | sort | \
    awk 'BEGIN{OFS="\t";
           while((getline < "rnaseq.samples") > 0){
             rna[$1];
           };
           print "Sample\thasRNAseq.V7"
         }{print $1,$1 in rna}' | \
    column -t > sample.hasRNAseq.v7.txt

rm rnaseq.samples
