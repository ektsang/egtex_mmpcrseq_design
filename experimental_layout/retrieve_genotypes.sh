#!/bin/bash

# author: Emily Tsang

set -o nounset -o pipefail

# Subset genotype files to the subset of individuals for which we have eGTEx samples.
# Get V6 OMNI imputed genotypes (all individuals) as well as V7 WES (80/86 inds) and WGS genotypes (79/86 inds).

## The variables in caps on the RHS should be set globally, e.g. in bashrc ##
omni=$GTEX_OMNIv6
wgs=$GTEX_WGSv7
exome=$GTEX_EXOMEv7
dir=${EGTEXDIR}/sampleInfo
#################################################################

indincl=${dir}/participant.info.txt

# output directory/files
outdir=${dir}/genotypes
omniid=${outdir}/omni.ids.txt
wgsid=${outdir}/wgs.ids.txt
exomeid=${outdir}/exome.ids.txt
wgsout=${outdir}/WGS_V7_79ind_SNP_PASS_ONLY.vcf.gz
exomeout=${outdir}/WES_V7_80ind_SNP_PASS_ONLY.vcf.gz
omniout=${outdir}/OMNI_V6_86ind_SNP_PASS_ONLY.vcf.gz

date

# function that makes file with IDs as they appear in the vcf file
getVCFids() {
    vcf=$1
    out=$2
    zcat $vcf | awk '{if(substr($1,1,2)!="##"){print; exit}}' | \
	awk -v indincl=$indincl 'BEGIN{
            while((getline<indincl)>0){
                IDs[$1]
            }
        }{
            for(i=10;i<=NF;i++){
                split($i,IDparts,"-"); 
                indID=IDparts[1]"-"IDparts[2]; 
                if(indID in IDs){print $i}
            }
        }' > $out
}
# actually get ids that are used in the various vcfs
getVCFids $omni $omniid
getVCFids $wgs $wgsid
getVCFids $exome $exomeid

echo "Subsetting VCFs..."

# in all cases, only keep SNPS with at least one observed allele
# the number of individuals with each data type determines the max-non-ref-ac-any (num.ind*2 - 1)
vcftools \
    --gzvcf $wgs \
    --keep $wgsid \
    --non-ref-ac-any 1 \
    --max-non-ref-ac-any 157 \
    --remove-indels \
    --recode \
    --stdout | bgzip > $wgsout &

vcftools \
    --gzvcf $exome \
    --keep $exomeid \
    --non-ref-ac-any 1 \
    --max-non-ref-ac-any 159 \
    --max-missing-count 15 \
    --remove-indels \
    --remove-filtered-all \
    --recode \
    --stdout | bgzip > $exomeout &

# recode-INFO-all to keep the imputation information
vcftools \
    --gzvcf $omni \
    --keep $omniid \
    --non-ref-ac-any 1 \
    --max-non-ref-ac-any 171 \
    --max-missing-count 15 \
    --remove-indels \
    --remove-filtered-all \
    --recode \
    --recode-INFO-all \
    --stdout | bgzip > $omniout &

wait
date

echo "indexing output VCFs..."

tabix -p vcf $wgsout &
tabix -p vcf $exomeout &
tabix -p vcf $omniout &

wait
date
