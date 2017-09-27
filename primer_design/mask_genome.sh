#!/bin/bash

# author: Emily Tsang

## Make the genome one file per chromosome and mask variant bases in these new fasta files
## Variants common (MAF > 0.05) in 1KG or at any frequency in any of the eGTEx individuals are masked

set -o nounset -o errexit -o pipefail

## The variables in caps on the RHS should be set globally, e.g. in bashrc ##
genomedir=$GENOMEDIR
kgdir=$KGDIR
wgs=$GTEX_WGSv7
omni=$GTEX_OMNIv6
wgsids=${EGTEXDIR}/sampleInfo/genotypes/wgs.ids.txt
omniids=${EGTEXDIR}/sampleInfo/genotypes/omni.ids.txt
#################################################################

outdir=${genomedir}/egtex_masked

# create output directory if it doesn't exits
if [ ! -d $outdir ]; then
    mkdir $outdir
fi

# split genome into individual chromosomes
python split_sequences.py ${genomedir}/hg19.fa ${genomedir}

# remove non-canonical chromosomes
rm ${genomedir}/*_*.fa


# get list of sites to mask
# function to process vcftools output
freq2bed() {
    awk 'NR > 1{
                split($5,splitcol,":");
                alleleLength = length(splitcol[1]);
                print "chr"$1"\t"$2-1"\t"$2-1+alleleLength
               }'
}
export -f freq2bed

# get list of common sites in 1KG in bed format
# remove output file if it exists
if [ -f ${outdir}/KGvars.bed ]; then
    rm -f ${outdir}/KGvars.bed
fi
# then populate it
for chrom in chr{1..22}
do
    vcftools --gzvcf ${kgdir}/ALL.${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --stdout \
	     --remove-filtered-all --maf 0.05 --freq | \
	freq2bed >> ${outdir}/KGvars.bed
done

# get list of sites in eGTEx individuals (using both WGS and genotyping data)
# WGS: 79 individuals
vcftools --gzvcf $wgs --keep $wgsids --non-ref-ac-any 1 --max-non-ref-ac-any 157 --stdout --freq | \
    freq2bed >> ${outdir}/eGTEX_WGSvars.bed
# OMNI: 86 individuals
vcftools --gzvcf $omni --keep $omniids --non-ref-ac-any 1 --max-non-ref-ac-any 171 \
	 --max-missing-count 15 --remove-filtered-all --stdout --freq | \
    freq2bed >> ${outdir}/eGTEX_OMNIvars.bed

# merge the variant lists
cat ${outdir}/KGvars.bed ${outdir}/eGTEX_WGSvars.bed ${outdir}/eGTEX_OMNIvars.bed | \
    sort -k1,1 -k2,2n | \
    bedtools merge -i stdin > ${outdir}/masking_vars.bed

# mask fasta with joint list
for chrom in chr{1..22}
do
    grep -F $chrom ${outdir}/masking_vars.bed | \
	bedtools maskfasta -fi ${genomedir}/${chrom}.fa -bed stdin -fo ${outdir}/${chrom}.fa
done

