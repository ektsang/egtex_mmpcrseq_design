#!/bin/bash

# author: Emily Tsang

# first get fasta from gene annotation
# then use fasta to make a blast database

set -o nounset -o errexit -o pipefail

## The variables in caps on the RHS should be set globally, e.g. in bashrc ##
dir=$BLASTDIR
genomedir=$GENOMEDIR
#################################################################

scriptdir=`dirname \$(readlink -f "\$0")`


bedtools getfasta -split -name -s -fi ${genomedir}/hg19.fa -bed ${scriptdir}/UCSC_TableBrowser_gencode19.bed -fo ${dir}/hg19.gencode19.fa

# remove duplicate sequences (some seqIDs are on both X and Y chromosomes with same name)
# uses a bit of a hack where I first make each ID, sequence pair be on one line
# then use that format to only keep unique pairs
# finally return it to the original format (with the ID on one line and the sequence on the next)
cat ${dir}/hg19.gencode19.fa | tr '\n' '\t' | sed 's/\t>/\n>/g' | sort | uniq -i | tr '\t' '\n' > tmp
mv tmp ${dir}/hg19.gencode19.fa

makeblastdb -in ${dir}/hg19.gencode19.fa -parse_seqids -dbtype nucl -title hg19.gencode19.blast -out ${dir}/hg19.gencode19.blast.db -logfile ${dir}/hg19.gencode19.blastdb.log

