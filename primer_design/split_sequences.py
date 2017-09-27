#!/usr/bin/env python

# author: Emily Tsang

# split a fasta file into multiple files
# create one file for each sequence (chromosome, in the case of a genome fasta)

from __future__ import print_function
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument('inputfile', help = 'input fasta')
parser.add_argument('outputdir', help = 'output directory')
args = parser.parse_args()

fasta = open(args.inputfile, 'r')

subfasta = None

for line in fasta:
    if line[0] == ">":
        # for all but the first sequence
        if subfasta is not None:
            subfasta.close()
        seqname = line.strip()[1:]
        outfile = args.outputdir + "/" + seqname + ".fa"
        subfasta = open(outfile, 'w')
    print(line, file = subfasta, end = "")
# close last sequence file
subfasta.close()
# close input fasta file
fasta.close()
