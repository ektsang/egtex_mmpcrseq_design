#!/usr/bin/env python

from __future__ import print_function
import sys

# add adaptor based on designed primers
# for this format, the fwd_adaptor must be the F primer

usage = "usage: python add_adaptor.py inputfile\n"
usage = usage + "inputfile format is pool_name locus_name forward_primer_sequence fp_temp reverse_primer_sequence rp_temp mapping info amp_length annotation"
if len(sys.argv)!= 2:
    print(usage, file = sys.stderr)
    sys.exit(0)

infile = open(sys.argv[1],'r')

fwd_adaptor = "GCGTTATCGAGGTC"
rev_adaptor = "GTGCTCTTCCGATCT"

for line in infile:
    line = line.strip().split()
    poolName = line[0]
    f_primer = fwd_adaptor + line[2]
    r_primer = rev_adaptor + line[4]
    f_name = line[1] + '_F'
    r_name = line[1] + '_R'
    outlines = '\t'.join([poolName, f_name, f_primer]) + '\n' + '\t'.join([poolName, r_name, r_primer])
    print(outlines)

infile.close()
