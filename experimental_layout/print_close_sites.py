#!/usr/bin/env python

# author: Emily Tsang

## takes as input sorted sites and prints out ones that may have overlapping amplicons
## assumes the following input format:
## site.name hgnc ensg strand chr pos0 [other columns not used]

from __future__ import print_function
import fileinput

prevChrom = ''
prevPos = -1
prevStrand = ''
prevLine = ''

for line in fileinput.input():
    line = line.strip()
    lineSplit = line.split()
    strand = lineSplit[3]
    chrom = lineSplit[4]
    pos = int(lineSplit[5])
    if chrom != prevChrom:
        prevChrom = chrom
        prevPos = pos
        prevStrand = strand
        prevLine = line
        continue
    # if criteria for closeness are met, print the pair of sites
    assert(pos >= prevPos)
    if (pos - prevPos <= 766 and prevStrand == "+" and strand == "-") or \
       (pos - prevPos <= 453 and prevStrand == strand) or \
       (pos - prevPos <= 140 and prevStrand == "-" and strand == "+"):
        print(prevLine)
        print(line)
        print("^ " + str(pos - prevPos) + " -------------------------------")
    # reset previous position
    prevChrom = chrom
    prevPos = pos
    prevStrand = strand
    prevLine = line

print("DONE!")

    
