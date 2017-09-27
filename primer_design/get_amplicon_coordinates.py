#!/usr/bin/env python

# author: Emily Tsang

# Take as input the output from the primer design, the gene annotation files,
# and the site information. Use it to infer the amplicon coordinates in genomic space.
# Output the sites and their amplicon coordinates.

# The annotation files are the ones used as input to the mmpcr primer design script.

# The output file is a bed-like format with columns:
# chr amplicon_start(0-based) amplicon_end(1-based) site_name transcript site_pos(1-based) single_exon

from __future__ import print_function
import argparse
import sys
import numpy as np

###################################################################################
# Classes
###################################################################################

class Annotation(object):
    def __init__(self, chrom, start, end, strand, sizeList, startList):
        self._chrom = chrom
        self._chromStart = int(start)
        self._chromEnd = int(end)
        self._strand = strand
        self._exonSizeList = map(int, sizeList)
        self._exonStartList = map(int, startList)
        assert(self._strand in ('+','-'))
        assert(len(self._exonSizeList) == len(self._exonStartList))
    
    def nexons(self):
        return len(self._exonSizeList)
        
    def getExonIndexAndOffset(self, index, pos, length):
        '''
        Helper function for getGenomicPositions.
        Gets the offset within the specified exon of the given position.
        Also transforms the index if the transcript is on the reverse strand.
        All input values are relative to the TSS, so will be from the "end"
        for transcripts on the reverse strand. The output values should be
        relative to the "start" (i.e., smaller coordinates).
            args:
                index: the number of the given exon
                pos: relative position for which we want the exon offset
                length: the length of the transcript up to and including the given exon
            return:
                tuple of the exon index (counting from the start of the list - small to large coordinates)
                and the offset within the exon (again from the smaller coordinate end)
        '''
        if self._strand == '-':
            index = self.nexons() - index - 1
            offset = length - pos + 1
        else:
            offset = self._exonSizeList[index] - (length - pos)
        return index, offset
        
    def getGenomicPositions(self, name, start, end, genomicPos):
        '''
        Translate the positions, which are relative to the annotation feature
        (transcript or exon), into genomic coordinates.
            args:
                name: site name, for error logging purposes
                start: the start coordinate relative to the feature this annotation object refers to
                end: relative end coordinate
                genomicPos: the position targeted, in genomic coordinates, as (chr, pos) tuple. 
                            it should fall within the start and end boundaries of the feature.
            return:
                tuple of the start and end coordinates in genomic space, and whether they are in the same exon
        '''
        # check that genomic position has same chromosome as transcript
        if genomicPos[0] != self._chrom:
            print(name, 'chromosome', genomicPos[0], 'does not match transcript chomosome',
                  self._chrom, file = sys.stderr)
        # start counting from back if on reverse strand
        if self._strand == '-':
            cumLengths = np.cumsum(self._exonSizeList[::-1])
        else:
            cumLengths = np.cumsum(self._exonSizeList)

        ## find the exon and the position in that exon of the provided start and end
        startFound = False
        endFound = False
        for i, length in enumerate(cumLengths):
            if not startFound and start <= length:
                startFound = True
                startExon, startOffset = self.getExonIndexAndOffset(i, start, length)
            if startFound and end <= length:
                endExon, endOffset = self.getExonIndexAndOffset(i, end, length)
                endFound = True
                break
        
        ## At least one amplicon had the end boundary larger than the transcript size.
        ## I checked on UCSC genome browser and the primer spans past the transcript by a few bases.
        ## It's unclear why that happened, but just setting the end position to be the transcript length in that case.
        firstCoord = self._chromStart + self._exonStartList[startExon] + startOffset
        if endFound:
            secondCoord = self._chromStart + self._exonStartList[endExon] + endOffset
        else:
            print(name, ': end position', end, 'outside transcript, which has length',
                  cumLengths[-1], file = sys.stderr)
            secondCoord = self._chromEnd - 1 if self._strand == '+' else self._chromStart
            endExon = self.nexons() - 1 if self._strand == '+' else 0

        ## flip the start and end positions if the transcript is on the reverse strand
        ## also add 1 to the end coordinate for bed formatting
        if self._strand == '+':
            genomicStart = firstCoord
            genomicEnd = secondCoord + 1
        else:
            genomicStart = secondCoord
            genomicEnd = firstCoord + 1

        # check that the genomic position provided falls between the start and the end
        if genomicPos[1] <= genomicStart or genomicPos [1]> genomicEnd:
            print(name, ': Given position', genomicPos, 'not between genomic start and end positions',
                  genomicStart, '-', genomicEnd, file = sys.stderr)
        
        sameExon = (startExon == endExon)
        return genomicStart, genomicEnd, sameExon
            

###################################################################################
# Functions
###################################################################################

def parseSites(filename):
    '''
    Create a dict mapping site names to positions.
        arg:
            filename: input file with site information in the format 
                      output by experimental_layout/get_sites_in_exons.py
        return:
            dict with the site name as key and 
            a tuple of the chromosome and the 1-based coordinate as value
    '''
    f = open(filename, 'r')
    site2pos = dict()
    # sanity checks - making sure the file has the expected format
    header = f.readline().strip().split()
    assert(header[0] == 'varname')
    assert(header[4] == 'var.chr')
    assert(header[6] == 'var.pos1')
    for line in f:
        line = line.strip().split()
        site2pos[line[0]] = (line[4], int(line[6]))
    f.close()
    return site2pos

def parseAnnotations(filename):
    '''
    Read through file and build a dict with transcript names as keys 
    and Annotation objects as values.
        arg:
            filename: File with the annotation information
        return:
            dict of Annotation objects keyed by transcript name
    '''
    transDict = dict()
    f = open(filename, 'r')
    for line in f:
        line = line.strip().split()
        chrom, start, end , name, dummy, strand = line[0:6]
        sizeList = line[10].strip(',').split(',')
        startList = line[11].strip(',').split(',')
        transDict[name] = Annotation(chrom, start, end, strand, sizeList, startList)
    f.close()
    return transDict

def extractAmpliconCoordinates(filename, site2pos, annoDict):
    f = open(filename, 'r')
    for line in f:
        line = line.strip().split()
        siteName = line[1]
        amplicon = line[6]
        
        ampliconInfoList = amplicon.split('_')[1].split(':')
        transcript = ampliconInfoList[0]
        relativeStart, relativeEnd = map(int, ampliconInfoList[1].split('-'))
        pos = site2pos[siteName]
        genomicStart, genomicEnd, sameExon = annoDict[transcript].getGenomicPositions(siteName, relativeStart, relativeEnd, pos)
        print('\t'.join(map(str, [pos[0], genomicStart, genomicEnd, siteName, transcript, pos[1], sameExon])))
    f.close()
    
###################################################################################
# Main
###################################################################################

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('sitefile', help = 'File with site information')
    parser.add_argument('primerfile', help = 'File with the designed primers.')
    parser.add_argument('annotationfile', help = 'File with gene annotation used for blast database (bed format).')
    args = parser.parse_args()

    # create the mapping of site names to positions
    site2posDict = parseSites(args.sitefile)

    # parse the annotation files into a dict keyed by transcript names
    annotationDict = parseAnnotations(args.annotationfile)

    # run through the primer file and output desired info one while doing so
    extractAmpliconCoordinates(args.primerfile, site2posDict, annotationDict)

if __name__ == '__main__':
    main()
