#!/usr/bin/env python

## author: Emily Tsang

## Make annotation files for the mmPCR-seq primer design for selected exons and transcripts in order.
## The primer design script requires that the entries in the annotation file be sorted by chromosome,
## but the ordering within a chromosome can be what we wish.
## Only the first transcript that contains the variant of interest will be used, so we must output the
## transcripts in the correct order.

## We make a separate file for the exons and the transcripts with the idea that primers will first
## be designed on the exons. Then sites that failed will attempt primer design on the transcripts.

from __future__ import print_function
from collections import defaultdict
import argparse
import os
import sys
import inspect
cmdFolder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))
sys.path.insert(0, cmdFolder + '/..')
from experimental_layout.retrieve_genes import Exon

###################################################################################
# FUNCTIONS
###################################################################################
def annotation2dict(annotationFile):
    """
    Parse the annotation file into a useful data structure.
        arg: 
            annotationFile: path to the annotation file
        return:
            dict with transcript name as key and the line from the annotation file as value
    """
    annoDict = dict()
    annoFile = open(annotationFile, 'r')
    header = annoFile.readline() # skip header
    assert(header == '#bin\tname\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\t' + \
           'exonStarts\texonEnds\tscore\tname2\tcdsStartStat\tcdsEndStat\texonFrames\n')
    for line in annoFile:
        line = line.strip() # remove new line
        transcript = line.strip().split()[1]
        chrom = line.strip().split()[2]
        assert(transcript not in annoDict or chrom == "chrX" or chrom == "chrY")
        annoDict[transcript] = line
    annoFile.close()
    return annoDict

def sites2dict(siteFile):
    """
    Parse the selected sites file into different chromosomes.
        arg:
            siteFile: path to the sites file
        return:
            dict with chromosome as key and (ordered) list of Exon objects as value
    """
    exDict = defaultdict(list)
    infile = open(siteFile, 'r')
    header = infile.readline() # skip header
    assert(header == 'varname\tHGNC\tENSG\tstrand\tvar.chr\tvar.pos0\tvar.pos1\tnhet\tgt\tcertainty\tVEP' + \
           '\tintra.exon.plausible\texon.chr\texon.pos0\texon.pos1\toverlap.same\toverlap.other\ttype\t' + \
           'transcript\trank\tmedian.TPM\tmad.TPM\tnon.zero.TPM\tmedian.perc\tmad.perc\tnon.zero.perc\n')
    for line in infile:
        line = line.strip().split()
        chrom = line[12]
        exStart = int(line[13])
        exEnd = int(line[14])
        transcriptPlus = (line[18], line[1], line[3])
        # IMPORTANT NOTE
        # keeping bed indexing, which is not the intended use of exon
        # also using a single transcript, gene, strand tuple instead of a list of transcript, transcript type pairs
        # won't be able to use all the exon functions, but will work for our purposes here
        newExon = Exon('', chrom, exStart, exEnd, transcriptPlus)
        if newExon not in exDict[chrom]:
            exDict[chrom].append(newExon)
    infile.close()
    return exDict

def outputSubsettedAnnotations(outPrefix, exDict, annoDict):
    """
    Produce correctly formatted annotation files with exon and transcript information.
    Process and output chromosomes in order.
    Within each chromosome, process the exons in their listed order.
        args:
            outPrefix: prefix out output file paths
            exDict: dict of exons keyed by chromosome
            annoDict: dict of annotation keyed by transcript name
        return:
            Nothing. Outputs results to two files. One with exons and one with transcripts.
            The file name suffixes are '_selected_exons' and '_selected_transcripts', respectively.
    """
    exonFile = open(outPrefix + '_selected_exons', 'w')
    transcriptFile = open(outPrefix + '_selected_transcripts', 'w')
    header = '#bin\tname\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\t' + \
             'exonStarts\texonEnds\tscore\tname2\tcdsStartStat\tcdsEndStat\texonFrames'
    # write headers
    print(header, file = exonFile)
    print(header, file = transcriptFile)
    
    transcriptsSeen = set() # to keep track of transcripts already processed
    sortedChroms = sorted(exDict.keys(), key = lambda x: int(x[3:]))
    for chrom in sortedChroms:
        for exon in exDict[chrom]:
            # process exon
            transcript, gene, strand = exon.getTranscripts()
            exChrom, exStart, exEnd = exon.getCoordinates()
            exonLine = '\t'.join(['.', transcript, chrom, strand, str(exStart), str(exEnd),
                                  '.', '.', '1', str(exStart), str(exEnd), '.', gene + '_exon', '.', '.', '.'])
            print(exonLine, file = exonFile)
            # process transcript
            if transcript in transcriptsSeen:
                continue
            if transcript not in annoDict:
                print("Transcript " + transcript + " not in annotation file to subset. Skipping.", file = sys.stderr)
                transcriptsSeen.add(transcript)
            else:
                print(annoDict[transcript], file = transcriptFile)
                transcriptsSeen.add(transcript)
    exonFile.close()
    transcriptFile.close()

###################################################################################
# MAIN
###################################################################################
def main():
    # parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('annotation', help = 'Annotation file in correct format for mmPCR-seq primer design.')
    parser.add_argument('sites', help = 'File in format output by experimental_layout/get_sites_in_exons.py')
    parser.add_argument('outputPrefix', help = 'Prefix of output files.')
    args = parser.parse_args()
    
    # read annotation into a dict keyed by transcript name with annotation line as value
    annotationDict = annotation2dict(args.annotation)

    # read sites into dict keyed by chromosome with a list of Exons as value
    exonDict = sites2dict(args.sites)

    # iterate through the sites and output transcripts and exons to file
    outputSubsettedAnnotations(args.outputPrefix, exonDict, annotationDict)

if __name__ == '__main__':
    main()
