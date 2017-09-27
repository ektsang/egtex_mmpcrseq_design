#!/usr/bin/env python

# author: Emily Tsang 

# From provided vcf, retrieve the number of hets and other information.

from __future__ import print_function
from __future__ import division
import sys
import argparse
import vcf
import posix

###################################################################################
# CLASSES
###################################################################################

class Variant(object):
    """
    Class to hold information relevant to variant sites.
    """
    def __init__(self, chrom, pos, hetList, imputed, certainty):
        self._chr = "chr" + chrom if len(chrom) < 3 else chrom
        self._pos = pos # int
        self._hetList = hetList # list of strings
        self._imputed = imputed # should be boolean
        self._certainty = certainty # should be float
    
    def __repr__(self):
        return self.prettyString() + "\n"
    
    def __eq__(self, other):
        """
        Checks for equality of two variants based on their chromosomes, positions, and hetLists.
            arg:
                other: object to compare to this variant
            return:
                True if other is a Variant with the same chr, position, and hetList.
        """
        return (isinstance(other, self.__class__) and \
                self._chr == other._chr and \
                self._pos == other._pos and \
                set(self._hetList) == set(other._hetList))
    
    def __hash__(self):
        return (self._chr, self._pos, '.'.join(self._hetList)).__hash__()

    def getChr(self):
        return self._chr

    def getPos(self):
        return self._pos

    def getHets(self):
        return self._hetList

    def nhet(self):
        return len(set(self._hetList))

    def combine(self, other):
        """
        Merges two variants together. Takes the union of their het lists and the minimum of the certainty, if available.
        If one of the two variants is imputed, will mark the combined variant as imputed.
        Only combines the variants if they have the same chr and pos.
            arg:
                other: Variant to combine with this one.
            return:
                Nothing. Just modifies this variant in place.
        """
        if self._chr != other._chr or self._pos != other._pos:
            print("Two variants have different chromosome and/or positions. Not combining.", file = sys.stderr)
            print(self._chr + ":"+ str(self._pos) + " and " + other._chr + ":" + str(other._pos), file = sys.stderr)
        self._hetList = set(self._hetList) | set(other._hetList)
        if self._imputed is not None and other._imputed is not None:
            self._imputed = self._imputed or other._imputed
            self._certainty = min(self._certainty, other._certainty)
    
    def prettyString(self, printImpute = True):
        """
        Print Variant in bed format.
            arg: 
                printImpute: True if want to include imputation and imputation certainty info. 
            return:
                string of the variant's information as a single tab-delimited line.
                Does not include the trailing new line.
        """
        value = "\t".join(map(str, [self._chr, self._pos - 1, self._pos, self.nhet()]))
        value = value + "\t" + ','.join(self._hetList)
        if printImpute:
            if self._imputed is None:
                value = value + "\tNA\tNA"
            else:
                imputedStr = "imp" if self._imputed else "gt"
                value = value + "\t" + imputedStr + "\t" + str(self._certainty)
        return value

###################################################################################
# FUNCTIONS
###################################################################################
def vcfRecordToVariant(record):
    """
    Turns a vcf record into a Variant object.
    Only includes imputed and certainty info if present.
    For multiallelic variants, include the minimum certainty.
        arg:
            record: vcf record from pyvcf
        return:
            Variant object based off that vcf record.
    """
    nhets = record.num_het
    if nhets == 0:
        return None
    imputed = None
    certainty = None
    if 'TYPE' in record.INFO and 'CERTAINTY' in record.INFO:
        imputed = True if record.INFO['TYPE'] == 0 else False
        certainty = min(record.INFO['CERTAINTY'])
    hets = record.get_hets()
    hetList = [h.sample for h in hets]
    # prune names to GTEX-XXXX
    hetList = ['-'.join(het.split('-')[0:2]) for het in hetList]
    variant = Variant(record.CHROM, record.POS, hetList, imputed, certainty)
    return variant

def vcfToVariantList(filename):
    """
    Reads a vcf file into a list of Variant objects.
    Uses pyvcf.
        arg:
            filename: the path the the vcf file
        return:
            List of Variant objects.
    """
    variantList = []
    n = 0
    prevVariant = None
    prevVarkey = "chr0:0"
    with open(filename, 'r') as vcfFile:
        vcfReader = vcf.Reader(vcfFile)
        for vcfRecord in vcfReader:
            variant = vcfRecordToVariant(vcfRecord)
            if variant is not None:
                varkey = variant._chr + ":" + str(variant._pos)
                if varkey == prevVarkey:
                    prevVariant.combine(variant)
                else:
                    variantList.append(prevVariant)
                    prevVariant = variant
                    prevVarkey = varkey
            n += 1
            if n % 10000 == 0:
                print("Parsed " + str(n) + " lines.", file = sys.stderr)
        # deal with last variant and remove first dummy variant
        variantList.append(prevVariant)
        variantList.pop(0)
    return variantList

def bedLineToVariantArgs(line):
    """
    Turn a line in bed format (as output by prettyString) to arguments to initialize a variant.
        arg:
            line: a list of strings from splitting a line in a bedfile
        return:
            A list of the arguments to initialize a variant.
    """
    imputed = None
    certainty = None
    if len(line) == 7:
        imputed = line[5] == "imp"
        certainty = float(line[6])
    variantArgs = [line[0], int(line[2]), line[4].split(','), imputed, certainty]
    return variantArgs

def bedToVariantDict(filename):
    """
    Reads a bed file (generated by this script).
    Transforms it into a variant dict where the key is chr:pos0-pos1.
        arg:
            filename: the path to the bed file
        return:
            Dict of Variant objects with key chr:pos0-pos1.
    """
    variantDict = dict()
    n = 0
    with open(filename, 'r') as bedFile:
        for line in bedFile:
            line = line.strip().split()
            varkey = line[0] + ":" + line[2]
            newvar = Variant(*bedLineToVariantArgs(line))
            n += 1
            if n % 100000 == 0:
                print("Parsed " + str(n) + " lines.", file = sys.stderr)
            if varkey in variantDict:
                variantDict[varkey].combine(newvar)
                continue
            variantDict[varkey] = newvar
    return variantDict

###################################################################################
# MAIN
###################################################################################
def main():
    # parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output', help = 'Output file. If omitted, outputs to stdout.')
    parser.add_argument('vcf', help = 'VCF file to parse.')
    args = parser.parse_args()
    
    # set output file to stdout if none was provided
    outputFile = args.output
    if outputFile is None:
        out = sys.stdout
    else:
        out = open(outputFile, 'w')
    
    print("Parsing vcf into Variant objects", file = sys.stderr)
    # parse vcf into variant object
    variantList = vcfToVariantList(args.vcf)
    # check whether or not to include imputation information based on the first variant
    inclImpute = (variantList[0]._imputed is not None)
    # print them
    for variant in variantList:
        print(variant.prettyString(inclImpute), file = out)

if __name__ == '__main__':
    main()
