#!/usr/bin/env python

# author: Emily Tsang

# takes as input gene list (gene name in first column, 0/1 exclusion/inclusion in second column)
# and gene annotation (Gencode GTF file)

from __future__ import print_function
from __future__ import division
import sys
import argparse
import pybedtools

###################################################################################
# CLASSES
###################################################################################

class Exon(object):
    """
    Exon definition.
    Coordinates work on a 1-based system and both end points are included.
    Note that transcript is expected to be a list of Tuples (transcript name, transcript type)
    """
    def __init__(self, id, chrom, start, end, transcript, overlapSameString = '', overlapDiffString = ''):
        self._name = id
        self._chr = "chr" + chrom if len(chrom) < 3 else chrom
        self._start = start
        self._end = end
        self._transcriptList = transcript
        if overlapSameString == '':
            self._overlapSameList = []
        else:
            self._overlapSameList = overlapSameString.split(";") # list of chr:start-end (keys in gene dict)
        if overlapDiffString == '':
            self._overlapDiffList = []
        else:
            self._overlapDiffList = overlapDiffString.split(";")
        self._variantList = []
        
    def __eq__(self, other):
        """
        Checks for equality of two exons based on their chr, start, and end positions.
            arg:
                other: object to compare to this exon
            return:
                True if other is an exon with the same chr, start, and end.
        """
        return (isinstance(other, self.__class__) and \
            self._chr == other._chr and \
            self._start == other._start and \
            self._end == other._end)
    
    def __ne__(self, other):
        return not self.__eq__(other)
    
    def __hash__(self):
        """ Returns a hash on the string composed of chr, start, and stop, so exons that are equal will have the same hash."""
        return hash("-".join([self._chr, self._start, self._end]))
    
    def __repr__(self):
        """ Returns a verbose representation of the exon."""
        value = "Exon " + str(self._name) + " => " 
        value = value + str(self._chr) + ":" + str(self._start) + "-" + str(self._end) + "\n"
        value = value + "Overlaps other exons in this gene " + str(self.numOverlapsSame()) + " times\n"
        value = value + "These exons are: " + str(self._overlapSameList) + "\n"
        value = value + "Overlaps other exons in other genes " + str(self.numOverlapsDiff()) + " times\n"
        value = value + "These exons are: " + str(self._overlapDiffList) + "\n"
        value = value + "Belongs to the " + str(self.numTranscripts()) + " following transcripts: "
        value = value + str(self._transcriptList) + "\n"
        value = value + "Variants:\n==\n"
        for variant in self._variantList:
            value = value + repr(variant)
        value = value + "==\n"
        return value
    
    def __str__(self):
        """ Returns a one-line representation **in bed format** of the exon (without headers), for easy printing."""
        valueList = [self._chr, self._start-1, self._end, self._name, \
                     ';'.join(map(lambda x: "(" + x[0] + "," + x[1] + ")", self._transcriptList)), \
                     ";".join(self._overlapSameList), ";".join(self._overlapDiffList)]
        return "\t".join(map(str, valueList))

    def getCoordinates(self):
        return (self._chr, self._start, self._end)
    
    def length(self):
        return (self._end - self._start + 1)
    
    def incrementOverlapSameGene(self, exon):
        """
        Increments overlap count of exon with other exons in same gene.
            arg: 
                exon: The Exon (object) with which the overlap was found.
            return:
                Nothing. Updates the exon directly.
        """
        self._overlapSameList.append(str(exon._chr) + ":" + str(exon._start) + "-" + str(exon._end))
    
    def incrementOverlapDiffGene(self, otherChr, otherStart, otherEnd):
        """ 
        Adds info of overlapping exon from another gene.
            args:
                otherChr: chromosome of the overlapping exon
                otherStart: start coordinate (0-indexed) of the overlapping exon as a string
                otherEnd: end coordinate (1-indexed) of the overlapping exon as a string
            return:
                True if the exon needs to be removed.
        """
        self._overlapDiffList.append(otherChr + ":" + otherStart + "-" + otherEnd)
        
    def numTranscripts(self):
        """ Returns the number of transcripts that contain this exon."""
        return len(self._transcriptList)

    def numOverlapsSame(self):
        """ Returns the number of exons in the same gene that this exon overlaps."""
        return len(self._overlapSameList)

    def numOverlapsDiff(self):
        """ Returns the number of exons in other genes that this exon overlaps."""
        return len(self._overlapDiffList)

    def addTranscript(self, transcript):
        """
        Adds a transcript to the transcript list.
            arg:
                transcript: the transcript to add to the exon
            return:
                Nothing. Modifies the exon's transcript list.
        """
        if transcript not in self._transcriptList:
            self._transcriptList.extend(transcript)
        else:
            print("Transcript " + transcript + " was already present. Not adding it.", file = sys.stderr)

    def getTranscripts(self):
        return self._transcriptList

    def getTranscriptNames(self):
        """Return a list of the transcript ids."""
        return [t[0] for t in self._transcriptList]

    def removeTranscripts(self, transcriptNameList):
        """Remove transcripts with names in the provided list."""
        self._transcriptList = [(tname, ttype) for tname, ttype in self._transcriptList if tname not in transcriptNameList]
    
    def addVariant(self, variant):
        """
        Adds a variant to the variant list.
            args:
                variant: Variant object to add to the exon
            return:
                Nothing. Modifies the exon's variant list.
        """
        self._variantList.append(variant)

    def nvariants(self):
        return len(self._variantList)
            
    def getVariants(self):
        """Return the Exon's variant list."""
        return self._variantList

    def getOverlapSameList(self):
        return self._overlapSameList
    
    def getOverlapDiffList(self):
        return self._overlapDiffList

    def removeVariant(self, variant):
        """Remove the given variant object. Must be in the list."""
        self._variantList.remove(variant)
            
    def overlaps(self, otherExon, minprop = 0):
        """
        Function to determine whether another exon overlaps with the current one.
        Function is symmetric so A.overlaps(B) == B.overlaps(A).
        Overlap proportion is with respect to the shorter exon
        args:
            otherExon: another exon object to compare to
            minprop: minimum proportion of length of shorter exon that must overlap
        return:
            True if at least prop of the shorter exon is overlapped.
        """
        # input check
        if minprop < 0 or minprop > 1:
            print("Invalid proportion. Must be between 0 and 1. Returning False.", file = sys.stderr)
            return False
        # check for equality
        if self == otherExon:
            return True
        # check chromosome       
        if otherExon is None or self._chr != otherExon._chr:
            return False
        if self._start > otherExon._end or self._end < otherExon._start:
            return False
        # overlap: check what proportion overlaps
        length = min(self.length(),otherExon.length())
        # ------------ OR ------
        #     ------        ---------
        if self._start <= otherExon._start:
            if self._end >= otherExon._end:
                return otherExon.length()/length >= minprop
            else: 
                return (self._end - otherExon._start)/length >= minprop
        #    ------------   OR   ---------
        # -----------------    ---------
        if self._start > otherExon._start:
            if self._end <= otherExon._end:
                return True
            else:
                return (otherExon._end - self._start)/length >= minprop


class Gene(object):
    """
    A gene is defined as a set of exons.
    This is stored as a dict where the key is "chr:start-end" and the value is an exon object.
    """
    def __init__(self, name, ensg, strand, geneType = 'Unspecified'):
        self._name = name
        self._ensg = ensg
        self._type = geneType
        self._strand = strand
        self._exons = dict()
    
    def __repr__(self):
        """
        Verbose gene representation which mostly relies on the verbose exon representations.
        """
        value = "Gene name: " + self._name
        value = value + "\nENSG: " + self._ensg
        value = value + "\ngene type: " + self._type
        value = value + "\nstrand: " + self._strand
        value = value + "\nExons:\n=====\n"
        for exon in self.getSortedExons():
            value = value + repr(exon)
        value = value + "=====\n"
        return value
    
    def __str__(self):
        """
        Returns string for printing as you might want it in a file, 
        one line per exon, with the gene information on every line.
        The exon information is printed first, in bed format.
        """
        valueList = [self._name, self._ensg, self._strand, self._type]
        geneString = '\t'.join(valueList)
        value = ''
        for exon in self.getSortedExons():
            value = value + str(exon) + '\t' + geneString + '\n'
        return value
    
    def numExons(self):
        """ Return the number of exons this gene contains."""
        return len(self._exons)

    def getName(self):
        return self._name

    def getENSG(self):
        return self._ensg

    def getStrand(self):
        return self._strand
    
    def getTranscriptNames(self):
        """ Return the set of transcript names associated with any of the gene's exons as a list."""
        names = []
        for exon in self._exons.values():
            names.extend(exon.getTranscriptNames())
        return list(set(names))
    
    def getVariants(self):
        """ Return Variants associated with any of the gene's exons as a single Variant list."""
        variantList = []
        for exon in self._exons.values():
            variantList.extend(exon.getVariants())
        return list(set(variantList))

    def getExons(self):
        return self._exons

    def getSortedExons(self):
        """ Returns list of exons in sorted order by exon name."""
        return sorted(self._exons.values(), key = lambda x: int(x._name.split("_")[1]))

    def removeExon(self, exonkey):
        """Remove an exon from the gene's exonDict based on the exon's key in the dict."""
        del self._exons[exonkey]
    
    def addExon(self, chrom, start, end, transcript):
        """
        Add an exon to the gene. 
        Checks that the exon doesn't already exist.
        If the exon already exists, adds transcript information to existing exon.
        Also checks for overlaps with other exons in the gene and updates that information if necessary.
            args:
                chrom: the chromosome of the exon
                start: the start position of the exon 
                end: the end position of the exon
                transcript: the name of the transcript this exon is associated with
            return:
                Nothing. Modifies the gene and or the gene's exons in place.
        """
        # check if the exon already exists in the gene
        # if so, modify it
        # if not, create it
        exonKey = str(chrom) + ":" + str(start) + "-" + str(end)
        if exonKey in self._exons:
            self._exons[exonKey].addTranscript(transcript)
        else:
            newExonID = "Exon_" + str(self.numExons() + 1)
            newExon = Exon(newExonID, chrom, start, end, transcript)
            # mark substantial overlap (>80%) with other exons already in the gene
            for prevExon in self._exons.values():
                if newExon.overlaps(prevExon):
                    newExon.incrementOverlapSameGene(prevExon)
                    prevExon.incrementOverlapSameGene(newExon)
            self._exons[exonKey] = newExon
    
    def addExonComplete(self, variant, chrom, start, end, name, transcriptsString, overlapSameString, overlapDiffString):
        """
        Build exon object from all the information provided as arguments and add it to the gene.
        All arguments can be provided as strings.
            args:
                variant: Variant object that overlaps with this exon.
                chrom: chromosome
                start: exon start, 0-based and included, will be converted to 1-based
                end: exon end, 1-based and included
                transcriptString: transcript list in string form, as output by str of Exon
                overlapSameString: list of overlapping exon coordinates in string form (semi-colon-delimited)
                overlapDiff: count of exons in different genes that this exon overlaps
            return: Nothing. Creates Exon and adds it to the gene, if it isn't already there.
                If the Exon exists, simply adds the variant to it.
        """
        # check if the exon already exists in the gene
        # if not, create it
        exonKey = str(chrom) + ":" + str(int(start)+1) + "-" + str(end)
        if exonKey not in self._exons:
            # turn transcript string into transcript list (of tuples)
            transcripts = [tuple(t.rstrip(')').lstrip('(').split(',')) for t in transcriptsString.split(";")]
            self._exons[exonKey] = Exon(name, chrom, int(start)+1, int(end), transcripts, overlapSameString, overlapDiffString)
        self._exons[exonKey].addVariant(variant)

###################################################################################
# FUNCTIONS
###################################################################################
def processGtfLine(geneDict, line):
    """
    Parse a line from a gtf and put the information into the gene dict.
    If the gene already exists, add/adjust exons of that gene.
    Otherwise, make the gene and add the exon.
        args:
            geneDict: dictionary where the keys are gene names and the values are Gene objects
            line: line from a gtf (as pybedtools Interval)
        return:
            Nothing. Modifies the geneDict object.
    """
    if line[2] != "exon":
        print("Warning: skipping the following line does not refer to an exon:\n" + str(line), file = sys.stderr)
    # from annotations, get ENSG, gene name, gene type, transcript ID, and transcript type
    annoDict = line.attrs # pybedtools automatically parses the attributes into a dict
    geneName = annoDict['gene_name']
    ensg = annoDict['gene_id']
    geneType = annoDict['gene_type']
    transcriptID = annoDict['transcript_id']
    transcriptType = annoDict['transcript_type']
    # if gene is not in the dict, add it
    if geneName not in geneDict:
        newGene = Gene(geneName, ensg, line.strand, geneType)
        geneDict[geneName] = newGene
    # update exon info
    geneDict[geneName].addExon(line.chrom, line.start+1, line.stop, [(transcriptID,transcriptType)]) # forcing 1-based indexing

def processIntersection(geneDict, line):
    """
    Updates exon information to include cases where it overlaps with exons in another gene.
    Changes the given exon to restrict it to the the non-overlapping section.
        args: 
            geneDict: dictionary where the keys are gene names and the values are Gene objects
            line: line from two intersected gtfs with the genes of interest listed first (as pybedtools Interval)
        return:
            Nothing: Modifies geneDict object.
    """
    # get exon from correct gene and increment its counter of overlaps in other genes
    gene = line.attrs['gene_name']
    exonKey = line.chrom + ":" + str(line.start+1) + "-" + str(line.stop)
    otherChr = line[9]
    otherStart = line[12]
    otherStop = line[13]
    geneDict[gene].getExons()[exonKey].incrementOverlapDiffGene(otherChr, otherStart, otherStop)

def extractSelectedGenes(filename):
    """
    Takes in a white-space delimited text file with HGNC gene names in the first column.
    Optional 0/1 exclusion/inclusion value in the second column.
    Makes and returns a list of the genes with a 1 in the second column (or all genes if only 1 column).
        arg:
            filename: the path to the gene list file
        return:
            A list of gene names
    """
    selectedGenes = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip().split()
            selectedGenes.append(line[0])
    return selectedGenes

def gtfToGenes(filename, genesToKeep):
    """
    Reads in a gtf file and a list of genes. And makes Gene objects that contain the exons of those genes.
    Also processes when other genes in the gtf overlap the genes of interest.
        args:
            filename: path to gtf file
            genesToKeep: list of genes to build into objects
        return:
            Dict of Gene objects where the key is the gene name and the value is the Gene object
    """
    print("Extracting exons from that set of provided genes.", file = sys.stderr)
    gtfExonSubset = pybedtools.BedTool(filename).filter(lambda line: line.attrs['gene_name'] in genesToKeep and line[2] == "exon").saveas()
    
    # create the dictionary to hold the genes
    # gene name : Gene object
    geneDict = dict()
    
    # Process returned gtf lines into Genes (and therefore Exons)
    print("Parsing subsetted gtf.", file = sys.stderr)
    for gtfLine in gtfExonSubset:
        processGtfLine(geneDict, gtfLine)
    
    # Intersect subsetted gtf with entire gtf and retrieve overlaps with exons in different genes
    # Update Genes (more precisely, their exons) with this information
    print("Adding information about overlap with exons in other genes.", file = sys.stderr)
    # run intersection and process output
    gtfsIntersected = gtfExonSubset.intersect(pybedtools.BedTool(filename), wa = True, wb = True)
    # remove cases where the intersect is a non-exonic feature or with the same gene
    gtfsIntersectedFiltered = gtfsIntersected.filter(lambda i: i[11] == "exon" and ('"' + i.attrs['gene_name'] + '";') not in i[17])
    for line in gtfsIntersectedFiltered:
        processIntersection(geneDict, line)
    
    return geneDict


###################################################################################
# MAIN
###################################################################################
def main():
    # parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output', help = 'Output file. If omitted, outputs to stdout.')
    parser.add_argument('gtf', help = 'Gene annotation in gtf format (e.g., Gencode file)')
    parser.add_argument('genelist', help = 'White-space delimited text file with HGNC gene names in the first column.')
    args = parser.parse_args()
    
    # set output file to stdout if none was provided
    outputFile = args.output
    if outputFile is None:
        out = sys.stdout
    else:
        out = open(outputFile, 'w')
    
    # get gene names from gene list file
    selectedGenes = extractSelectedGenes(args.genelist)
    print("Read in " + str(len(selectedGenes)) + " genes.", file = sys.stderr)
    
    geneDict = gtfToGenes(args.gtf, selectedGenes)
    
    for gene in geneDict.values():
        print(str(gene), file = out, end = "")
    
    out.close()

if __name__ == '__main__':
    main()
