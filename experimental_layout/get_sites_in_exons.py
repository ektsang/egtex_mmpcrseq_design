#!/usr/bin/env python

# author: Emily Tsang

# find the variants that fall within the exons of selected genes
# and that pass various filters

from __future__ import print_function
from __future__ import division
from collections import defaultdict
import sys
import copy
import re
import itertools
import argparse
import pybedtools
import pysam
import retrieve_genes
import retrieve_sites

###################################################################################
# CLASSES
###################################################################################

class VariantPlus(retrieve_sites.Variant):
    """
    Variant that also contains extra information 
    (e.g., orthogonal genotyping information and VEP annotation).
    """
    def __init__(self, variantArgs, hetListExtra, imputedExtra, certaintyExtra, vep = ''):
        super(VariantPlus, self).__init__(*variantArgs)
        self._hetListExtra = hetListExtra # None if there was no extra entry for this variant
        self._imputedExtra = imputedExtra
        self._certaintyExtra = certaintyExtra
        self._vep = vep
        
    def __eq__(self, other):
        """ Note: checks equality of coordinates and of hets, but not of vep string."""
        eq = super(VariantPlus, self).__eq__(other)
        if self._hetListExtra is None and other._hetListExtra is None:
            return eq
        elif self._hetListExtra is not None and other._hetListExtra is not None:
            return eq and set(self._hetListExtra) == set(other._hetListExtra)
        else:
            return False # exactly one of the hetListExtras is None
        
    def __hash__(self):
        extra = [] if self._hetListExtra is None else self._hetListExtra
        return (self._chr, self._pos, '.'.join(self._hetList), '.'.join(extra)).__hash__()
    
    def hasExtra(self):
        """ Return True if the variant has the extra data from a second genotyping source."""
        return (self._hetListExtra is not None)
    
    def getHets(self):
        """ Return the combined list of hets from both the main and the extra het lists."""
        if self._hetListExtra is not None:
            hets = list(set(self._hetList) | set(self._hetListExtra))
        else:
            hets = self._hetList
        return hets
    
    def nhet(self):
        """ Return the combined number of hets from both the main and the extra het lists. """
        return len(self.getHets())
    
    def nhetPrimary(self):
        return super(VariantPlus, self).nhet()
    
    def nhetExtra(self):
        return len(set(self._hetListExtra)) if self._hetListExtra is not None else 0
    
    def addVEP(self, vepstring):
        if self._vep != '':
            print("For variant: " + str(self), file = sys.stderr)
            print("Replacing non-empty VEP string from " + self._vep + " to " + vepstring + "!!", file = sys.stderr)
        self._vep = vepstring
        
    def prettyString(self, printImpute = True):
        if self._vep == '':
            print("Missing VEP string for variant:" + super(VariantPlus, self).prettyString(printImpute), file = sys.stderr)
            vep = 'NA'
        else:
            vep = self._vep
        value = super(VariantPlus, self).prettyString(printImpute)
        # assuming that the extra data is definitely genotyping data if it exists
        # and therefore that the imputed and certainty fields are not None
        if self.hasExtra():
            value = value + "\t" + ",".join(self._hetListExtra)
            imputedStr = "imp" if self._imputedExtra else "gt"
            value = value + "\t" + imputedStr + "\t" + str(self._certaintyExtra) + "\t" + vep
        else:
            value = value + "\tNA\tNA\t" + vep
        return value
    
    def prettyStringLimited(self):
        vep = self._vep if self._vep != '' else 'NA'
        value = self._chr + "\t" + "\t".join( map(str, [self._pos-1, self._pos, self.nhet()]))
        if self.hasExtra():
            imputedStr = "imp" if self._imputedExtra else "gt"
            value = value + "\t" + imputedStr + "\t" + str(self._certaintyExtra) + "\t" + vep
        else:
            value = value + "\tNA\tNA\t" + vep
        return value

class Transcript(object):
    """ Basic information about transcript abundance and ordering."""
    def __init__(self, name, medianTPM, madTPM, nnzTPM, medianPerc, madPerc, nnzPerc, rank):
        self.name = name
        self.medianTPM = medianTPM
        self.madTPM = madTPM
        self.nNonZeroTPM = nnzTPM
        self.medianPerc = medianPerc
        self.madPerc = madPerc
        self.nNonZeroPerc = nnzPerc
        self.rank = rank
        
    def __repr__(self):
        value = map(str, [self.name, self.rank, self.medianTPM, self.madTPM, self.nNonZeroTPM,\
                          self.medianPerc, self.madPerc, self.nNonZeroPerc])
        return "\t".join(value)
    
    def __str__(self):
        return repr(self)

    
###################################################################################
# FUNCTIONS
###################################################################################
def parseIntersection(intersection, varDict):
    """
    Parse the gene-variant intersection BedTool into gene/variant/exon objects.
        args:
            intersection: BedTool object of the intersection
            varDict: dict of orthogonal variant information. keys are "chr:pos".
            annoDict: dict of vep annotations. dict of dicts with "chr:pos0:pos1" as first key and gene name as second key.
        return:
            dict of Genes, which contain Exons and Variants
    """
    # build objects as iterate through intervals
    geneDict = dict()
    for interval in intersection:
        # first make gene if necessary
        geneName = interval[7]
        if geneName not in geneDict:
            geneDict[geneName] = retrieve_genes.Gene(*interval[7:11])
        # build up the variant object
        variantArgs = retrieve_sites.bedLineToVariantArgs(interval[11:])
        # get extra variant information
        varkey = interval[11] + ":" + interval[13]
        if varkey not in varDict:
            extraArgs = [None, None, None]
        else:
            extra = varDict[varkey]
            extraArgs = [extra._hetList, extra._imputed, extra._certainty]
        # actually instantiate variant
        variant = VariantPlus(variantArgs, *extraArgs)
        # add exon (called function takes care of checking whether exon already exists in the gene)
        geneDict[geneName].addExonComplete(variant, *(interval[0:7]))
    return geneDict

def parseTranscripts(transcriptFile):
    """
    Read transcript file into a dict of transcripts where the key is the transcript id.
        args:
            transcriptFile: path to file with transcript information.
        return:
            dict with transcript objects keyed by the transcript id.
    """
    transDict = dict()
    tfile = open(transcriptFile, 'r')
    # skip header
    assert(tfile.readline() == 'gene_id\ttranscript_id\tmedian.tpm\tmad.tpm\ttpm.gt0\tmedian.perc\tmad.perc\tperc.gt0\tHGNC\trank\n')
    # process body
    for line in tfile:
        line = line.strip().split()
        tid = line[1]
        if tid in transDict:
            print("Already processed transcript " + tid + ". Skipping.", file = sys.stderr)
            continue
        transDict[tid] = Transcript(tid, float(line[2]), float(line[3]), int(line[4]), float(line[5]), \
                                    float(line[6]), int(line[7]), int(line[9]))
    tfile.close()
    return transDict

def exportSiteCountsByFreqPerGene(geneDict, outfile, breaks = None, names = None):
    """
    Extracts from geneDict info to create a file with 3 columns:
    gene hetGroup #sites
    where hetGroup is one of 1, 2-10, 11-20, 21-40, 41-100 (or specified as argument)
        args:
            geneDict: dict of genes with necessary information (constructed by parseIntersection)
            outfile: file to which output will be written
            breaks: definition of intervals (modeling the cut function in R)
                for the default param the intervals are (0,1], (1,10], (10,20], (20,40], (40,100]
            names: string representation of the breaks. Should be provided if breaks is provided.
        return:
            Nothing. Outputs directly to file.
    """
    if breaks is None:
        breaks = [0,1,10,20,40,100]
    if names is None:
        names = ['1', '2-10', '11-20', '21-40', '>40']
    out = open(outfile, 'w')
    # print header
    print("gene\thetGroup\tcount", file = out)
    assert len(names) == (len(breaks) - 1)
    for gene in geneDict.values():
        bins = [0] * len(names)
        variantList = gene.getVariants()
        for variant in variantList:
            nhet = variant.nhet()
            assigned = False
            for i, (l,u) in enumerate(zip(breaks[:-1], breaks[1:])):
                if l < nhet <= u:
                    bins[i] += 1
                    assigned = True
                    break
            if not assigned:
                print('Warning: number of hets (' + str(nhet) + ') not in breaks.', file = sys.stderr)
        # print out info for that gene
        for name, count in zip(names, bins):
            print("\t".join([gene.getName(), name, str(count)]), file = out)
    out.close()
    return

def getVariantsCoveringIndividuals(variantList, maxVariants = 30, inds = None):
    """
    For the given gene, retrieve a set of sites such that the largest possible number of individuals are het in at least one site.
    Return up to the maximum number of sites.
    Only add sites when they add new het individuals.
        args:
            variantList: the list of variants to choose from, make sure this is a copy becuase it will be mutated
            maxVariants: the maxinal number of Variants to return
            inds: the list of individuals already selected
        return:
            A list of selected variants, in order, where the first is the first selected
            and the set of individuals heterozygous for at least one site.
    """
    if inds is None:
        inds = set()
    else:
        inds = set(inds)
    chosenList = []
    # select the site that covers the most new individuals
    # repeat until no new individuals are covered or all variants are selected
    chosenIndex = -1
    numAdded = 0
    while(len(variantList) > 0 and len(chosenList) < maxVariants):
        change = False
        # find site, remove it once added
        for i, variant in enumerate(variantList):
            newInds = set(variant.getHets()) - inds
            if len(newInds) > numAdded:
                chosenIndex = i
                numAdded = len(newInds)
                change = True
        # once new variants are no longer being added, sort by number of hets
        if not change:
            chosenList.extend(sorted(variantList, key = lambda v: -len(v.getHets())))
            break
        chosenVariant = variantList.pop(chosenIndex)
        chosenList.append(chosenVariant)
        # process this variant
        inds = inds.union(set(chosenVariant.getHets()))
        numAdded = 0
        chosenIndex = -1
    return chosenList, inds

def exportVariantsCoveringIndividuals(geneDict, outfile):
    """
    For each gene, retrieve a set of sites such that the largest possible number of individuals are het in at least one site.
    Write this information to the provided file with 6 columns:
    gene, chr of site, pos of site, the order position in which the site was picked, 
    the cumulative number of people with hets once that site is picked, the number of het individuals at that site.
    Uses a greedy approach to select sites.
        args:
            geneDict: dict of genes with necessary information (constructed by parseIntersection)
            outfile: file to which output will be written
        return:
            Nothing. Outputs directly to file.
    """
    out = open(outfile, 'w')
    # print header
    print("gene\tsite.chr\tsite.pos\tpick.order\tcum.num.people\tnum.people", file = out)
    for gene in geneDict.values():
        varList = gene.getVariants()
        chosen, hets = getVariantsCoveringIndividuals(varList)
        cumulative = []
        for index, variant in enumerate(chosen):
            hets = variant.getHets()
            cumulative.extend(hets)
            cumulative = list(set(cumulative))
            print("\t".join([gene.getName(), variant.getChr(), str(variant.getPos()), str(index+1), \
                             str(len(cumulative)), str(len(hets))]), file = out)
    out.close()
    return

def isAcceptableTranscript(transcript, transDict, maxRank = 6, minPerc = 10, minPosTPM = 8, minPosPerc = 8):
    """
    Decides whether or not a given transcript is common enough to be used to design on.
    Requires that the transcript be ranked well, be used often enough on average,
    have TPM > 0 in a minimum number of tissues and be used (by 50% of individuals) in a minimum number of tissues.
        args:
            transcript: tuple of the name of the transcript and the transcript type
            transDict: the dict keyed by transcript names that contains transcript objects
            maxRank: the maximum rank of the transcript, where transcripts are rank by median usage across tissues
            minPerc: the minimum (median + MAD) percent usage of transcript across tissues
            minPosTPM: the minimum number of tissues with a nonzero median TPM across individuals
            minPosPerc: the minimum number of 
        return:
            True if the transcript passes all criteria.
    """
    passes = False
    transName = transcript[0]
    if transName in transDict:
        transcript = transDict[transName]
        if (transcript.rank <= maxRank) and (transcript.medianPerc + transcript.madPerc >= minPerc) and \
           (transcript.nNonZeroTPM >= minPosTPM) and (transcript.nNonZeroPerc >= minPosPerc):
            passes = True
    return passes

def removeUnsuitableVariantsAndExons(geneDict, transDict, siteMinHet = 6, geneMinHet = 10):
    """
    Run through the dict and remove exons that are in unused transcripts. 
    Remove sites with too few het individuals.
    Remove exons for which all variants are removed and genes for which all exons are removed.
    Remove genes for which the are two few het individuals overall
        args:
            geneDict: the data structure that contains all the gene, exon, variant information
            transDict: the dict of transcript information
            siteMinHet: the minimal number of het individuals for a variant to be kept
            geneMinHet: the minimal number of het individuals overall variants in a gene
        return:
            A dict of the removed genes, 
            where the key is the gene name and the value is a list of the reasons for removal.
    """
    # initialize dict
    removedGenes = dict()
    for genekey, gene in geneDict.items():
        removedExons = []
        reasons = set()
        for exonkey, exon in gene.getExons().items():
            # remove any exon that is not part of a commonly used transcript
            trans = exon.getTranscripts()
            acceptable = [isAcceptableTranscript(t, transDict) for t in trans]
            if not any(acceptable):
                removedExons.append(exonkey)
                reasons.add('uncommon transcript')
                continue
            # even if the exon had some acceptable transcripts, remove the offending ones
            else:
                removedTranscripts = [t[0] for t, accept in zip(trans, acceptable) if not accept]
                exon.removeTranscripts(removedTranscripts)
            # remove any variants that aren't frequent enough or that overlap another gene
            removedVariants = []
            overlapsDiff = exon.numOverlapsDiff() > 0
            if overlapsDiff:
                # first get a list of the exon coordinates as ('chr', pos0, pos1) tuples
                overlapping = [re.compile("[:-]").split(e) for e in exon.getOverlapDiffList()]
                overlapping = [(triplet[0], int(triplet[1]), int(triplet[2])) for triplet in overlapping]
            for variant in exon.getVariants():
                if variant.nhet() < siteMinHet:
                    removedVariants.append(variant)
                elif overlapsDiff:
                    varPos = variant.getPos()
                    if any([variant.getChr() == chrom and varPos > pos0 and varPos <= pos1 for chrom,pos0,pos1 in overlapping]):
                        removedVariants.append(variant)
            for v in removedVariants:
                exon.removeVariant(v)
            # remove the exon if all its variants have been removed
            if exon.nvariants() == 0:
                removedExons.append(exonkey)
                reasons.add('no common unique-to-this-gene variants')
        # actual exon removal
        for ekey in removedExons:
            gene.removeExon(ekey)
        # remove gene if all of its exons were removed
        if gene.numExons() == 0:
            removedGenes[genekey] = list(reasons)
            continue
        # remove gene if it too few het individuals overall by picking the 5 best sites
        selected, hetInds = getVariantsCoveringIndividuals(gene.getVariants(), maxVariants = 5)
        if len(hetInds) < geneMinHet:
            removedGenes[genekey] = ['too few hets overall'] + list(reasons)
    # actual gene removal
    for gkey in removedGenes.keys():
        del geneDict[gkey]
    return removedGenes

def addVepAnnotation(geneDict, vepFile):
    """
    Add VEP variant annotations to variants nested in the gene dict.
        args:
            geneDict: the data structure that contains all the gene, exon, variant information. 
            vepFile: the path to the (tabixed) vcf file with the vep annotation.
        return:
            Nothing. Modifes the gene dict in place.
    """
    regionTabix = pysam.Tabixfile(vepFile,'r')
    for geneName, gene in geneDict.items():
        for exon in gene.getSortedExons():
            for variant in exon.getVariants():
                picked = False
                vepStrings = []
                chrom = str(variant.getChr()[3:]) # str so that it's not unicode
                pos = variant.getPos()
                transcripts = exon.getTranscriptNames()
                # remove everything after the period in the transcript names
                transcripts = [t.split(".")[0] for t in transcripts]
                lines = regionTabix.fetch(chrom, pos-1, pos)
                # iterate through matching records until the one that matches the SNP exactly
                for line in lines:
                    line = line.strip().split()
                    if line[0] == chrom and line[1] == str(pos) and len(line[3]) == 1 \
                       and all([len(allele)==1 for allele in line[4].split(',')]):
                        if picked:
                            print("Already came across a line that matched the variant. Skipping.", file = sys.stderr)
                            continue
                        picked = True
                        vepgroups = line[-1].split(";")[-1].split("=")[1].split(",")
                        # VEP annotation is the 2nd field, gene name is the 4th (1-indexed)
                        for vepgrp in vepgroups:
                            vepfields = vepgrp.split('|')
                            if vepfields[3] == geneName and vepfields[6] in transcripts:
                                vepStrings.append(vepfields[1])
                        vepStringCombined = ",".join(set(vepStrings))
                        variant.addVEP(vepStringCombined)

def orderAndExportSites(geneDict, transDict, outfile):
    """
    Prioritize sites in each gene first by transcript abundance then by number of new hets they add.
    Print sites to provided output file.
        args:
            geneDict: data structure with the gene, exon, variant information.
            transDict: dict with the transcript information.
            outfile: path to the output file
        return:
            Nothing.
    """
    out = open(outfile, 'w')
    # print header
    header = "varname\tHGNC\tENSG\tstrand\tvar.chr\tvar.pos0\tvar.pos1\tnhet\tgt\tcertainty\tVEP\tintra.exon.plausible"
    header = header + "\texon.chr\texon.pos0\texon.pos1\toverlap.same\toverlap.other"
    header = header + "\ttype\ttranscript\trank\tmedian.TPM\tmad.TPM\tnon.zero.TPM\tmedian.perc\tmad.perc\tnon.zero.perc"
    print(header, file = out)
    for geneName, gene in geneDict.items():
        varCounter = 0
        # get transcript names with ranks in dict
        transcripts = {t:transDict[t].rank for t in gene.getTranscriptNames()}
        # get order in which to process exons
        rankDict = defaultdict(list)
        exons = gene.getExons()
        for exonkey, exon in exons.items():
            rank = reduce(min, [(transcripts[t],t) for t in exon.getTranscriptNames()]) # best transcript for that exon
            rankDict[rank].append(exonkey)
        # keep track of which variants have been processed and which individuals have hets
        processedVariants = set()
        coveredIndividuals = set()
        # for each rank, process the given exons together
        for r, tname in sorted(rankDict.keys()):
            exonList = [exons[ekey] for ekey in rankDict[r, tname]]
            variantList = list(itertools.chain.from_iterable([e.getVariants() for e in exonList]))
            selected, hetInds = getVariantsCoveringIndividuals(variantList, maxVariants = len(variantList), inds = coveredIndividuals)
            coveredIndividuals = coveredIndividuals.union(hetInds)
            # print chosen, taking care to identify their exons
            # only printing variants that aren't already in the processed variants list
            for v in selected:
                if v in processedVariants:
                    continue
                varCounter = varCounter + 1
                # gene info
                outLine = "\t".join([geneName + "_" + str(varCounter), geneName, gene.getENSG(), gene.getStrand()])
                # variant info
                outLine = outLine + "\t" + v.prettyStringLimited()
                # exon info
                found = False
                for ex in exonList:
                    if v in ex.getVariants():
                        if found: # means there were two exons with the variant - shouldn't happen
                            print("Variant found in more than one exon of the same transcript:", file = sys.stderr)
                            print(v, file = sys.stderr)
                            print("Not processing extra occurrence.", file = sys.stderr)
                            continue
                        found = True
                        chrom, start, end = ex.getCoordinates()
                        ttype = [t[1] for t in ex.getTranscripts() if t[0] == tname]
                        assert(len(ttype) == 1)
                        ttype = ttype[0]
                        outLine = outLine + "\t" + ("1" if plausibleVariant(v, ex) else "0")
                        outLine = outLine + "\t" + chrom + "\t" + str(start - 1) + "\t" + str(end)
                        same = ex.getOverlapSameList()
                        diff = ex.getOverlapDiffList()
                        if len(same) > 0:
                            outLine = outLine + "\t" + ";".join(same)
                        else:
                            outLine = outLine + "\tNA"
                        if len(diff) > 0:
                            outLine = outLine + "\t" + ";".join(diff)
                        else:
                            outLine = outLine +"\tNA"
                        outLine = outLine + "\t" + ttype
                # transcript info
                transcript = transDict[tname]
                assert(transcript.rank == r)
                outLine = outLine + "\t" + str(transcript)
                print(outLine, file = out)
            # update the processed variants list with the new set of variants
            processedVariants = processedVariants.union(set(selected))
        # sanity check that we printed the right number of variants (each variant once)
        assert(varCounter == len(gene.getVariants()))
    out.close()

def plausibleVariant(variant, exon, primerLength = 24, minAmpliconLength = 150):
    """
    A quick evalutation of whether the variant could have primers its exon.
        args:
            variant: Variant object to check
            exon: Exon that contains this variant
            primerLength: Assumed length of the primer
            minAmpliconLength: Minimum distance between the beginning of the forward primer and the end of the reverse primer
        return:
            True if the variant is sufficiently far from the exon boundaries and the exon is large enough
    """
    plausible = False
    varPos = variant.getPos()
    exonChr, exonStart, exonEnd = exon.getCoordinates()
    if min(abs(varPos - exonEnd), abs(varPos - exonStart)) >= primerLength and exon.length() >= minAmpliconLength:
        plausible = True
    return plausible

###################################################################################
# MAIN
###################################################################################
def main():
    ## Process input
    # parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--extravar', help = 'Additional processed variant list, in bed format, as created by retrieve_variants.py.')
    parser.add_argument('--vep', help = 'VCF file with VEP annotation.')
    parser.add_argument('-o','--output', help = 'Output prefix. If omitted, defaults to working directory.')
    parser.add_argument('variants', help = 'Processed variant list, in bed format, as created by retrieve_variants.py.')
    parser.add_argument('genes', help = 'Processed gene/exon information, in bed format, as created by retrieve_genes.py.')
    parser.add_argument('transcripts', help = 'Processed transcript information, as created by process_transcripts.R.')
    args = parser.parse_args()

    # get output file prefix
    outPrefix = args.output
    if outPrefix is None:
        outPrefix = ''

    # read in input files
    print("Reading input files and intersecting them...", file = sys.stderr)
    variantBed = pybedtools.BedTool(args.variants)
    geneBed = pybedtools.BedTool(args.genes)

    intersected = geneBed.intersect(variantBed, wa = True, wb = True)

    # make transcript dict
    transcriptDict = parseTranscripts(args.transcripts)
    
    # make extra variant annotation dict
    # leave them empty if those files were not provided
    variantDict = dict()
    if args.extravar is not None:
        variantDict = retrieve_sites.bedToVariantDict(args.extravar)
        
    # parse intersection back into objects, linking them together
    print("Parsing intersection into objects...", file = sys.stderr)
    selectedGeneDict = parseIntersection(intersected, variantDict)
    # remove variantDict now that the information has been included in geneDict
    del variantDict
    
    ## This is for some diagnostic plotting. Kept so as to not break those downstream scripts.
    print("Outputting processed files...", file = sys.stderr)
    # write file with het counts per gene, binned into frequency groups
    exportSiteCountsByFreqPerGene(selectedGeneDict, outPrefix + ".hetcounts.byGene.txt")
    # write file with sites that maximize the number of people with at least 1 het
    exportVariantsCoveringIndividuals(selectedGeneDict, outPrefix + ".orderedSites.byGene.txt")
    
    ## Do some filtering of genes/sites based of various criteria and other postprocessing
    # keep track of genes that drop out
    print("Filtering out genes with problematic exons or too few hets...", file = sys.stderr)
    genesDropped = removeUnsuitableVariantsAndExons(selectedGeneDict, transcriptDict)
    # print dropped genes
    print("Removed genes:", file = sys.stderr)
    for g in sorted(genesDropped.keys()):
        print(g + "\t" + ",".join(genesDropped[g]), file = sys.stderr)

    # add vep annotation
    addVepAnnotation(selectedGeneDict, args.vep)

    # order sites for each gene
    # use transcript information to select and order sites
    print("Selecting sites by transcript and number of hets...", file = sys.stderr)
    orderAndExportSites(selectedGeneDict, transcriptDict, outPrefix + ".selected.variants.ordered.txt")

if __name__ == '__main__':
    main()
