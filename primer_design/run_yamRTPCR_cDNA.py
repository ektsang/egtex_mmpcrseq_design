#!/usr/bin/env python

# author: Emily Tsang
# loosely based on script from Boxiang Liu
# adapted for my purposes and moved some data structure to be classes instead of dicts

# To design multiple multiplex PCR primer pools

from __future__ import print_function
import argparse
import sys, os
import timeit, datetime
import shlex, subprocess
import uuid
import random
import copy

###################################################################################
# Classes
###################################################################################
class PrimerPool(object):
    
    def __init__(self, name, primerPairList):
        self._name = name
        self._primerPairList = primerPairList
    
    def __repr__(self):
        result = ''
        for pair in self._primerPairList:
            result = result + self._name + '\t' + repr(pair) + '\n'
        result = result.strip() # strip last newline
        return result
    
    def size(self):
        return len(self._primerPairList)
    
    def addPrimerPairs(self, primerPairs):
        self._primerPairList = self._primerPairList + primerPairs
    
    def getPrimerPairs(self):
        return self._primerPairList

    def nPrimerPairs(self):
        return len(self._primerPairList)


class PrimerPair(object):
    
    def __init__(self, primerString, annotation):
        self._targetName, self._forward, self._ftemp, \
            self._reverse, self._rtemp, self._amplicon, self._length = primerString.strip().split()
        self._annotation = annotation
    
    def __repr__(self):
        fields = [self._targetName, self._forward, self._ftemp, self._reverse, \
                  self._rtemp, self._amplicon, self._length, self._annotation]
        result = '\t'.join(fields)
        return result
    
    def getTargetName(self):
        return self._targetName
    
    def getPrimers(self):
        """Return forward and reverse primer sequences."""
        return([self._forward, self._reverse])


class PoolCounter(object):
    """Simple object to maintain pool count and return simple pool names based on current count."""
    
    def __init__(self, count):
        self._count = count # this is last count used, next pool will be count + 1
        
    def getPoolName(self):
        self._count += 1
        return "Pool" + str(self._count)

    def decrement(self):
        self._count -= 1
    
###################################################################################
# FUNCTIONS
###################################################################################
def parsePools(filename):
    """
    Parse pools into a list of PrimerPool objects.
        args:
            filename: path to the file with the pools
        return:
            list of PrimerPools
    """
    infile = open(filename, 'r')
    resultList = []
    primerPairList = []
    prevPool = ''
    maxCount = 0
    
    for line in infile:
        line = line.strip().split()
        poolName = line[0]
        primerPair = PrimerPair('\t'.join(line[1:-1]), line[-1])
        if poolName != prevPool and prevPool != '':
            resultList.append(PrimerPool(prevPool, primerPairList))
            poolNumber = int(prevPool[4:])
            if poolNumber > maxCount:
                maxCount = poolNumber
            primerPairList = [primerPair]
        else:
            primerPairList.append(primerPair)
        prevPool = poolName
    
    # process last line
    resultList.append(PrimerPool(poolName, primerPairList))
    return resultList, PoolCounter(maxCount)


def parseSites(filename):
    """
    Parse sites into a dict with the site name as key and the site line (newline stripped) as value.
        args:
            filename: path to the file with the site information
                      has 5 tab-delimited columns (name, chr, start, end, strand) with 1-based indexing
        return:
            dict with site name as key
    """
    infile = open(filename, 'r')
    resultDict = dict()
    for line in infile:
        line = line.strip()
        splitLine = line.split()
        if line[0] in resultDict:
            print("Duplicate site name " + line[0] + ". Skipping later occurrence.", file = sys.stderr)
        resultDict[splitLine[0]] = line
    infile.close()
    return resultDict


def generateFileNames(directory):
    """
    Creates unique file names for input sites and included primers.
        args:
            directory: the start of the path for the file names
        returns:
            a tuple with the input file name and the included primers file name
    """
    fileSuffix = str(uuid.uuid4())
    inputName = directory + '/input.' + fileSuffix
    includedName = directory + '/included.' + fileSuffix
    return (inputName, includedName)


def callPrimerDesign(sitesList, includedPrimerList, annotation, directory, nCandidates):
    """
    Calls yamRTPCR_cDNA.pl script to get a new pool of multiplex primers.
        args:
            sitesList: list of the sites to include in primer design
            includedPrimerList: list of primer pairs that have already been included in the pool
            annotation: path to file to use as annotation
            directory: the directory in which to write the input files
            nCandidates: number of candidate primers to try per site
        return:
            tuple of stdout and stderr
    """
    assert(sitesList > 0) # sanity check
    ## prepare the command
    siteFilename, includedFilename = generateFileNames(directory) # generate a unique file names the input files
    command = "perl yamRTPCR_cDNA.pl --targetFile='" + siteFilename + "'"
    command = command + " --numberCandidatePrimers=" + str(nCandidates)
    # write sites to file
    with open(siteFilename, 'w') as siteFile:
        print('\n'.join(sitesList), file = siteFile)
    # deal with included primers if not an empty list
    if len(includedPrimerList) > 0:
        command = command + " --includedPrimers='" + includedFilename + "'"
        # make included primer file
        includedList = []
        for i in includedPrimerList:
            includedList.extend(i.getPrimers())
        with open(includedFilename, 'w') as includedFile:
            print('\n'.join(includedList), file = includedFile)
    # deal with non default annotation if provided
    if annotation != "default":
        command = command + " --annotation='" + annotation + "'"

    ## run the command
    print('-'*40)
    print("Running " + command)
    args = shlex.split(command)
    result = subprocess.Popen(args, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
    ## output entire stderr and stdout of primer design script to stdout
    print(result[1])
    print(result[0])
    print('-'*40)
    return (result[0], result[1])


def parseOutput(output, annotation):
    """
    Parse stdout from primer design script.
        args:
            output: the string of stdout from running the primer design script
            annotation: name of the annotation associated with this primer pair
        return:
            list of PrimerPairs
    """
    if output == '':
        return []
    output = output.strip().split('\n')
    pairsList = [PrimerPair(pairString, annotation) for pairString in output]
    return pairsList


def parseError(error):
    """
    Parse stderr from primer design script.
        args:
            error: the string of stderr from running the primer design script
        return:
            tuple of lists (1) sites that were found in no transcript,
                           (2) sites removed because of no primers, 
                           (3) primers removed because they failed filters,
                           (4) primers removed because of bad interactions
    """
    error = error.strip().split('\n')
    noTranscriptList = []
    noPrimerList = []
    filterFailedList = []
    badInteractList = []
    
    for line in error:
        if "No transcript containing" in line:
            noTranscriptList.append(line.split(',')[0].split()[-1])
        elif "Primer3 couldn't find any primer for" in line:
            noPrimerList.append(line.split(',')[0].split()[-1])
        elif "input sequence for primer3 is shorter" in line:
            noPrimerList.append(line.split()[0])
        elif "No good candidate primers could be found" in line:
            filterFailedList.append(line.split()[-1].split('!')[0])
        elif "ejected" in line:
            badInteractList.append(line.split()[1])
    
    return (noTranscriptList, noPrimerList, filterFailedList, badInteractList)


def completePool(siteDict, pool, poolSize, annotation, directory, nCandidates = 10, multFactor = 1.5):
    """
    Tries to complete the provided pool.
    Will not shrink the pool size if initially larger than the max pool size.
        args:
            siteDict: dict of sites to add to pools
            pool: primer pool to complete, may be empty
            poolsSize: maximal pool size
            annotation: path to annotation file
            directory: path to directory for output
            nCandidates: number of candidate primers per site
            multFactor: number of times more sites to try than spots to fill in the pool
        return:
            tuple of (1) the designed PrimerPool, 
                     (2) the set of sites for which no primers could be designed,
                     (3) the set of other sites not used
    """
    availableSites = set(siteDict.keys())
    noPrimerSites = set()
    interactingSites = set()
    initialSize = pool.size()
    
    # return unchanged objects if the pool is already big enough
    if initialSize >= poolSize:
        print("WARNING: completePool() called with pool bigger than requested size. Returning.", file = sys.stderr)
        return (pool, dict(), siteDict)
    
    # try to design new primers, stop when the pool is big enough
    print("current pool size: " + str(pool.size()))
    print("number remaining sites: " + str(len(availableSites)))
    ntries = 0
    while pool.size() < poolSize and len(availableSites) > 0:
        ntries += 1
        numNeeded = poolSize - pool.size()
        chosen = random.sample(availableSites, min(int(numNeeded * multFactor), len(availableSites)))
        sitesList = [siteDict[site] for site in chosen]
        # run and process output from primer design script
        output, error = callPrimerDesign(sitesList, pool.getPrimerPairs(), annotation, directory, nCandidates)
        newPairs = parseOutput(output, annotation)
        noTrans, noPrimer, filtered, badInteract = parseError(error)
        if len(noTrans) > 0:
            print("", file = sys.stderr)
        failed = noTrans + noPrimer
        notThisPool = filtered + badInteract
        assert(len(sitesList) == (len(newPairs) + len(failed) + len(notThisPool))) # sanity check that all cases were handled
        # deal with sites for which design is not possible in this pool (or at all)
        noPrimerSites.update(failed)
        interactingSites.update(notThisPool)
        availableSites.difference_update(failed + notThisPool)
        # deal with newly designed sites
        if pool.size() + len(newPairs) > poolSize:
            newPairs = newPairs[:numNeeded]
        pool.addPrimerPairs(newPairs)
        availableSites.difference_update([pair.getTargetName() for pair in newPairs])
    
    # put together the set of sites to try on other pools
    availableSites.update(interactingSites)
    assert((pool.size() + len(availableSites) + len(noPrimerSites)) == (len(siteDict) + initialSize)) # sanity check
    print("Number of calls to primer design script: " + str(ntries), file = sys.stderr)
    return pool, noPrimerSites, availableSites


def completePrimerPools(siteDict, poolList, poolSize, annotation, directory):
    """
    Tries to complete the provided pools.
    Will not shrink the pool size of existing pools that are larger than the max pool size.
        args:
            siteDict: dict of sites to add to pools
            poolList: existing set of pools, some of which may be compete
            poolsSize: maximal pool size
            annotation: path to annotation file
            directory: path to directory for output
        return:
            tuple of the updated pool list and the dict of the remaining sites
    """
    completedPools = []
    startedPools = []
    remainderDict = copy.deepcopy(siteDict)
    noPrimerSites = set()
    
    # sort completed pools from started ones
    for pool in poolList:
        if pool.size() >= poolSize:
            completedPools.append(pool)
        else:
            startedPools.append(pool)
    
    # try to complete incomplete pools
    for spool in startedPools:
        # once we run out of sites, keep pools as they were
        if len(remainderDict) == 0:
            completedPools.append(spool)
        updatedPool, npSites, remainderSites = completePool(remainderDict, spool, poolSize,
                                                            annotation, directory, nCandidates = 30, multFactor = 1)
        completedPools.append(updatedPool)
        noPrimerSites.update(npSites)
        # reduce remainder dict to the remaining sites
        remainderDict = {k: remainderDict[k] for k in remainderSites}
    
    # combine remainderDict and no primer sites before returning
    noPrimerDict = {k: siteDict[k] for k in noPrimerSites}
    remainderDict.update(noPrimerDict)
    assert(len(completedPools) == len(poolList)) # sanity check
    return (completedPools, remainderDict)


def makePrimerPools(siteDict, poolCounter, poolSize, numPools, annotation, directory, minSize = 4):
    """
    Tries to make up to the specified number of pools.
    Discards pools that are too small be be good starter pools.
        args:
            siteDict: dict of sites to add to pools
            poolCounter: PoolCounter object used to get the next pool name
            poolsSize: maximal pool size
            numPools: max number of pools to make
            annotation: path to annotation file
            directory: path to directory for output
            minSize: minimum pool size for it to be a decent starter pool.
        return:
            tuple of the pool list and the dict of the remaining sites
    """
    newPoolList = []
    remainderDict = copy.deepcopy(siteDict)
    noPrimerSites = set()
    numTimesTooSmall = 0 # track the number of times in a row the pool comes back too small
    i = 0
    while i < numPools and len(remainderDict) > 0 and numTimesTooSmall <= 5:
        newPool, npSites, remainderSites = completePool(remainderDict, PrimerPool(poolCounter.getPoolName(), []),
                                                        poolSize, annotation, directory)
        if newPool.nPrimerPairs() >= minSize:
            newPoolList.append(newPool)
            numTimesTooSmall = 0
            i += 1
        else:
            numTimesTooSmall += 1
            if newPool.nPrimerPairs() > 0:
                sitesSmallPool = [pair.getTargetName() for pair in newPool.getPrimerPairs()]
                remainderSites.update(sitesSmallPool) # add in sites from pools that were too small
            else:
                i = numPools # to break out of the loop after completing this iteration
            poolCounter.decrement()
        noPrimerSites.update(npSites)
        # reduce remainder dict to the remaining sites
        remainderDict = {k: remainderDict[k] for k in remainderSites}

    # combine remainderDict and no primer sites before returning
    noPrimerDict = {k: siteDict[k] for k in noPrimerSites}
    remainderDict.update(noPrimerDict)
    return (newPoolList, remainderDict)


def outputPools(poolList, filename):
    """
    Output pool list to the given file.
        args:
            poolList: list of PrimerPools to write to file
            filename: path to output file
        return:
            nothing
    """
    with open(filename, 'w') as outfile:
        for pool in poolList:
            print(pool, file = outfile)


def outputSites(siteDict, filename):
    """
    Output sites to the given file.
        args:
            siteDict: dict of sites to write to file
            filename: path to the output file
        return:
            nothing
    """
    with open(filename, 'w') as outfile:
        for site in sorted(siteDict.keys()):
            print(siteDict[site], file = outfile)


###################################################################################
# MAIN
###################################################################################

def main():
    # for reproducibility
    random.seed(1)
    
    # start timer
    tic = timeit.default_timer()

    # change directory to the script's directory
    # (from http://stackoverflow.com/questions/1432924/python-change-the-scripts-working-directory-to-the-scripts-own-directory)
    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir(dname)
    
    # process command line input
    parser = argparse.ArgumentParser("\nIMPORTANT NOTE: all paths should be absolute or relative to the script's directory")
    parser.add_argument('--poolSize', type = int, default = 12,
                        help = 'Target pool size.')
    parser.add_argument('--numberPools', type = int, default = 48,
                        help = 'Maximum number of pools.')
    parser.add_argument('--geneAnnotation', nargs = '+',
                        help = 'Ordered list of gene annotation files with subset of exons or transcripts of interest.\
                        If not provided, the annotation file hard-coded in the primer design script will be used.\
                        Make sure that these are a subset of the GENE_ANNOTATION_FULL in the primer design script.')
    parser.add_argument('--existingPools',
                        help = 'Path to file with existing pools. In format output by outputPools().')
    parser.add_argument('--counterStart', type = int, default = 0,
                        help = 'Number at which to start naming pools. Will be ignored if existingPools is provided.')
    parser.add_argument('--directory',
                        help = 'Directory in which to write the input files for the primer design script.')
    parser.add_argument('--inputSiteFile', required = True,
                        help = 'Input file with variants for which we want primers.')
    args = parser.parse_args()

    # set gene annotation to "default" if not provided
    if not args.geneAnnotation:
        annotations = ["default"]
    else:
        annotations = args.geneAnnotation

    # if existing pools were provided, read them in
    if args.existingPools:
        poolList, counter = parsePools(args.existingPools)
    else:
        poolList = []
        counter = PoolCounter(args.counterStart)
    
    # read in sites
    siteDict = parseSites(args.inputSiteFile)

    for annotation in annotations:
        # try to complete existing pools, if any
        if len(poolList) > 0 and len(siteDict) > 0:
            poolList, siteDict = completePrimerPools(siteDict, poolList,
                                                     args.poolSize, annotation, args.directory)

        # make new pools if need more pools and try to complete them
        if len(poolList) < args.numberPools and len(siteDict) > 0: 
            newPools, siteDict = makePrimerPools(siteDict, counter, args.poolSize,
                                                 args.numberPools - len(poolList),
                                                 annotation, args.directory)
            newPools, siteDict = completePrimerPools(siteDict, newPools, args.poolSize,
                                                     annotation, args.directory)
            poolList = poolList + newPools

    # output pools
    outputPools(poolList, args.directory + "/pools.txt")
    outputSites(siteDict, args.directory + "/remaining.sites.txt")

    # stop timer
    toc = timeit.default_timer()
    elapsed = str(datetime.timedelta(seconds = toc - tic))
    
    print("Total time elapsed: " + elapsed)

if __name__ == '__main__':
    main()
