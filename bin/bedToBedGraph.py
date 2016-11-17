#!/usr/bin/env python

import os
import sys
import getopt
import csv
import subprocess
import collections
from os.path import basename
import resource
import datastructures.cache as cache

def do_chrm(bedfile, outfile, millionsReads):
    
    with open(bedfile) as ifile:
    
        infile = csv.reader(ifile, delimiter="\t")
        
        chr = None
        
        previousPosition = None
        currentValue = 0
        
        waitingDecrements = collections.defaultdict(int)
        
        rowCount = 0
        
        for row in infile:
            
            rowCount += 1
            
            chr = row[0]
            start = int(row[1])
            stop = int(row[2])
        
            assert start <= stop, "Check for Start after End"
        
            # preface - go through and output a line up to here
            
            for decrementPosition in sorted(waitingDecrements.keys()):
                if decrementPosition <= start:
                    assert currentValue>0
                    assert previousPosition<decrementPosition
                    outfile.writerow([chr, previousPosition, decrementPosition, currentValue/millionsReads])
                    previousPosition = decrementPosition
                    currentValue -= waitingDecrements[decrementPosition]
                    
                    if currentValue == 0:
                        previousPosition = None
                    
                    del waitingDecrements[decrementPosition]
                else:
                    break # no need to keep looping        
    
            # if we have a previous position and we've just read in a new value then we need to output what our current status is
            # with exceptions:
            # * if we dont have a previous position yet (i.e. we are just starting a series of reads) and there may be more that start at
            #   the same position
            # * if this read has the same position as the previous read (i.e. it starts in the same place)
            assert previousPosition<=start, "Previous position is somehow after this one, is file sorted?: "+bedfile 
            if previousPosition!=None and previousPosition<start:
                assert currentValue>0
                outfile.writerow([chr, previousPosition, start, currentValue/millionsReads])
                
            # this indicates we have another read so we should increment our current value
            currentValue += 1
            
            previousPosition = start
            
            # and add the decrement position to the list
            waitingDecrements[stop]+=1
            
        # we have done almost all of them at this point but probably still have some decrement positions to do
        for decrementPosition in sorted(waitingDecrements):
            
            assert currentValue>0
            assert previousPosition<decrementPosition
            
            outfile.writerow([chr, previousPosition, decrementPosition, currentValue/millionsReads])
            
            previousPosition = decrementPosition
            currentValue -= waitingDecrements[decrementPosition]
            del waitingDecrements[decrementPosition]
            
        return rowCount

def file_len(fname):
    p = subprocess.Popen(['wc', '-l', fname], stdout=subprocess.PIPE,
                                              stderr=subprocess.PIPE)
    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    return int(result.strip().split()[0])

# GNU sort is much faster than Python sort.  Especially for big files.
def sortAndUniq(bedfile):
    sortchr = subprocess.Popen(args=["sort","-u", "-k", "1,1", "-k", "2,2n", "-k", "3,3n", "-k", "6,6"],  stdin=open(bedfile),  stdout=open(bedfile+".sorted", "w"))
    sortchr.wait()
    os.remove(bedfile)

# GNU sort is much faster than Python sort.  Especially for big files.
def sortFile(bedfile):
    sortchr = subprocess.Popen(args=["sort", "-k", "1,1", "-k", "2,2n", "-k", "3,3n"],  stdin=open(bedfile),  stdout=open(bedfile+".sorted", "w"))
    sortchr.wait()
    os.remove(bedfile)
    
def moveFile(bedfile):
    os.rename(bedfile,bedfile+".sorted")

import fcntl
def get_open_fds():
    fds = []
    soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    for fd in range(0, soft):
        try:
            flags = fcntl.fcntl(fd, fcntl.F_GETFD)
        except IOError:
            continue
        fds.append(fd)
    return fds

# our split file approach wants to create a new file for each chromosome,
# if there are a huge number of chromosomes (i.e. in unfinished assemblies)
# this would have resulted in too many files for python to open (default 1024)
# so we keep a cache of the files we have been using so we arent constantly opening and
# closing file descriptors
# caching strategy depends on what the underlying file looks like
# if it's sorted or mostly sorted we want to use LRU cache
# if it's completely random unsorted we want to use LFU cache
def splitFile(bedFileLoc,fh_cache = "lfu"):
    
    assert fh_cache in ["singleton","lfu","lru"]
    
    def closeCSV(f):
        _, fh = f
        fh.close()
    
    # defined in the function so the RLIMIT data is up to date
    # max open files in this = max allowed by python - open files already - 10 for breathing space (debugger attachment etc)
    @cache.lru_cache(maxsize=min(resource.getrlimit(resource.RLIMIT_NOFILE))-10-len(get_open_fds()), destructor=closeCSV)
    def lru_getSplitFileChrm(fileName):
        f = open(fileName,"a")
        return (csv.writer(f,delimiter="\t"),f)
    
    @cache.lfu_cache(maxsize=min(resource.getrlimit(resource.RLIMIT_NOFILE))-10-len(get_open_fds()), destructor=closeCSV)
    def lfu_getSplitFileChrm(fileName):
        f = open(fileName,"a")
        return (csv.writer(f,delimiter="\t"),f)
    
    @cache.lru_cache(maxsize=1, destructor=closeCSV)
    def singleton_getSplitFileChrm(fileName):
        f = open(fileName,"a")
        return (csv.writer(f,delimiter="\t"),f)
    
    with open(bedFileLoc,"r") as bedFile:
        basefilename = basename(bedFileLoc)
    
        infile = csv.reader(bedFile, delimiter="\t")
        for row in infile:
            chrm = row[0]
            
            # different caching is better for different types of splitting the file
            if fh_cache == "singleton": 
                outcsvfile,_ = singleton_getSplitFileChrm("temp-"+chrm+"-"+basefilename)
            elif fh_cache == "lfu":
                outcsvfile,_ = lfu_getSplitFileChrm("temp-"+chrm+"-"+basefilename)
            elif fh_cache == "lru":
                outcsvfile,_ = lru_getSplitFileChrm("temp-"+chrm+"-"+basefilename)
            else:
                assert False, "Unknown fh cache type"
            outcsvfile.writerow(row)
        
        # clear caches and close out the files
        if fh_cache == "singleton": 
            singleton_getSplitFileChrm.clear()
        elif fh_cache == "lfu":
            lfu_getSplitFileChrm.clear()
        elif fh_cache == "lru":
            lru_getSplitFileChrm.clear()
        else:
            assert False, "Unknown fh cache type"

def doFile(bedFileLoc, outfileLoc, withDupes, normalise, cores, fullySorted, chrmSorted):
    
    # split file by chromosome
    splitFile(bedFileLoc,"singleton" if chrmSorted else "lfu")

    chrs = []
    files = os.listdir(os.getcwd())
    for bedfile in files:
        if bedfile.startswith("temp-"):
            chrs.append(bedfile)
    
    # sort each chr (if necessary)
    from multiprocessing import Pool
    p = Pool(cores)
    if fullySorted:
        p.map(moveFile,chrs)
    elif withDupes:
        p.map(sortFile,chrs)
    else:
        p.map(sortAndUniq,chrs)
    
    reads = file_len(bedFileLoc)
    print reads
    
    accumulatedreads = 0

    with open(outfileLoc, "w") as outfile:
        outcsvfile = csv.writer(outfile, delimiter="\t", lineterminator="\n")
        for bedfile in chrs:
            # do each chromosome
            if normalise:
                accumulatedreads += do_chrm(bedfile+".sorted", outcsvfile, reads/1000000.0)
            else:
                accumulatedreads += do_chrm(bedfile+".sorted", outcsvfile, 1.0)
            os.remove(bedfile+".sorted")
            
    assert accumulatedreads == reads, "we somehow skipped some reads"

if __name__ == "__main__":
    try:
        opts, args = getopt.getopt(sys.argv[1:], "i:o:c:dns", ["bed=","bedGraph=","cores=","keepDupes","normalise","fullySorted","chrmSorted"])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)

    withDupes = False
    cores = 4
    normalise = False
    
    fullySorted = False
    chrmSorted = False

    infileLoc = ""
    outfileLoc = ""
    
    for o, a in opts:
        # need a vstep width
        if o in ["-i","--bed"]:
            infileLoc = a
        elif o in ["-o","--bedGraph"]:
            outfileLoc = a
        elif o in ["-c","--cores"]:
            cores = int(a)
        elif o in ["-d","--keepDupes"]:
            withDupes = True
        elif o in ["-n","--normalise"]:
            normalise = True
        elif o in ["--fullySorted"]:
            fullySorted = True
            chrmSorted = True
        elif o in ["--chrmSorted"]:
            chrmSorted = True

    if infileLoc == "" or outfileLoc == "":
        print "missing file argument!"
        sys.exit(2)

    print "Executing:",infileLoc,"to",outfileLoc

    doFile(infileLoc, outfileLoc, withDupes, normalise, cores, fullySorted, chrmSorted)
