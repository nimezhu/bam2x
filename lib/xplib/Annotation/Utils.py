#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 03 Oct 2012 14:28:26
import os,sys,argparse,types
from xplib.Annotation import Bed 
from xplib import TableIO

'''
obsoleted
please using Struct.binindex instead.


BinIndex Data Structure

Data is hash table. 
Key is Chromosome
Data["chr1"] is a list of feature list
Data["chr"][bin_number] is a list of features that have bin index.

The Most Useful Function is:

Read features into a data structure:
    Example:
        f=TableIO.parse("filename","genebed")
        data=readIntoBinIndex(f)
    Example:
        f=TableIO.parse("filename","vcf")
        data=readIntoBinIndex(f)
    Example:
        data=readIntoBinIndex(bedlist)
    

Query the overlap features in the data structure:
    Usage:
        from xplib.Annotation.Utils import *
        for i in iterOverlapFeature(bed,data):
            print i

This BinIndex Data Structure was wrapped into DBI.init() and DBI.query():
    Usage Example:
        dbi=DBI.init("file.bed","bed")
        x=Bed(chr="chr1",start=1,stop=200)
        for i in dbi.query(x):
            print i
TO DO:
    Now BinIndex only have item inherent from Bed Class, which means it need to have at least three attributes, chr, start and stop. and start is 0-index.  
    For mergeBinIndex()
    it use the fourth attribute, id.
    Need to improve to read any table that contain this information but not not in a pre defined format. 
'''
binOffsets=(512+64+8+1,64+8+1,8+1,1,0)
binFirstShift=17
binNextShift=3
def appendIntoBinIndex(data,bed):
    '''
    data is a BinIndex Data Structure
    bed is an anntation feature. ( bed vcf or genebed etc.)
    append an annotation into data
    '''
    a=bed
    bin=binFromRangeStandard(a.start,a.stop)
    if not data.has_key(a.chr):
        data[a.chr]=[[] for row in range(4096+512+64+8+1)]
    data[a.chr][bin].append(a)

def deleteFromBinIndex(data,bed):
    '''
    delete bed entry from a BinIndex Data Structure data.
    '''
    a=bed
    bin=binFromRangeStandard(a.start,a.stop)
    for i,x in enumerate(data[a.chr][bin]):
        if a.start==x.start and a.stop==x.stop and a.id==x.id:
            del data[a.chr][bin][i]

def mergeBinIndex(data,other):
    '''
    merge two binindex data structure
    and only ratain one if duplicates.
    '''
    for i in iterBinIndex(other):
        flag=1
        for j in iterOverlapFeature(data):
            if j.start==i.start and j.stop==i.stop and j.id==i.id:
                flag=0
        if flag:
            appendIntoBinIndex(data,i)

    
def iterBinIndex(data):
    '''
    iter all the features in binindex structure
    Usase:
        data is a binindex data structure

        for i in iterBinIndex(data):
            print i
    '''
    chrs=data.keys()
    chrs.sort()
    for chr in chrs:
        for binindex in range(4096+512+64+8+1):
            for i in data[chr][binindex]:
                yield i

def readIntoBinIndex(handle):
    '''
    Read list  or iterator into BinIndex Data Structure
    Usage:
    Example:
        f=TableIO.parse("filename","genebed")
        data=readIntoBinIndex(f)
    Example:
        f=TableIO.parse("filename","vcf")
        data=readIntoBinIndex(f)
    Example:
        data=readIntoBinIndex(bedlist)
    '''
    data={}
    
    for i in handle:
        a=i
        if type(i)==type([]) or type(i)==type((1,2,3)):
            a=Bed(i)
        appendIntoBinIndex(data,a)
    return data
    
def binFromRangeStandard(start,end):
    '''
        Convert (start,end) to bin index
        
        Usage:
           data=[[] for row in range(4096+512+64+8+1)]
           for i in list:
               bin=binFromRangeStandard(i.start,i.stop)
               data[bin].append(i)
    '''
    startBin=start
    endBin=end-1 
    startBin >>= binFirstShift
    endBin >>= binFirstShift
    for i in range(len(binOffsets)):
        if startBin==endBin:
            return binOffsets[i]+startBin
        startBin >>= binNextShift
        endBin >>= binNextShift
    print >>sys.stderr,"start %d , end %d out of range"%(start,end) 
    return 0
def iterRangeOverlapBins(start,end):
    '''
        Iterate the possible overlap bin index
        the core function of query overlap features.
        Usage:
            for i in iterRangeOverlapBins(start,end):
                for j in data[i]:
                    if start > j.stop and stop < j.start:
                        print j
    '''
    startBin=start
    endBin=end-1
    startBin >>= binFirstShift
    endBin >>= binFirstShift
    for i in range(len(binOffsets)):
        for j in range(startBin,endBin+1):
            yield j+binOffsets[i]
        startBin >>= binNextShift
        endBin >>= binNextShift

def iterOverlapFeature(bed,data,**kwargs):
    '''
    iterator the bed overlap features in data
    Usage:
        from xplib.Annotation.Utils import *
        for i in iterOverlapFeature(bed,data):
            print i

    **kwargs for further extension for annotation with no pre define format
    '''

    if type(bed)==type((1,2,3)) or type(bed)==([1,2,3]):
        bed=Bed(bed[0:3])  # guess (chrome,start,stop)
    if not data.has_key(bed.chr):
         raise StopIteration 
    D=data[bed.chr]
    for bin in iterRangeOverlapBins(bed.start,bed.stop):
        for f in D[bin]:
            if f.start < bed.stop and f.stop > bed.start:
                yield f
    raise StopIteration

def getOverlapFeatures(bed,data):
    '''
        get a list of features that overlap with bed.
    '''
    a=[]
    for f in iterOverlapFeature(bed,data):
        a.append(f)
    return a


def test():
    if len(sys.argv)==1:
        print >>sys.stderr,"Usage: Utils.py file.bed"
        exit()
    a=TableIO.parse(sys.argv[1],'bed')
    data=readIntoBinIndex(a)
    bed=Bed(["chr1",100000,200000,".",".","."])
    g=getOverlapFeatures(bed,data)
    print "Overlap with",bed
    for i in g:
        print i
    
    
    
if __name__=="__main__":
    test()


