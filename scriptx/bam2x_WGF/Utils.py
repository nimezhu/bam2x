#!/usr/bin/python
# Programmer : zhuxp, <nimezhu@gmail.com>
# Date: 
# Last-modified: 22 Nov 2012 13:10:10

import os,sys,argparse,types,string
from xplib.Annotation import Bed 
from xplib import TableIO


binOffsets=(512+64+8+1,64+8+1,8+1,1,0)
binFirstShift=17
binNextShift=3

def iterBinIndex(data):
    pass

def readIntoBinIndex(handle):
    '''
    Read list  or iterator into BinIndex Data Structure
    Usage:
    Example:
        f=TableIO.parse("filename","genebed")
        data=readIntoBinIndex(f)
    Example:
        data=readIntoBinIndex(bedlist)
    '''
    data={}
    
    for i in handle:
        a=i
        if type(i)==type([]) or type(i)==type((1,2,3)):
            a=Bed(i)
        if not data.has_key(a.chr):
            data[a.chr]=[[] for row in range(4096+512+64+8+1)]
        bin=binFromRangeStandard(a.start,a.stop)
        data[a.chr][bin].append(i)
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

def iterOverlapFeature(bed,data):
    '''
    iterator the bed overlap features in data
    Usage:
        from xplib.Annotation.Utils import *
        for i in iterOverlapFeature(bed,data):
            print i
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

def Classify_Overlap(bed, overlap_gene_list=[]):
    '''
        input: 
            bed of interest; 
            overlap_gene_list is the output of the function 'getOverlapFeatures' above
        output:
            return a dictionary:
            dict_ovlp ={ 'antisense':[], 'intronic':[], 'sense':[] }
        usage:
            gene = readIntoBinIndex( genelist_genePred )
            for bed in bedlist:
                OverlapGene    = getOverlapFeatures(bed, gene)
                Overlap_dict   = Classify_Overlap(bed, OverlapGene)
                overlap_string = ''
                for k,v in Overlap_dict.iteritems():
                    if v:
                        overlap_string += "".join([ str(k+'_'+each)+';' for each in v])
                print bed, overlap_string

    '''
    dict_ovlp = {}
    dict_ovlp['antisense'] = []
    dict_ovlp['intronic']  = []
    dict_ovlp['sense']     = []

    for g in overlap_gene_list:
        if g.strand in ['+', '-'] and g.strand is not bed.strand:
            dict_ovlp['antisense'].append(g.id)
        else:
            flag = True
            for intron in g.Introns():
                if intron.start <= bed.start < bed.stop <= intron.stop:
                    dict_ovlp['intronic'].append(intron.id)
                    flag = False
            if flag:
                dict_ovlp['sense'].append(g.id)
    return dict_ovlp
                
def test():
    if len(sys.argv)==1:
        print >>sys.stderr,"Usage: Utils.py file.bed"
        exit()
    a=TableIO.parse(sys.argv[1],'genebed')
    data=readIntoBinIndex(a)
    bed=Bed( ["chr12", 54380000, 54392000, "HOXC", 0, "+"] )
    g=getOverlapFeatures(bed,data)
    Overlap_dict = Classify_Overlap(bed, g)
    overlap_string = ''
    for k, v in Overlap_dict.iteritems():
        if v:
            overlap_string += "".join([ str(k+'_'+each)+';' for each in v])
    print bed, overlap_string
    
    
if __name__=="__main__":
    test()

