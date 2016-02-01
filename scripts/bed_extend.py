#!/usr/bin/env python
# Programmer : zhuxp
# Date: 
# Last-modified: 02-01-2016, 13:53:56 EST
from __future__ import print_function
VERSION="0.1"
import os,sys,argparse
from bam2x.Annotation import BED3,BED6
from bam2x import TableIO,Tools
from bam2x import IO
import time
from collections import defaultdict
def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency : bam2x')
    p.add_argument('-v','--version',action='version',version='%(prog)s '+VERSION)
    p.add_argument('-i','--input',dest="input",default="stdin",type=str,help="input file DEFAULT: STDIN")
    p.add_argument('-I','--input_format',dest="format",default="bed3",type=str,help="input file format")
    p.add_argument('-o','--output',dest="output",type=str,default="stdout",help="output file DEFAULT: STDOUT")
    p.add_argument('-b','--bp',dest="bp",type=int,default=10000,help="extend bp number on both end")
    p.add_argument('-n','--name',dest="name",type=str,default="noname",help="merged peak name")
    p.add_argument('-g','--genome',dest="genomeSizes",type=str,help="genome sizes file")
    return p.parse_args()
'''
extend bed region 
 example: find nearby region.
          combine use with bed_union.py

'''
def Main():
    '''
    IO TEMPLATE
    '''
    global args,out
    args=ParseArg()
    fin=IO.fopen(args.input,"r")
    out=IO.fopen(args.output,"w")
    '''
    END OF IO TEMPLATE 
    '''
    beds=[]
    bp=args.bp
    '''
    TODO add chromosizes into extend consideration.
    '''
    for i in TableIO.parse(fin,args.format):
        start=i.start-bp
        if start < 0:
            start = 0
        end=i.stop+bp
        chr=i.chr
        beds.append(BED3(chr,start,end))
    for i in merge_bed6_set(beds,0,args.name):
        print(str(i),file=out)


def merge_bed6_set(beds,cutoff=0,id="noname"):
    '''
    beds list to merge.
    1.assign beds to chromosome
    2.report the coord in bed format
    '''
    k=0;
    mbed=assign_bed_to_chrom(beds)
    for i in mbed.keys():
        b=mbed[i]
        for j in _merge_bed6(b,cutoff):
            yield BED6(i,j[0],j[1],id+"_"+str(k),cutoff+1,".")
            k+=1
    
def initarray():
    return []
def assign_bed_to_chrom(beds):
    mbed=defaultdict(initarray);
    for b in beds:
        mbed[b.chr].append(b);
    return mbed

def _merge_bed6(beds,cutoff=0):
    '''
    a simple turing state change
    given a bedlist in same chromosome
    report the region > cutoff
    > 0 union
    > 1 intesect
    > 2 ( coverage > 3 region )
    '''
    l=[]
    strand=beds[0].strand
    for i in beds:
        l.append((i.start,-1))
        l.append((i.stop,1))
        if strand!=i.strand: strand="."
    chr=beds[0].chr
    l.sort()
    switch=0
    start=l[0][0]
    end=l[-1][0]
    state=0
    assert i[0][1] > 1
    last_pos=l[0][0]
    state=0
    switch=0
    for i in l[0:]:
        state-=i[1]
        if switch==1 and state==cutoff:
            if last_pos!=i[0]:
                yield last_pos,i[0]
            switch=0
        elif switch==0 and state>cutoff:
            last_pos=i[0]
            switch=1
    assert state==0
    #blockCount=len(blockSizes)
    #return BED12(chr,start,end,id,0.0,strand,start,start,"0,0,0",blockCount,blockSizes,blockStarts) 




    
if __name__=="__main__":
    Main()





