#!/usr/bin/env python
# Programmer : zhuxp
# Date: 
# Last-modified: 01-31-2016, 10:10:50 EST
from __future__ import print_function
VERSION="0.1"
import os,sys,argparse
from bam2x.Annotation import BED6
from bam2x import TableIO,Tools
from bam2x import IO
import time
from collections import defaultdict
def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency : bam2x')
    p.add_argument('-v','--version',action='version',version='%(prog)s '+VERSION)
    p.add_argument('-i','--input',dest="input",nargs="*",help="input bed6 or bed5(macs output) files")
    p.add_argument('-o','--output',dest="output",type=str,default="stdout",help="output file DEFAULT: STDOUT")
    p.add_argument('-c','--cutoff',dest="cutoff",type=int,default=0,help="cutoff( 0:union ) (1: intesect)")
    p.add_argument('-n','--name',dest="name",type=str,default="noname",help="region name DEFAULT: %(default)s")
    if len(sys.argv)==1:
        print(p.print_help(),file=sys.stderr)
        exit(0)
    return p.parse_args()
    '''
    bed union merge bed files and report the union or intersect region (generate the new bed region) .
    '''
def Main():
    '''
    IO TEMPLATE
    '''
    global args,out
    args=ParseArg()
    out=IO.fopen(args.output,"w")
    beds=[]
    for f in args.input:
        for a in TableIO.parse(IO.fopen(f),"bed3"):
            beds.append(a)
    '''
    a=BED6("chr01",100,200,"id",0,".")
    b=BED6("chr03",150,250,"id",0,".")
    b2=BED6("chr02",150,250,"id",0,".")
    b3=BED6("chr01",150,250,"id",0,".")
    c=BED6("chr01",255,300,"id",0,".")
    d=BED6("chr01",275,350,"id",0,".")
    beds=[a,b,c,d,b2,b3]
    '''
    for i in merge_bed6_set(beds,args.cutoff,args.name):
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





