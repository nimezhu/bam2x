#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 28 Sep 2012 18:33:32

import os,sys,argparse
import pysam
import random
from math import log,sqrt
from xplib import TableIO
from xplib import DBI
from xplib.Annotation import *
def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency : pysam')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.2')
    p.add_argument('--bamA',nargs='*',dest="bamA",action="store",default=[],help="add bam to list A")
    p.add_argument('--bamB',nargs='*',dest="bamB",action="store",default=[],help="add bam to list B")
    p.add_argument('--bamlistA',dest="bamlistA",action="store",default=None,type=str,help="add filename in file to list A")
    p.add_argument('--bamlistB',dest="bamlistB",action="store",default=None,type=str,help="add filename in file to list B")
    p.add_argument('-o','--out',dest="out",type=str,action="store",default="stdout",help="output file")
    p.add_argument('--min_coverage',dest="min_coverage",type=int,action="store",default=10,help="not count the chi-square if either bam coverage is less than this value [%(default)i]")
    p.add_argument('--min_minor_cov',dest="min_minor",type=int,action="store",default=5,help="not count the chi-square if minor allele coverage (add coverage in two bam) is less than this value [%(default)i]")
    p.add_argument('-r','--region',dest="region",type=str,action="store",default=None,help="only report region example: chr1:1-1000")
    p.add_argument('-a','--bed',dest="beds",type=str,action="store",default=None,help="beds file")
    p.add_argument('-g','--chromsize',dest="chromsize",type=str,action="store",default="/data/zhuxp/Data/hg19.chrom.25.sizes",help="get chromosome sizes")
   
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()

hNtToNum={'a':0,'A':0,
          'c':1,'C':1,
          'g':2,'G':2,
          't':3,'T':3
         }
Nt=['A','C','G','T']
def parseRegion(region):
    a=region.split(":")
    chrom=a[0]
    (start,stop)=a[1].split("-")
    start=int(start)-1
    stop=int(stop)
    return Bed([chrom,start,stop,".",".","."])
def print_header():
    print >>out,"# Compare The SNP preferences between files:"
    print >>out,"# Group A:"
    for f in args.bamA:
        print >>out,"#",f
    if args.bamlistA:
        f=open(args.bamlistA)
        for line in f:
            line=line.strip()
            print >>out,"#",line
        f.close()
    print >>out,"# Group B:"
    for f in args.bamB:
        print >>out,"#",f
    if args.bamlistB:
        f=open(args.bamlistB)
        for line in f:
            line=line.strip()
            print >>out,"#",line
        f.close()
    print >>out,"# Parameters:"
    print >>out,"# Minimal Coverage:          ",args.min_coverage
    print >>out,"# Minimal Minor SNP Coverage:",args.min_minor
    if args.region:
        print >>out,"# Report Region: ",args.region
    elif args.beds:
        print >>out,"# Report Regions in Bedfile:",args.beds
    elif args.chromsize:
        print >>out,"# Report Chromosome in ChrSizes File:",args.chromsize

def binaryFilter(aps):
    A_coverage=sum(aps.A_nt_dis)
    B_coverage=sum(aps.B_nt_dis)
    minor_coverage=aps.A_nt_dis[hNtToNum[aps.minor_allele]]+aps.B_nt_dis[hNtToNum[aps.minor_allele]]
    if A_coverage < args.min_coverage : return False
    if B_coverage < args.min_coverage : return False
    if minor_coverage < args.min_minor : return False
    return True

    
def Main():
    global args,chrs,lengths,out
    args=ParseArg()
    if args.out=="stdout":
        out=sys.stdout
    else:
        try:
            out=open(args.out,"w")
        except IOError:
            print >>sys.stderr,"can't open file",args.out,"to write, using stdout instead"
            out=sys.stdout

    if args.bamlistA:
        dbi_A=DBI.init(args.bamlistA,"bamlist")
    else:
        dbi_A=DBI.init(args.bamA,"bamlist")
    if args.bamlistB:
        dbi_B=DBI.init(args.bamlistB,"bamlist")
    else:
        dbi_B=DBI.init(args.bamB,"bamlist")

    print_header()
    if args.region:
        i=parseRegion(args.region)
        for aps in QueryBed(i,dbi_A,dbi_B):
                print >>out,aps
    elif args.beds:
        for i in TableIO.parse(args.beds,"bed"):
            for aps in QueryBed(i,dbi_A,dbi_B):
                print >>out,aps
    elif args.chromsize:
        for x in TableIO.parse(args.chromsize):
            (chr,size)=x
            binsize=1000000
            chr=chr.strip()
            for i in xrange(0,size,binsize):
                start=i
                stop=i+binsize
                if stop>size: stop=size
                bed=Bed([chr,start,stop,".",".","."])
                for aps in QueryBed(bed,dbi_A,dbi_B):
                    print >>out,aps


        


def QueryBed(i,dbi_A,dbi_B):
    print >>sys.stderr,"Query Region:",i.chr,i.start,i.stop,"                                    \r",
    chr=i.chr
    offset=i.start
    As=[]
    Bs=[]
    for x in dbi_A.query(i):
        As.append(x)
    for x in dbi_B.query(i):
        Bs.append(x)
    for j in range(i.length()):
        aps=OddsRatioSNP(A=As[j],B=Bs[j],chr=chr,start=offset+j)
        if binaryFilter(aps):
            yield aps
            
    
   
if __name__=="__main__":
    Main()


