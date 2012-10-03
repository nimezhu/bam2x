#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 10-02-2012, 16:58:24 CDT
'''
This program calculate the Allele Preference Score for Specific SNPs or Specific Regions or whole genome.

Input: bams file or bamlist file
       chromsizefile or bed region file or vcf file or aps file

For whole genome:
    xbams2APS.py --bamA Case.Rep1.bam Case.Rep2.bam--bamB Control.Rep1.bam Control.Rep2.bam -g chrom.sizes -o Case_Control.APS

For a region: 
    xbams2APS.py --bamA Case.Rep1.bam Case.Rep2.bam--bamB Control.Rep1.bam Control.Rep2.bam -r chr1:1-10000 -o Case_Control.APS [-n]
    if choose -n , report all the sites
    if not, report sites that have min_coverage in both case and control and have at least min_minor_cov minor allele in total samples.

For some regions:
    xbams2APS.py --bamA Case.Rep1.bam Case.Rep2.bam--bamB Control.Rep1.bam Control.Rep2.bam -a file.bed -A bed -o Case_Control.APS [-n]
    
    xbams2APS.py --bamA Case.Rep1.bam Case.Rep2.bam--bamB Control.Rep1.bam Control.Rep2.bam -a file.vcf -A vcf -o Case_Control.APS 
    
the output example:
chr1    564514  T/G     2.75045262342   ( 10 1 332 5 )  [0, 0, 1, 10] vs       [4, 2, 5, 332]
the six columns are:

1. chromosome
2. position or start(0-index) : the first base is 0
3. major allele/minor allele
4. APS score
5. ( major in case, minor in case, major in control, minor in control)
6. [A,C,G,T] in case vs [A,C,G,T] in control

'''
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
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.3')
    p.add_argument('--bamA',"-H",nargs='*',dest="bamA",action="store",default=[],help="Case Bam Files ")
    p.add_argument('--bamB',"-C",nargs='*',dest="bamB",action="store",default=[],help="Control Bam Files")
    p.add_argument('--bamlistA',dest="bamlistA",action="store",default=None,type=str,help="case bam filenames file, if choose this option , --bamA are ignored")
    p.add_argument('--bamlistB',dest="bamlistB",action="store",default=None,type=str,help="control bam filenames file, if choose this option, --bamB are ignored")
    p.add_argument('-o','--out',dest="out",type=str,action="store",default="stdout",help="output file")
    p.add_argument('--min_coverage',dest="min_coverage",type=int,action="store",default=10,help="not count the chi-square if either bam coverage is less than this value [%(default)i]")
    p.add_argument('--min_minor_cov',dest="min_minor",type=int,action="store",default=5,help="not count the chi-square if minor allele coverage (add coverage in two bam) is less than this value [%(default)i]")
    p.add_argument('-r','--region',dest="region",type=str,action="store",default=None,help="Only report this region. example: chr1:1-1000. if choose this option, -a and -g are ignored")
    p.add_argument('-a','--annotation',dest="annotations",type=str,action="store",default=None,help="Only Query Regions in This Annotation File, if choose this option, -g are ignored")
    p.add_argument('-A','--annotation_format',dest="annotation_format",type=str,action="store",default="bed",help="Region Annotation File Format {bed,vcf,genebed}")
    p.add_argument('-g','--chromsize',dest="chromsize",type=str,action="store",default="/data/zhuxp/Data/hg19.chrom.25.sizes",help="get chromosome sizes")
    p.add_argument('-n','--no_filter',dest="no_filter",action="store_true",default=False,help="No Filter, Report all the sites, ignore the coverage cutoff")
   
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
    '''
    parse string like "chr1:222-444"
    into 
    Bed(chr="chr1",start=221,stop=444)
    '''
    a=region.split(":")
    chrom=a[0]
    (start,stop)=a[1].split("-")
    start=int(start)-1
    stop=int(stop)
    return Bed([chrom,start,stop,".",".","."])
def print_header():
    '''
    print the output header
    '''
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
    if not args.no_filter:
        print >>out,"# \tMinimal Coverage:          ",args.min_coverage
        print >>out,"# \tMinimal Minor SNP Coverage:",args.min_minor
    if args.region:
        print >>out,"# \tReport Region: ",args.region
    elif args.annotations:
        print >>out,"# \tReport Regions in Annotation File:",args.annotations
    elif args.chromsize:
        print >>out,"# \tReport Chromosome in ChrSizes File:",args.chromsize

def binaryFilter(aps):
    '''
    Filter the sites that have less coverage or have less minor allele coverage.
    '''
    if args.no_filter: return True
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
    '''
    Priority:
    Region > Annotations > chromSize
    '''
    if args.region:
        '''
        Query Only Region
        '''
        i=parseRegion(args.region)
        for aps in QueryBed(i,dbi_A,dbi_B):
                print >>out,aps
    elif args.annotations:
        '''
        Query Regions in Bed file or VCF file etc.
        '''
        for i in TableIO.parse(args.annotations,args.annotation_format):
            for aps in QueryBed(i,dbi_A,dbi_B):
                print >>out,aps
    elif args.chromsize:
        '''
        Query Whole Genome
        Chromsize File Example:
        chr1    249250621
        chr2    243199373
        chr3    198022430
        .
        .
        .

        '''
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
    else:
        print >>sys.stderr," at least one of the options -r,-g,-a are required"
        exit(0)


        


def QueryBed(i,dbi_A,dbi_B):
    '''
    Query the region nt distribution use DBI.query()
    '''
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


