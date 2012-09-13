#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 10 Sep 2012 23:08:56

import os,sys,argparse
from xplib.Annotation import Bed
from xplib import TableIO
import pysam
def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -i file.snp -V file.vcf.gz -o output.file', epilog='Library dependency : pysam xplib')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('-i','--input',dest="input",type=str,help="input file")
    p.add_argument('-V','--VCF',dest="vcf",type=str,help="vcf tabix file")
    p.add_argument('-o','--output',dest="output",type=str,default="stdout",help="output file")
    p.add_argument('-1',dest="index1",action='store_true',default=False,help="1 index, default: 0 index")
    p.add_argument('-c','--ichr',dest="ichr",type=int,default=1,help='chrom column number, default: 1')
    p.add_argument('-p','--ipos',dest="ipos",type=int,default=2,help='position column number, default: 2')
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()
def Main():
    global args,out
    args=ParseArg()
    if args.output=="stdout":
        out=sys.stdout
    else:
        try:
            out=open(args.output,"w")
        except IOError:
            print >>sys.stderr,"can't open file ",args.output,"to write. Using stdout instead"
            out=sys.stdout
    ichr=args.ichr-1
    ipos=args.ipos-1
    tabixfile=pysam.Tabixfile(args.vcf)
    index=0
    if args.index1: index=1
    for x in TableIO.parse(args.input,"simple"):
        print x
        chr=x[ichr].strip()
        chr=chr.replace("chr","")
        pos=x[ipos]-index #convert to 0 index
        for vcf in tabixfile.fetch(chr,pos,pos+1):
            print vcf





    
if __name__=="__main__":
    Main()




