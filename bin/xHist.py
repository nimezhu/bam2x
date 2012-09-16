#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 16 Sep 2012 14:32:45
import numpy
import os,sys,argparse
from xplib import TableIO
def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Draw LogRatio Result Figures', epilog='Library dependency : xplib')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('-i','--input',dest="input",type=str,required=True,help="input file")
    p.add_argument('-L','--chr_length',dest="chr_length_file",type=str,help="chromosome length file",default="/data/zhuxp/Data/hg19.chrom.25.sizes")
    p.add_argument('-o','--output',dest="output",type=str,default="stdout",help="output file")
    p.add_argument('--shiftsize',dest="shiftsize",type=int,default=20,help="bin shiftsize, default: bin = position>>%(default)i")
    p.add_argument('--format',dest="format",type=str,default="oddsratiosnp",help="annotation format default: %(default)s")
    
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(1)
    return p.parse_args()
def parseAnnotationFile(fn):
    format=args.format
    if format=="bam":
        format="bam2bed"
    for x in TableIO.parse(fn,format):
        if not data.has_key(x.chr):
            print >>sys.stderr,"ignore",x
            print >>sys.stderr,"since this chromosome size is not in ",args.chr_length_file
            continue
        bin_start=x.start>>SHIFTSIZE
        bin_stop=x.stop>>SHIFTSIZE
        for bin in range(bin_start,bin_stop+1):
            data[x.chr][bin]+=1


    
def Main():
    global args,out,SHIFTSIZE,data
    args=ParseArg()
    if args.output=="stdout":
        out=sys.stdout
    else:
        try:
            out=open(args.output,"w")
        except IOError:
            print >>sys.stderr,"can't open file ",args.output,"to write. Using stdout instead"
            out=sys.stdout


    SHIFTSIZE=args.shiftsize
    data={}
    for x in TableIO.parse(args.chr_length_file,"simple"):
        data[x[0].strip()]=[0 for row in range((long(x[1]>>SHIFTSIZE)+1))] 
    parseAnnotationFile(args.input)
    for x in TableIO.parse(args.chr_length_file,"simple"):
        chrom=x[0].strip()
        length=long(x[1])
        for i,bin in enumerate(data[chrom]):
            start=i<<SHIFTSIZE
            stop=(i+1)<<SHIFTSIZE
            if stop > length: stop=length
            print >>out,chrom+"\t"+str(start)+"\t"+str(stop)+"\t"+str(bin)
    
if __name__=="__main__":
    Main()


