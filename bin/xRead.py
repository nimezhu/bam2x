#!/usr/bin/env python
# Programmer : zhuxp
# Date: 
# Last-modified: 10-24-2013, 18:06:00 EDT
import os,sys,argparse
from xplib import TableIO,Tools
from xplib.Tools import IO
def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency : pysam')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('-i','--input',dest="input",type=str,help="input file")
    p.add_argument('-I','--converter',dest="format",type=str, default="guess", choices=TableIO.FormatToIterator.keys(),help="input format, converter: bed , bam2bed, bam2bed12, vcf ....")
    p.add_argument('-o','--output',dest="output",type=str,default="stdout",help="output file")
    
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()
def Main():
    '''
    This program is a test for TableIO.parse(file.bam,"bam2bed")

    '''
    global args,out
    args=ParseArg()
    fin=IO.fopen(args.input,"r")
    out=IO.fopen(args.output,"w")
    if args.format=="guess":
        args.format=Tools.guess_format(args.input)
    s=TableIO.parse(args.input,args.format)
    for i in s:
        print >>out,i
if __name__=="__main__":
    Main()




