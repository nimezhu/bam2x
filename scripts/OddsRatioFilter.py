#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 12 Sep 2012 14:25:35

import os,sys,argparse
from xplib.Annotation import Bed
from xplib import TableIO
def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Filter High Coverage SNPs in OddsRatioSNP output file Example: %(prog)s -i file.OddsRatio.out -t 1000 -o high.out > low.out', epilog='Library dependency : xplib')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('-i','--input',dest="input",type=str,help="input file")
    p.add_argument('-o','--output',dest="output",type=str,default="stdout",help="output file (the high coverage file) ")
    p.add_argument('-t','--threshold',dest="t",type=int,default=1000,help="threshold default: %(default)i")
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()
def sum(x):
    s=0
    for i in x:
        s+=i
    return s
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
    h=0
    l=0
    for i in TableIO.parse(args.input,"oddsratiosnp"):

        if sum(i.A_nt_dis) > args.t and sum(i.B_nt_dis) > args.t:
            h+=1
            print >>out,i
        else:
            l+=1
            print i
    print >>sys.stderr,"High Coverage:",h
    print >>sys.stderr,"Low Coverage:",l
    
if __name__=="__main__":
    Main()




