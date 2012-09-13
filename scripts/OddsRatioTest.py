#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 12 Sep 2012 10:35:37

import os,sys,argparse
from xplib.Annotation import Bed
from xplib import TableIO
def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency : pysam')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('-i','--input',dest="input",type=str,help="input file")
    p.add_argument('-o','--output',dest="output",type=str,default="stdout",help="output file")
    
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
    for i in TableIO.parse(args.input,"oddsratiosnp"):
        print i


    
if __name__=="__main__":
    Main()




