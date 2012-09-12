#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 10 Aug 2012 15:10:52

import os,sys,argparse
from xplib.Annotation import Bed
from xplib import TableIO
def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -i file.bed -o file_extend.bed -e 1000')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('-i','--input',dest="input",type=str,help="input file")
    p.add_argument('-o','--output',dest="output",type=str,default="stdout",help="output file")
    p.add_argument('-e','--extend',dest="extend",type=int,default=5000,help="extend bp")
    
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
    extend=args.extend
    f=open(args.input)
    for bed in TableIO.parse(f,"bed"):
        bed.start-=extend
        bed.stop+=extend
        bed.id+=str("_extend_"+str(extend)+"bp")
        print >>out,bed


    
if __name__=="__main__":
    Main()



