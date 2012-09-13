#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 11 Sep 2012 10:28:11

import os,sys,argparse
from xplib.Annotation import Bed
from xplib import TableIO
import math
def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -i file.oddratio.peak -o file.bed -p cell_line_name', epilog='Library dependency : xplib')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('-i','--input',dest="input",type=str,help="input file")
    p.add_argument('-o','--output',dest="output",type=str,default="stdout",help="output file")
    p.add_argument('-p','--prefix',dest="prefix",type=str,default="",help="prefix to ID")
    
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()
def Main():
    global args,out
    MAX_SCORE=200
    args=ParseArg()
    if args.output=="stdout":
        out=sys.stdout
    else:
        try:
            out=open(args.output,"w")
        except IOError:
            print >>sys.stderr,"can't open file ",args.output,"to write. Using stdout instead"
            out=sys.stdout
    i=0
    x=args.input.split("/")
    name=x[-1]
    name=name.replace(".OddsRatio.peaks","")
    name=name.replace(".LogR.Peaks","")
    name=name.replace(".out","")
    name=args.prefix+name

    for x in TableIO.parse(args.input,"simple"):
        if x[0]=="REGION":
            i+=1
            if x[4]==0:
                score=MAX_SCORE
            else:
                score=-10*math.log(x[4],10)
            if score > MAX_SCORE:
                score=MAX_SCORE
            ID=name+"_ORP_"+str(i)
            print >>out,x[1]+"\t"+str(x[2])+"\t"+str(x[3])+"\t"+ID+"\t",
            print >>out,"%.2f"%score        



    
if __name__=="__main__":
    Main()




