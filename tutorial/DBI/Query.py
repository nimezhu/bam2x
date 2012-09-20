#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 20 Sep 2012 00:44:23

import os,sys,argparse
from xplib.Annotation import Bed
from xplib import TableIO
import pysam

        
from xplib import DBI


def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -i file.snp -V file.vcf.gz -o output.file', epilog='Library dependency : pysam xplib')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.2')
    p.add_argument('-i','--input',dest="input",type=str,help="input file")
    p.add_argument('-I','--format',dest="input_format",type=str,help="input file format",default="bed")
    p.add_argument('-A','--dbformat',dest="dbformat",type=str,help="input file database format",default="bed")
    p.add_argument('-o','--output',dest="output",type=str,default="stdout",help="output file")
    p.add_argument('-a','--annotations',dest="db",type=str,default="",required=True,help="query annotation files")
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
    dbi=DBI.init(args.db,args.dbformat)
    for x in TableIO.parse(args.input,args.input_format):
        print "QR\t",x
        for j in DBI.query(x,dbi):
                print "HT\t",j
        x.chr=x.chr.replace("chr","")
        for j in DBI.query(x,dbi):
                print "HT\t",j
        



    
if __name__=="__main__":
    Main()




