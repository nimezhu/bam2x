#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 28 Sep 2012 14:00:36

import os,sys,argparse
from xplib.Annotation import Bed
from xplib.Annotation import OddsRatioSNP
from xplib import TableIO
import pysam
from xplib import DBI
import signal
signal.signal(signal.SIGPIPE,signal.SIG_DFL)
def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = '', epilog='Library dependency : pysam xplib')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('-i','--input',dest="input",type=str,default="stdin",help="input file")
    p.add_argument('-I','--format',dest="input_format",type=str,help="input file format: vcf or oddsratiosnp",default="vcf")
    p.add_argument('-o','--output',dest="output",type=str,default="stdout",help="output file")
    p.add_argument('-H','--HM',dest="HM",default=[],nargs="+",help="Histone Modification BAM files",required=True)
    p.add_argument('-C','--CTRL',dest="CT",default=[],nargs="+",help="Control BAM files",required=True)
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
    HM_dbi=DBI.init(args.HM,"bamlist")
    CT_dbi=DBI.init(args.CT,"bamlist")
    hits=0
    query=0
    if args.input=="stdin":
        input=sys.stdin
    else:
        input=args.input
    print >>out,"# Query VCFs:",args.input
    print >>out,"# CASE BAMs:"
    for i in args.HM:
        print >>out,"#\t",i
    print >>out,"# CONTROL BAMs:"
    for i in args.CT:
        print >>out,"#\t",i

    j=0
    for x in TableIO.parse(input,args.input_format):
        #print >>out,"QR\t",x
        hms=DBI.query(x,HM_dbi)
        for i in hms:
            hm=i
        cts=DBI.query(x,CT_dbi)
        for i in cts:
            ct=i
        aps=OddsRatioSNP(A=hm,B=ct,chr=x.chr,start=x.start)
        print >>out,aps
        j+=1
        if j%1000==0:
            print >>sys.stderr,"query ",j," entries\r",
    
        



    
if __name__=="__main__":
    Main()




