#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 28 Sep 2012 13:51:09

import os,sys,argparse
from xplib.Annotation import Bed
from xplib import TableIO
import pysam
from xplib import DBI
import signal
signal.signal(signal.SIGPIPE,signal.SIG_DFL)

def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -i file.snp -a file.vcf.gz -A tabix -o output.file', epilog='Library dependency : pysam xplib')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.2')
    p.add_argument('-i','--input',dest="input",type=str,default="stdin",help="input file")
    p.add_argument('-I','--format',dest="input_format",type=str,help="input file format",default="bed")
    p.add_argument('-A','--dbformat',dest="dbformat",type=str,help="input file database format. {bed|genebed|tabix|bam}",default="bed")
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
    hits=0
    query=0
    if args.input=="stdin":
        input=sys.stdin
    else:
        input=args.input

    for x in TableIO.parse(input,args.input_format):
        print "QR\t",x
        hit=0
        query+=1
        for j in DBI.query(x,dbi):
            print "HT\t",j
            hit=1
        if args.dbformat=="tabix":
            x.chr=x.chr.replace("chr","")
            for j in DBI.query(x,dbi):
                print "HT\t",j
                hit=1
        hits+=hit
    print >>out,"# Query:",query,"\n# Hits:",hits
        



    
if __name__=="__main__":
    Main()




