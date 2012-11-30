#!/usr/bin/env python
# Programmer : zhuxp
# Date: 
# Last-modified: 30 Nov 2012 13:59:10
'''
xQuery.py is an example program for using xplib.DBI interface
it reports the overlap features or data from the query region.

the query file format can be:
bed,vcf,genebed etc.

the database or data file can be:
bam,tabix,vcf,bed,genebed etc.

for tabix ant other genome annotation file ,
    it yield the overlap annotations in this region
for bam file
    it yield the Nucleotides Distribution in each site of this region.
Example:
    xQuery.py -i file.bed -a file.bam
'''
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

    query_length=0
    hits_number=0
    for x in TableIO.parse(input,args.input_format):
        print >>out,"QR\t",x
        hit=0
        query+=1
        query_length+=len(x)
        for j in dbi.query(x):
            print >>out,"HT\t",j
            hit=1
            hits_number+=1

        if args.dbformat=="tabix":
            x.chr=x.chr.replace("chr","")
            for j in dbi.query(x):
                print >>out,"HT\t",j
                hit=1
                hits_number+=1
        hits+=hit
    print >>out,"# Query Number:",query,"\n# Query Have Hits:",hits
    print >>out,"# Query Length:",query_length
    print >>out,"# Hits Number:",hits_number
        



    
if __name__=="__main__":
    Main()




