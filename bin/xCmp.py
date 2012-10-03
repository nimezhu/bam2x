#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 02 Oct 2012 17:14:41

import os,sys,argparse
from xplib.Annotation import *
from xplib import TableIO
from xplib import DBI

'''
this program Compare overlap features or data in the query regions.
The query regions could be in the format below:
    vcf,bed,genebed, [private format: aps] etc.  DEFAULT: bed
The database or data file format could be:
    tabix,bam,vcf,bed,genebed, [private format : aps], etc. DEFAULT: bed

The Simple Example:
    xCmp.py -i query.bed -a db.bed
    or
    xCmp.py -i query.bed -a db1.bed db2.bed
    or 
    xCmp.py -i query.bed -a db1.bed db2.bed -m
    the result will be table format if -m option is chosen.
    or
    xCmp.py -i query.vcf -I vcf -a db1.vcf db2.vcf -A VCF
    or
    xCmp.py -i query.bed -a file1.bam file2.bam -A bam
    or
    cat query.bed | xCmp.py -a db.bed

Most of these queries are dependent on xplib.DBI

'''
def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -i [query.bed|query.genebed|query.OddsRatio.out] --input-format [bed|genebed|oddsratiosnp] -a file1.genebed file2.genebed --db_format genebed', epilog='Library dependency : xplib')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.3')
    p.add_argument('-i','--input',dest="input",type=str,default="stdin",help="input file")
    p.add_argument('-I','--input_format',dest="input_format",action="store",default="bed",help="input file format")
    p.add_argument('-o','--output',dest="output",type=str,default="stdout",help="output file")
    p.add_argument('-a','--annotations',dest="db",action="store",default=[],help="feature annotation files",nargs="+")
    p.add_argument('-A','--annotations_format',dest="db_format",action="store",default=[],help="annotation files format [bed|genebed|tabix|vcf]",nargs="+")
    p.add_argument('-m',dest="m",action="store_true",default=False,help="print in table format")
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
    if args.input=="stdin":
        input=sys.stdin
    else:
        input=open(args.input,"r")

    db_format=args.db_format
    if len(db_format)==0:
        db_format=["bed" for i in range(len(args.db))]
    if len(db_format)==1:
        db_format=[db_format[0] for i in range(len(args.db))]
    if len(db_format)!=len(args.db):
        print >>sys.stderr,"the number of annotation files is not same with the number of annotation formats"
        print >>sys.stderr,"db format",db_format
        print >>sys.stderr,"db ",args.db
        exit(0)
       
    print >>out,"# Input:",args.input
    dbis=[]
    hits=[]  #count the hits
    hcode={}
    for i,f in enumerate(args.db):
        print >>out,"# Features File No."+str(i+1),":",f
        dbis.append(DBI.init(f,db_format[i]))
        hits.append(0)
    query_num=0
    for bed in TableIO.parse(input,args.input_format):
        if not args.m:
            print >>out,"QR\t",bed
        query_num+=1
        code="@"
        for i,dbi in enumerate(dbis):
            flag=0
            for hit in dbi.query(bed):
                if not args.m:
                    print >>out,"\tDB"+str(i+1)+" HT\t",hit
                flag=1
            hits[i]+=flag
            code+=str(flag)
        if hcode.has_key(code):
            hcode[code]+=1
        else:
            hcode[code]=1
        if not args.m:
            print >>out,"CD "+code,"\t",bed
            print >>out,""
            print >>out,""
        else:
            print >>out,bed,"\t","CD "+code

    for i,x in enumerate(hits):
        print >>out,"#",x,"/",query_num,"overlap with No."+str(i+1),args.db[i]

    for key in sorted(hcode.keys()):
        print >>out,"# code:"+key,"\t",hcode[key]




    
if __name__=="__main__":
    Main()




