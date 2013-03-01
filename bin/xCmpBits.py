#!/usr/bin/env python
# Programmer : zhuxp
# Date: 
# Last-modified: 09 Oct 2012 12:51:47

import os,sys,argparse
from xplib.Annotation import *
from xplib import TableIO
from xplib import DBI

'''
this program Compare overlap nucleitide.
The query regions could be in the format below:
    vcf,bed,genebed, [private format: aps] etc.  DEFAULT: bed
The database or data file format could be:
    vcf,bed,genebed, [private format : aps], etc. DEFAULT: bed

The Simple Example:
    xCmpBits.py -i query.bed -a db.bed


'''
from bitarray import bitarray
def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -i [query.bed|query.genebed|query.OddsRatio.out] --input-format [bed|genebed|oddsratiosnp] -a file1.genebed file2.genebed --db_format genebed', epilog='Library dependency : xplib')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.3')
    p.add_argument('-i','--input',dest="input",type=str,default="stdin",help="input file")
    p.add_argument('-I','--input_format',dest="input_format",action="store",default="bed",help="input file format")
    p.add_argument('-o','--output',dest="output",type=str,default="stdout",help="output file")
    p.add_argument('-a','--annotations',dest="db",action="store",help="feature annotation file",type=str)
    p.add_argument('-A','--annotations_format',dest="db_format",action="store",help="annotation files format [bed|genebed|tabix|vcf]",type=str,default="bed")
    p.add_argument('-g',dest="genomesize",action="store",default="/data/zhuxp/Data/hg19.chrom.25.sizes",help="chrom size file")
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

    db_format=args.db_format

    qchrs={}
    pchrs={}
    for x in TableIO.parse(args.genomesize):
        qchrs[x[0]]=bitarray(x[1])
        qchrs[x[0]].setall(False)
        pchrs[x[0]]=bitarray(x[1])
        pchrs[x[0]].setall(False)
    chrs=sorted(qchrs.keys())
    print >>sys.stderr,"Reading Query File",args.input
    for i,x in enumerate(TableIO.parse(args.input,args.input_format)):
        if i%100==0:
            print >>sys.stderr,"reading ",i," entries          \r",
        if qchrs.has_key(x.chr):    
            qchrs[x.chr][x.start:x.stop]=True
    print >>sys.stderr,"Reading database File",args.db
    for i,x in enumerate(TableIO.parse(args.db,args.db_format)):
        if i%100==0:
            print >>sys.stderr,"reading ",i," entries          \r",
        if pchrs.has_key(x.chr):    
            pchrs[x.chr][x.start:x.stop]=True     
    
    
    q=0
    p=0
    qp=0
    genome=0
    print >>out,"# Query File: ",args.input
    print >>out,"# DB File: ",args.db
    print >>out,"chrom\tquery\tdb\toverlap\texpectation\tfoldchange"

    for i in chrs:
        qbits=qchrs[i].count(42)
        pbits=pchrs[i].count(42)
        qpbits=(pchrs[i] & qchrs[i]).count(42)
        chrom=len(pchrs[i])
        print >>out,i,"\t",chrom,"\t",qbits,"\t",pbits,"\t",qpbits,"\t",
        expect=(float(qbits)/chrom) * float(pbits)
        if expect!=0:
            fold=float(qpbits)/expect
            print >>out,"%.2f"%expect,"\t%.2f"%fold
        else:
            fold=None
            print >>out,"%.2f"%expect,"\tNone"
        q+=qbits
        p+=pbits
        qp+=qpbits
        genome+=chrom
    print >>out,"total",genome,"\t",q,"\t",p,"\t",qp,"\t",
    expect=float(q)/genome*float(p)
    fold=qp/expect
    print >>out,"%.2f"%expect,"\t","%.2f"%fold
    qr=float(q)/genome
    print >>out,"# ",args.input,":\t",qr*100,"%"
    pr=float(p)/genome
    print >>out,"# ",args.db,":\t",pr*100,"%"
    qpr=float(qp)/genome
    print >>out,"#  ovarlap :\t",qpr*100,"%"



if  __name__=="__main__":
    Main()






