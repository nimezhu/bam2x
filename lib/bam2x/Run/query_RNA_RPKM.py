import os
import sys
import logging
import argparse
from bam2x import TableIO,IO,DBI
from bam2x.Tools import seq_wrapper,compatible,overlap,_translate_to_meta,translate,compatible_with_transcript
def help():
    return "query RNA Seq bam file , and report RPKM ( or unique RPKM ) only compatible reads are counted. \n input file is a Bed12 format."
def set_parser(parser):
    parser.add_argument("-b","--bam",type=str,dest="bam",help="bam file")
    parser.add_argument("-s","--strand",type=str,dest="strand",choices=["read1","read2"],default="read2",help="strand of RNA SEQ library [read1 or read2] , default %(default)s")
    parser.add_argument("-u",dest="uniq",action="store_true",default=False,help="only report unique mapped reads")
    
def run(args):
    #logging.basicConfig(level=logging.DEBUG)
    dbi=DBI.init(args.bam,"bam")
    mapped=dbi.mapped
    out=IO.fopen(args.output,"w")
    for i in TableIO.parse(IO.fopen(args.input,"r"),"bed12"):
        print >>out,i.id,"\t",
        s=0
        l=i.cdna_length()
        for j in dbi.query(i,method="bam1",strand=args.strand,uniq=args.uniq):
            if compatible_with_transcript(j,i):
                s+=1
        rpkm=float(s)*(1000000.0/mapped)*(1000.0/float(l))
        print >>out,rpkm

if __name__=="__main__":
    from bam2x.IO import parser_factory
    p=parser_factory(description=help())
    set_parser(p)
    if len(sys.argv)==1:
        print(p.print_help())
        exit(0)
    run(p.parse_args())

 
