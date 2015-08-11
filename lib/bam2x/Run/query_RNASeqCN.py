from __future__ import print_function
import os
import sys
import logging
import argparse
from bam2x import TableIO,IO,DBI
from bam2x.Tools import seq_wrapper,compatible,overlap,_translate_to_meta,translate,compatible_with_transcript
def help():
    return "query RNA Seq bam file , and report it in Transcript coordinates. \n input file is a Bed12 format."
def set_parser(parser):
    parser.add_argument("-b","--bam",type=str,dest="bam",help="bam file")
    parser.add_argument("-s","--strand",type=str,dest="strand",choices=["read1","read2"],default="read2",help="strand of RNA SEQ library [read1 or read2] , default %(default)s")
    
def run(args):
    #logging.basicConfig(level=logging.DEBUG)
    dbi=DBI.init(args.bam,"bam")
    out=IO.fopen(args.output,"w")
    print("id\tcompatible\toverlap_not_compatible",file=out)
    for i in TableIO.parse(IO.fopen(args.input,"r"),"bed12"):
        ov=0
        cp=0
        #print("QR",i)
        for j in dbi.query(i,method="bam1",strand=args.strand):
            #print("HT",j)
            if compatible_with_transcript(j,i):
                cp+=1
            else:
                ov+=1
        print("{}\t{}\t{}".format(i.id,cp,ov),file=out)

if __name__=="__main__":
    from bam2x.IO import parser_factory
    p=parser_factory(description=help())
    set_parser(p)
    if len(sys.argv)==1:
        print(p.print_help())
        exit(0)
    run(p.parse_args())

 
