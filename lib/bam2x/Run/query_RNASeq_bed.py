from __future__ import print_function
import os
import sys
import logging
import argparse
from bam2x import TableIO,IO,DBI
from bam2x.Annotation import BED12
from bam2x.Tools import seq_wrapper,compatible,overlap,_translate_to_meta,translate,compatible_with_transcript
def help():
    return "query RNA Seq bed file (or in tabix format) , and report it in Transcript coordinates. \n input file is a Bed12 format."
def set_parser(parser):
    parser.add_argument("-b","--bed",type=str,dest="bed",help="bed file")
    parser.add_argument("-m",dest="hit",action="store_true",default=False,help="mininum report, only report comptible reads within transcript")
    
def run(args): 
    if (os.path.isfile(args.bed+".tbi")):
        dbi=DBI.init(args.bed,"tabix",cls=BED12)
    else:
        dbi=DBI.init(args.bed,"binindex",cls=BED12)
    out=IO.fopen(args.output,"w")
    for i in TableIO.parse(IO.fopen(args.input,"r"),"bed12"):
        print ("QR\t",i,file=out)
        for j in dbi.query(i):
            if compatible_with_transcript(j,i):
                print ("HT\t{}".format(_translate_to_meta(i,j)),file=out)
            elif not args.hit:
                print ("OP\t{}".format(j),file=out)
        print ("",file=out)

if __name__=="__main__":
    from bam2x.IO import parser_factory
    p=parser_factory(description=help())
    set_parser(p)
    if len(sys.argv)==1:
        print(p.print_help())
        exit(0)
    run(p.parse_args())

 
