from __future__ import print_function
import os
import sys
import logging
import argparse
from bam2x import TableIO,IO,DBI
from bam2x.Tools import seq_wrapper
def help():
    return "get sequence from twobit genome file"
def set_parser(parser):
    parser.add_argument("-m",type=str,choices=("seq","cDNA","cdna","cds","utr5","utr3"),dest="method",default="seq")
    parser.add_argument("-g","--genome",type=str,dest="genome",help="chromsome.2bit file")
    parser.add_argument("-b","--bed_column_number",type=int,choices=[3,6,12],dest="bed_column_number",default=12) 
    
def run(args):
    bedformat="bed"+str(args.bed_column_number)
    dbi=DBI.init(args.genome,"genome")
    out=IO.fopen(args.output,"w")
    for i in TableIO.parse(IO.fopen(args.input,"r"),bedformat):
        seq=dbi.query(i,method=args.method)
        if(len(seq)==0): continue
        print (">{}".format(i.id+"_"+args.method),file=out)
        print (seq_wrapper(seq),file=out,end="")


if __name__=="__main__":
    from bam2x.IO import parser_factory
    p=parser_factory(description=help())
    set_parser(p)
    if len(sys.argv)==1:
        print(p.print_help())
        exit(0)
    run(p.parse_args())

 
