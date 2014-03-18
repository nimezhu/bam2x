from __future__ import print_function
import os
import sys
import logging
import argparse
from bam2x import TableIO,IO,DBI
from bam2x.Tools import seq_wrapper
def help():
    return "query bam file"
def set_parser(parser):
    parser.add_argument("-m",type=str,choices=("fetch","pileup","seq_and_qual"),dest="method",default="fetch")
    parser.add_argument("-b","--bam",type=str,dest="bam",help="bam file")
    parser.add_argument("-c","--bed_column_number",type=int,choices=[3,6,12],dest="bed_column_number",default=12) 
    
def run(args):
    bedformat="bed"+str(args.bed_column_number)
    dbi=DBI.init(args.bam,"bam")
    out=IO.fopen(args.output,"w")
    for i in TableIO.parse(IO.fopen(args.input,"r"),bedformat):
        print("QR",i,file=out)
        for j in dbi.query(i,method=args.method):
            print("HT",j,file=out)
        print("",file=out)

if __name__=="__main__":
    from bam2x.IO import parser_factory
    p=parser_factory(description=help())
    set_parser(p)
    if len(sys.argv)==1:
        print(p.print_help())
        exit(0)
    run(p.parse_args())

 
