import os
import sys
import logging
import argparse
from bam2x import TableIO,IO,DBI
from bam2x.Tools import seq_wrapper
def help():
    return "query bam file"
def set_parser(parser):
    parser.add_argument("-m",type=str,choices=("fetch","pileup"),dest="method",default="fetch")
    parser.add_argument("-b","--bam",type=str,dest="bam",help="bam file")
    parser.add_argument("-c","--bed_column_number",type=int,choices={3,6,12},dest="bed_column_number",default=12) 
    
def run(args):
    bedformat="bed"+str(args.bed_column_number)
    dbi=DBI.init(args.bam,"bam")
    out=IO.fopen(args.output,"w")
    for i in TableIO.parse(IO.fopen(args.input,"r"),bedformat):
        print >>out,"QR",i
        for j in dbi.query(i,method=args.method):
            print >>out,"HT",j
        print >>out,""

