import os
import sys
import logging
import argparse
from bam2x import TableIO,IO,DBI
from bam2x.Tools import seq_wrapper
def help():
    return "get sequence from twobit genome file"
def set_parser(parser):
    parser.add_argument("-m",type=str,choices=("seq","cDNA","cdna","cds","utr5","utr3"),dest="method",default="cdna")
    parser.add_argument("-g","--genome",type=str,dest="genome",help="chromsome.2bit file")
    parser.add_argument("-b","--bed_column_number",type=int,choices={3,6,12},dest="bed_column_number",default=12) 
    
def run(args):
    bedformat="bed"+str(args.bed_column_number)
    dbi=DBI.init(args.genome,"genome")
    out=IO.fopen(args.output,"w")
    for i in TableIO.parse(IO.fopen(args.input,"r"),bedformat):
        print >>out,">",i.id+"_"+args.method
        print >>out,seq_wrapper(dbi.query(i,method=args.method))
