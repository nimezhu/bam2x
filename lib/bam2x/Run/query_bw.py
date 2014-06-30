from __future__ import print_function
import os
import sys
import logging
import argparse
from bam2x import TableIO,IO,DBI
from bam2x.Tools import seq_wrapper
def help():
    return "query bigwig file"
def set_parser(parser):
    parser.add_argument("-m",type=str,choices=("DNA","cDNA"),dest="method",default="DNA",help="method default: %(default)s")
    parser.add_argument("-b","--bigwig",type=str,dest="bw",help="bigwig file")
    parser.add_argument("-I","--format",type=str,choices=("bed3","bed6","bed12","vcf"),dest="format",default="bed12",help="input format : %(defualt)s") 
    
def run(args):
    dbi=DBI.init(args.bw,"bigwig")
    out=IO.fopen(args.output,"w")
    for i in TableIO.parse(IO.fopen(args.input,"r"),args.format):
        ht=[ j for j in dbi.query(i,method=args.method) ]
        print("QR",i,file=out)
        print("HT",ht,file=out)

if __name__=="__main__":
    from bam2x.IO import parser_factory
    p=parser_factory(description=help())
    set_parser(p)
    if len(sys.argv)==1:
        print(p.print_help())
        exit(0)
    run(p.parse_args())

 
