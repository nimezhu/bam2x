#!/usr/bin/env python
# Programmer : zhuxp
# Date: 
# Last-modified: 07-08-2014, 10:28:41 EDT
from __future__ import print_function
import os,sys,argparse
from bam2x import TableIO,Tools
from bam2x import IO
def set_parser(p):
    ''' This Function Parse the Argument '''
    p.add_argument('-s','--strand',dest="strand",type=str, default="read1", choices=("read1","read2"),help="which read represent the RNA strand")
def help():
    return "read bam file and only report the spliced reads."
    
def run(args):
    out=IO.fopen(args.output,"w")
    s=TableIO.parse(args.input,"bam2bed12",strand=args.strand)
    for i in s:
        if i.blockCount>1:
            print(i,file=out)



if __name__=="__main__":
    from bam2x.IO import parser_factory
    p=parser_factory(description=help())
    set_parser(p)
    if len(sys.argv)==1:
        print(p.print_help())
        exit(0)
    run(p.parse_args())

 

