#!/usr/bin/env python
# Programmer : zhuxp
# Date: 
# Last-modified: 03-17-2014, 15:14:48 EDT
from __future__ import print_function
import os,sys,argparse
from bam2x import TableIO,Tools
from bam2x import IO
def set_parser(p):
    ''' This Function Parse the Argument '''
    p.add_argument('-I','--converter',dest="format",type=str, default="guess", choices=TableIO.FormatToIterator.keys(),help="input format")
def help():
    return "read bed or bam files."
    
def run(args):
    fin=IO.fopen(args.input,"r")
    out=IO.fopen(args.output,"w")
    if args.format=="guess":
        args.format=IO.guess_format(args.input)
    s=TableIO.parse(args.input,args.format)
    for i in s:
        print(i,file=out)


if __name__=="__main__":
    from bam2x.IO import parser_factory
    p=parser_factory(description=help())
    set_parser(p)
    if len(sys.argv)==1:
        print(p.print_help())
        exit(0)
    run(p.parse_args())

 

