#!/usr/bin/env python
# Programmer : zhuxp
# Date: 
# Last-modified: 02-12-2014, 15:49:33 EST
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
        print >>out,i
if __name__=="__main__":
    Main()




