#!/usr/bin/env python
# Programmer : zhuxp
# Date: 
# Last-modified: 03-17-2014, 15:15:32 EDT
from __future__ import print_function
import os,sys,argparse
from bam2x import TableIO,Tools
from bam2x import IO
import logging
def set_parser(p):
    ''' This Function Parse the Argument '''
    p.add_argument('-I','--input_format',dest="format",type=str, default="guess", choices=TableIO.hclass.keys(),help="input format")
def help():
    return "sort genome annotation files. input file is genome annotation file."
    
def run(args):
    fin=IO.fopen(args.input,"r")
    out=IO.fopen(args.output,"w")
    if args.format=="guess":
        args.format=IO.guess_format(args.input)
    s=TableIO.parse(args.input,args.format)
    l=[]
    for i,x in enumerate(s):
        if i/10000==0:
            logging.info("reading %s entrys in %s",i,args.input)
        l.append(x)
    logging.info("begin sorting")
    l.sort()
    logging.info("sorting done")
    for i in l:
        print(i,file=out)
    logging.info("completed")

if __name__=="__main__":
    from bam2x.IO import parser_factory
    p=parser_factory(description=help())
    set_parser(p)
    if len(sys.argv)==1:
        print(p.print_help())
        exit(0)
    run(p.parse_args())

 

