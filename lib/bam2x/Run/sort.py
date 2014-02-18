#!/usr/bin/env python
# Programmer : zhuxp
# Date: 
# Last-modified: 02-15-2014, 18:08:41 EST
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
        print >>out,i
    logging.info("completed")
if __name__=="__main__":
    run()




