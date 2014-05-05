from __future__ import print_function
import os
import sys
import logging
import argparse
import pysam
from bam2x import IO,TableIO
from os.path import splitext
from string import digits
from pysam import tabix_index
def help():
    return "tabix index a bed or vcf file"
def set_parser(parser):
    parser.add_argument("-S","--sorted",dest="sorted",action="store_true",help="if file is sorted")
    parser.add_argument("-I","--format",dest="format",type=str,choices=["bed12","bed6","bed3","vcf"],default="bed12")
    #parser.add_argument("-s","--seq",type=int,dest="seq",default=None)
    #parser.add_argument("-b","--begin",type=int,dest="begin",default=None)
    #parser.add_argument("-e","--end",type=int,dest="end",default=None)
    #parser.add_argument("-Z","--zero_index",action="store_true",dest="zero",default=True)
    
    
def run(args):
    fin=IO.fopen(args.input,"r")
    outfile=args.input
    if not args.sorted:
        l = [ i for i in TableIO.parse(fin,args.format) ]
        l.sort()
        name=splitext(args.input)
        outfile = "{name[0]}.sorted{name[1]}".format(name=name)
        out = IO.fopen(outfile,"w")
        for i in l:
            print(i,file=out)
        out.close()
    format=args.format.translate(None,digits)
    tabix_index(outfile,preset=format)
    

if __name__=="__main__":
    from bam2x.IO import parser_factory
    p=parser_factory(description=help())
    set_parser(p)
    if len(sys.argv)==1:
        print(p.print_help())
        exit(0)
    run(p.parse_args())







