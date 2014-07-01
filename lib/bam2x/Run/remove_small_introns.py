from __future__ import print_function
import os
import sys
import logging
import argparse
from  bam2x import IO,TableIO,DBI
import itertools
def help():
    return "remove small introns, input is bed12 format file. useful for blat psl output"
def set_parser(parser):
    parser.add_argument("-c","--cutoff",type=int,default=10,dest="cutoff",help="remove intron less than Default:%(default)s bp")
    
    
def run(args):
    fin=IO.fopen(args.input,"r")
    out=IO.fopen(args.output,"w")
    for i in TableIO.parse(fin,"bed12"):
        print(remove_small_introns(i,args.cutoff),file=out)
def remove_small_introns(bed,cutoff):
    state=0
    newBlockStarts=[bed.blockStarts[0]]
    newBlockSizes=[bed.blockSizes[0]]
    i=0
    for start,size in itertools.izip(bed.blockStarts[1:],bed.blockSizes[1:]):
        intronLen=start-newBlockStarts[-1]-newBlockSizes[-1]
        if intronLen < cutoff:
            newBlockSizes[-1]=start-newBlockStarts[-1]+size
        else:
            newBlockStarts.append(start)
            newBlockSizes.append(size)
    newBlockCount=len(newBlockStarts)
    return bed._replace(id=bed.id+"_rsi",blockCount=newBlockCount,blockSizes=newBlockSizes,blockStarts=newBlockStarts)
    



    

if __name__=="__main__":
    from bam2x.IO import parser_factory
    p=parser_factory(description=help())
    set_parser(p)
    if len(sys.argv)==1:
        print(p.print_help())
        exit(0)
    run(p.parse_args())







