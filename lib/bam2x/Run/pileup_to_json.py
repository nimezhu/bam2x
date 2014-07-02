#!/usr/bin/env python
# Programmer : zhuxp
# Date: 
# Last-modified: 07-02-2014, 15:35:20 EDT
from __future__ import print_function
VERSION="0.1"
import os,sys,argparse
import json
from bam2x import IO,TableIO
def help():
    return "convert bam2x pileup_compatible_reads output file to json file"  
def set_parser(parser):
    #parser.add_argument("-m",type=str,choices=("seq","cDNA","cdna","cds","utr5","utr3"),dest="method")
    parser.add_argument('-q','--query',dest="query",default="all",type=str,help="query RNA name")
    
def run(args):
    fin=IO.fopen(args.input,"r")
    out=IO.fopen(args.output,"w")
    r=[]
    m=0
    ideograms=[]
    qr=""
    for i in TableIO.parse(fin,sep=","):
        if len(i)==1:
            a=i[0].split("\t")
            if len(a)==2:
                if args.query=="all" or args.query==qr:
                    if m > 0:
                        ideograms.append({"id":qr,"length":m})
                qr=a[1].strip()
        else:
            if args.query=="all" or args.query==qr:
                r.append({"chr":qr,"start":i[0],"length":i[2],"value":i[1]})
                m=int(i[0])+int(i[2])
    if args.query=="all" or args.query==qr:
        if m > 0:
           ideograms.append({"id":qr,"length":m})
    j={
        "ideograms":ideograms,
        "tracks":
        [
         {
             "name":args.input,
             "type":"bedgraph",
             "values":r
         }
        ]
    }
    print(json.dumps(j,indent=4),file=out)



if __name__=="__main__":
    from bam2x.IO import parser_factory
    p=parser_factory(description=help())
    set_parser(p)
    run(p.parse_args())

