from __future__ import print_function
import os
import sys
import logging
import argparse
from bam2x import TableIO,IO,DBI
from bam2x.Tools import seq_wrapper
from bam2x.TableIO import hclass
def help():
    return "query tabix file"
def set_parser(parser):
    parser.add_argument("-t","--tabix",type=str,dest="db",help="tabix db file")
    parser.add_argument("-A","--annotation_type",type=str,dest="type",choices=hclass.keys(),default=None,help="tabix entry type , default is None")
    parser.add_argument("-I","--format",type=str,choices=hclass.keys(),dest="format",default="bed12",help="input format : %(default)s") 
    
def run(args):
    out=IO.fopen(args.output,"w")
    cls=None
    if hclass.has_key(args.type):
        cls=hclass[args.type]
        dbi=DBI.init(args.db,"tabix",cls=cls)
    else:
        dbi=DBI.init(args.db,"tabix")
    for i in TableIO.parse(IO.fopen(args.input,"r"),args.format):
        print("QR",i,file=out)
        
        for j,ht in enumerate(dbi.query(i)):
            print("HT_{k}\t{ht}".format(k=j+1,ht=ht),file=out)

if __name__=="__main__":
    from bam2x.IO import parser_factory
    p=parser_factory(description=help())
    set_parser(p)
    if len(sys.argv)==1:
        print(p.print_help())
        exit(0)
    run(p.parse_args())

 
