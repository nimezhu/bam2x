from __future__ import print_function
import os
import sys
import logging
import argparse
from  bam2x import IO,TableIO,DBI
from  bam2x.Annotation import BED12
import logging

def help():
    return "Processing Query RNASeq Result, and count the number of compatible reads and overlap&not compatible reads. normalized by number of hits. Input is the output of query_RNASeq."
def set_parser(parser):
    #parser.add_argument("-m",type=str,choices=("seq","cDNA","cdna","cds","utr5","utr3"),dest="method")
    
    pass

def run(args):
    logging.basicConfig(level=logging.INFO)
    fin=IO.fopen(args.input,"r")
    out=IO.fopen(args.output,"w")
    print("GENE\tHIT\tOVERLAP",file=out);
    for qr,hits,overlap in iterate(fin):
        print("{qr_name}\t{hits_num:10.2f}\t{overlap_num:10.2f}".format(qr_name=qr.id,hits_num=count(hits),overlap_num=count(overlap)),file=out);
def count(bedlist):
    s=0.0
    for i in bedlist:
        x=float(i.itemRgb.split(",")[0])
        if x==0.0: x=1.0
        s+=1.0/x
    return s


def iterate(fin):
    buf=fin.next();
    x=buf.split("\t")[1:]
    qr=BED12._make(BED12._types(x))
    hits=[]
    overlap=[]
    i=0
    for x in TableIO.parse(fin):
        if x[0]=="QR":
            if i%100==0: 
                logging.info("processing "+str(i)+"  genes");
            i+=1
            yield qr,hits,overlap
            qr=BED12._make(BED12._types(x[1:]))
            hits=[]
            overlap=[]
        elif x[0]=="HT":
            hits.append(BED12._make(BED12._types(x[1:])))
        elif x[0]=="OP":
            overlap.append(BED12._make(BED12._types(x[1:])))
    yield qr,hits,overlap
        

if __name__=="__main__":
    from bam2x.IO import parser_factory
    p=parser_factory(description=help())
    set_parser(p)
    if len(sys.argv)==1:
        print(p.print_help())
        exit(0)
    run(p.parse_args())







