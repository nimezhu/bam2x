from __future__ import print_function
import os
import sys
import logging
import argparse
from  bam2x import IO,TableIO,DBI
from  bam2x.Annotation import BED12
import logging

def help():
    return "Processing Query RNASeq Result, and count the number of compatible reads and overlap&not compatible reads. normalized by number of mapping hits. Input is the output of query_RNASeq."
def set_parser(parser):
    #parser.add_argument("-m",type=str,choices=("seq","cDNA","cdna","cds","utr5","utr3"),dest="method")
    pass

def run(args):
    logging.basicConfig(level=logging.INFO)
    fin=IO.fopen(args.input,"r")
    out=IO.fopen(args.output,"w")
    for qr,hits,overlap in iterate(fin):
        l=qr.cdna_length()
        pileup=[0.0 for i in xrange(l)]
        for i in hits:
            for j in xrange(i.start,i.stop):
                pileup[j]+=nh(i)
        print("QR\t{id}".format(id=qr.id),file=out)
        print("PILEUP\n{value}".format(value=rep(pileup)),file=out)
def rep(list):
    s=""
    last=list[0]
    step=0
    offset=0
    for i in list:
        if i!=last:
            s+="{offset},{value},{step}\n".format(offset=offset,value=last,step=step)
            last=i
            offset+=step
            step=1
        else:
            step+=1
    s+="{offset},{value},{step}".format(offset=offset,value=last,step=step)
    return s
            
def nh(bed): 
    '''
    Return 1.0/NH ( number of hits )  
    '''
    x=float(bed.itemRgb.split(",")[0])
    if x==0.0: x=1.0
    return 1.0/x

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
    run(p.parse_args())







