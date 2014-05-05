from __future__ import print_function
import os
import sys
import logging
import argparse
from  bam2x import IO,TableIO,DBI
from  bam2x.Annotation import BED6
from bam2x.Tools import _translate_to_meta,merge_bed,gini_coefficient
from math import log
from bam2x.gmm import model_str,fit_two_peaks_EM,bayes_p2
from numpy import array
import time
VERSION="0.2"
def help():
    return "merged tss or tts (Result from bam2ends.py)"

def set_parser(parser):
    #parser.add_argument("-m",type=str,choices=("seq","cDNA","cdna","cds","utr5","utr3"),dest="method")
    parser.add_argument("--gap",type=int,dest="gap",default=10,help="min gap default:%(default)s")
    parser.add_argument("--min",type=int,dest="min_reads_number",default=0,help="reads number threshold, default:%(default)s")
    parser.add_argument("--gene",type=str,dest="gene",default=0,help="known gene bed12 file")
    parser.add_argument("--tts",action="store_true",dest="tts",help="tts default: tss")
    parser.add_argument("--prefix",dest="prefix",help="prefix for bed name",default="meta")
    

def run(args):
    logging.basicConfig(level=logging.INFO)
    def process():
        if len(buff)==1: return 0
        max_score=0.0
        total_score=0.0
        e=[]
        for i in buff:
            total_score+=i.score
            e.append(i.score)
        e=[i/total_score for i in e]
        gini=gini_coefficient(e)
        if total_score < args.min_reads_number:
            return 0
        record={}
        meta=BED6(buff[0].chr,buff[0].start,buff[-1].stop,args.prefix+"."+str(group_id),total_score,buff[0].strand)
        peak=max(buff,key=lambda x:x.score)
        record["peak"]=peak._replace(score=peak.score/total_score)
        record["meta"]=meta._replace(strand=peak.strand)
        record["gini"]=gini
        records.append(record)
        return 1
    
    
    def simple_output():
        print("# formats: bayes_prob_model2, gini, [ region bed, score is total reads], [peak bed , score is proportion ]",file=out)
        for i,x in enumerate(records):
            print("{p2}\t{gini}\t{meta}\t".format(p2=p2[i],meta=x["meta"],gini=x["gini"]),end="",file=out)
            print(x["peak"],file=out)
    def bed12_output():
        print("# formats: bed12 , [R,G,B] are corresponding to [ TTS_GINI_PVALUE*200, TSS_GINI_PALUE*200, PROPORTION_OF_PEAK*200 ]",file=out)
        for i,x in enumerate(records):
            if args.tts:
                g=0
                r=int(p2[i]*200)
            else:
                g=int(p2[i]*200)
                r=0
            b=int(x["gini"]*200)
            if p2[i]>0.5:
                meta=x["meta"]._replace(id=x["meta"].id+".end")
            else:
                meta=x["meta"]
            rgb="{r},{g},{b}".format(r=r,g=g,b=b)
            print("{bed6}\t{thickStart}\t{thickEnd}\t{itemRgb}\t{blockCount}\t{blockSizes}\t{blockStarts}".format(bed6=meta,thickStart=x["peak"].start,thickEnd=x["peak"].end,itemRgb=rgb,blockSizes=x["meta"].stop-x["meta"].start,blockCount=1,blockStarts=0),file=out)

    
    records=[]
    GAP=args.gap
    fin=IO.fopen(args.input,"r")
    out=IO.fopen(args.output,"w")
    iterator=TableIO.parse(fin,"bed6")
    last=iterator.next()
    last_stop=last.stop
    group_id=0
    buff=[last]
    last_chr=last.chr
    for x,i in enumerate(iterator):
        if x%10000==0: logging.info("processing {x} reads".format(x=x));
        if i.chr!=last_chr or i.start-last_stop > GAP:
            group_id+=process()
            buff=[i]
            last_chr=i.chr
            last_stop=i.stop
        else:
            buff.append(i)
            if i.stop>last_stop:
                last_stop=i.stop

    process()
    gini=array([i["gini"] for i in records])
    model=fit_two_peaks_EM(gini)
    p2=bayes_p2(gini,model)
    print("# Date: ",time.asctime(),file=out)
    print("# Program Version ",VERSION,file=out)
    print("# The command line is :",file=out)
    print("#\t"," ".join(sys.argv),file=out)
    print("# learning model:",file=out)
    print("#",model_str(model),file=out)
    #simple_output()
    bed12_output()
        

    
    
    
        

if __name__=="__main__":
    from bam2x.IO import parser_factory
    p=parser_factory(description=help())
    set_parser(p)
    if len(sys.argv)==1:
        print(p.print_help())
        exit(0)
    run(p.parse_args())







