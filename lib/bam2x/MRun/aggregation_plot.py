from __future__ import print_function
import os
import sys
import logging
import argparse
from bam2x import TableIO,IO,DBI
from bam2x.Tools import compatible_with_transcript,_translate
import math
import logging
import itertools
import multiprocessing as mp
def help():
    return "aggregation plot, input is bed12 file and bam file, output the aggregation plot number ( and the entropy )"
def set_parser(parser):
    parser.add_argument("-b","--bam",type=str,dest="bam")
    parser.add_argument("-B","--bin_num",type=int,default=100,dest="bin_num")
def run(args):
    logging.basicConfig(level=logging.INFO)
    global bam
    bam=DBI.init(args.bam,"bam")
    fin=IO.fopen(args.input,"r")
    out=IO.fopen(args.output,"w")
    bin_sum=[0 for i in xrange(args.bin_num)]
    bin_e=[0.0 for i in xrange(args.bin_num)]
    #neg_bin_sum=[0 for i in xrange(args.bin_num)]
    bin_dis=[[] for i in xrange(args.bin_num)]
    #neg_bin_dis=[[] for i in xrange(args.bin_num)]
    p=mp.Pool(processes=args.num_cpus)
    beds_list=[[] for i in xrange(args.num_cpus)]
    for i0,bed in enumerate(TableIO.parse(fin,"bed12")):
        beds_list[i0%args.num_cpus].append(bed)
    results = p.map(count_list_star,itertools.izip(beds_list,itertools.repeat(args.bin_num)))
    
    bin_sum=[0 for i in xrange(args.bin_num)]
    bin_dis=[[] for i in xrange(args.bin_num)]
    for result in results:
        for i in xrange(args.bin_num):
            bin_sum[i]+=result[0][i]
            bin_dis[i].extend(result[1][i])
    for i in xrange(args.bin_num):
        bin_e[i]=dis2entropy(bin_dis[i])
    print("bin\taggregation\tentropy",file=out)
    for i in xrange(args.bin_num):
        print("{bin}\t{aggregation}\t{E}".format(bin=i,aggregation=bin_sum[i],E=bin_e[i]),file=out)




    

def count_list_star(a_b):
    return count_list(*a_b)
def count_list(beds,bin_num):        
    bin_sum=[0 for i in xrange(bin_num)]
    bin_dis=[[] for i in xrange(bin_num)]
    for bed in beds:
        bed_bin=[0 for i in xrange(bin_num)]
        for read in bam.query(bed,"bam1"):
            cdna_length=bed.cdna_length()
            if compatible_with_transcript(read,bed):
                gene_bed=_translate(bed,read)
                bin_start=gene_bed.start*bin_num/cdna_length
                bin_stop=gene_bed.stop*bin_num/cdna_length
                for i in xrange(bin_start,bin_stop):
                    bed_bin[i]+=1
        for i in xrange(bin_num):
            bin_sum[i]+=bed_bin[i]
            bin_dis[i].append(bed_bin[i])
    return bin_sum,bin_dis

    

def dis2entropy(iterator):
    s=0
    h={}
    for i in iterator:
        if h.has_key(i):
            h[i]+=1
        else:
            h[i]=1
        s+=1
    e=0.0
    for i in h.values():
        f=float(i)/s
        if f!=0.0:
            e-=f*math.log(f)
    return e


if __name__=="__main__":
    from bam2x.IO import mparser_factory
    p=mparser_factory(description=help())
    set_parser(p)
    if len(sys.argv)==1:
        print(p.print_help())
        exit(0)
    run(p.parse_args())

