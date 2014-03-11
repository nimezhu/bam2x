from __future__ import print_function
import os
import sys
import logging
import argparse
from bam2x import TableIO,IO,DBI
from bam2x.Tools import compatible_with_transcript,_translate,reverse_strand,translate_coordinates
import math
import logging
import itertools
import multiprocessing as mp
def help():
    return "aggregation plot, input is bed12 file and bam file, output the aggregation plot number ( and the entropy )"
def set_parser(parser):
    parser.add_argument("-b","--bam",type=str,dest="bam")
    parser.add_argument("-B","--bin_num",type=int,default=100,dest="bin_num")
    parser.add_argument("--bp",type=int,default=100,dest="bp",help="upstream and downstream bp number")
    parser.add_argument("--strand",type=str,choices=["read1","read2"],default="read2",dest="strand",help="strand: read1 or read2")

def output(results,bin_num,gene_num,prefix=""):
    bin_sum=[0 for i in xrange(bin_num)]
    bin_e=[0.0 for i in xrange(bin_num)]
    neg_bin_e=[0.0 for i in xrange(bin_num)]
    neg_bin_sum=[0 for i in xrange(bin_num)]
    bin_dis=[[] for i in xrange(bin_num)]
    neg_bin_dis=[[] for i in xrange(bin_num)]
    bin_dis=[[] for i in xrange(bin_num)]
    for result in results:
        for i in xrange(bin_num):
            bin_sum[i]+=result[0][i]
            bin_dis[i].extend(result[1][i])
            neg_bin_sum[i]+=result[2][i]
            neg_bin_dis[i].extend(result[3][i])
    for i in xrange(bin_num):
        bin_e[i]=dis2entropy(bin_dis[i])
        neg_bin_e[i]=dis2entropy(neg_bin_dis[i])
    for i in xrange(bin_num):
        print("{bin}\t{mean}\t{E}\t{neg_mean}\t{neg_E}".format(bin=prefix+str(i),mean=float(bin_sum[i])/gene_num,E=bin_e[i],neg_mean=float(neg_bin_sum[i])/gene_num,neg_E=neg_bin_e[i]),file=out)

def run(args):
    logging.basicConfig(level=logging.INFO)
    global bam,out
    bam=DBI.init(args.bam,"bam")
    fin=IO.fopen(args.input,"r")
    out=IO.fopen(args.output,"w")
    p=mp.Pool(processes=args.num_cpus)
    beds_list=[[] for i in xrange(args.num_cpus)]
    for i0,bed in enumerate(TableIO.parse(fin,"bed12")):
        beds_list[i0%args.num_cpus].append(bed)
    gene_num=i0+1
    print("bin_id\tmean\tentropy\treverse_strand_mean\treverse_strand_entropy",file=out)
    up_results=p.map(count_flank_star,itertools.izip(beds_list,itertools.repeat(args.bp),itertools.repeat(args.strand),itertools.repeat(True)))
    output(up_results,args.bp,gene_num,"UP")
    results = p.map(count_list_star,itertools.izip(beds_list,itertools.repeat(args.bin_num),itertools.repeat(args.strand)))
    output(results,args.bin_num,gene_num,"TR")
    down_results=p.map(count_flank_star,itertools.izip(beds_list,itertools.repeat(args.bp),itertools.repeat(args.strand),itertools.repeat(False)))
    output(down_results,args.bp,gene_num,"DN")

    

def count_list_star(a_b_c):
    return count_list(*a_b_c)
def count_list(beds,bin_num,strand):        
    bin_sum=[0.0 for i in xrange(bin_num)]
    neg_bin_sum=[0.0 for i in xrange(bin_num)]
    bin_dis=[[] for i in xrange(bin_num)]
    neg_bin_dis=[[] for i in xrange(bin_num)]
    for bed in beds:
        neg_bed=bed._replace(strand=reverse_strand(bed.strand))
        bed_bin=[0.0 for i in xrange(bin_num)]
        neg_bed_bin=[0.0 for i in xrange(bin_num)]
        for read in bam.query(bed,"bam1",strand=strand):
            cdna_length=bed.cdna_length()
            bin_step=float(cdna_length)/bin_num
            if compatible_with_transcript(read,bed):
                gene_bed=_translate(bed,read)
                bin_start=gene_bed.start*bin_num/cdna_length
                bin_stop=(gene_bed.stop-1)*bin_num/cdna_length
                if bin_stop==bin_start:
                    bed_bin[bin_start]+=float(gene_bed.stop-gene_bed.start)/bin_step
                else:
                    bed_bin[bin_start]+=bin_start+1-float(gene_bed.start)/bin_step
                    bed_bin[bin_stop]+=float(gene_bed.stop)/bin_step-bin_stop
                    for i in xrange(bin_start+1,bin_stop):
                        bed_bin[i]+=1.0  
            elif compatible_with_transcript(read,neg_bed):
                gene_bed=_translate(bed,read)
                bin_start=gene_bed.start*bin_num/cdna_length
                bin_stop=(gene_bed.stop-1)*bin_num/cdna_length
                if bin_stop==bin_start:
                    neg_bed_bin[bin_start]+=float(gene_bed.stop-gene_bed.start)/bin_step
                else:
                    neg_bed_bin[bin_start]+=bin_start+1-float(gene_bed.start)/bin_step
                    neg_bed_bin[bin_stop]+=float(gene_bed.stop)/bin_step-bin_stop
                    for i in xrange(bin_start+1,bin_stop):
                        neg_bed_bin[i]+=1.0
        for i in xrange(bin_num):
            bin_sum[i]+=bed_bin[i]
            neg_bin_sum[i]+=neg_bed_bin[i]
            bin_dis[i].append(int(bed_bin[i]))
            neg_bin_dis[i].append(int(neg_bed_bin[i]))
    return bin_sum,bin_dis,neg_bin_sum,neg_bin_dis

def count_flank_star(l):
    return count_flank(*l)
def count_flank(beds,bp,strand,upstream=True):
    '''
    flank seq aggregation ( no splicing)
    if upstream is False , count downstream
    '''
    pos_sum=[0 for i in xrange(bp)]
    pos_dis=[[] for i in xrange(bp)]
    neg_sum=[0 for i in xrange(bp)]
    neg_dis=[[] for i in xrange(bp)]
    for bed in beds:
        if upstream:
            flank=bed.upstream(bp)
            offset=bp-flank.cdna_length()  # in case upstream is less than bp 
        else:
            flank=bed.downstream(bp)
            offset=0
        pos=[0 for i in xrange(bp)]
        neg=[0 for i in xrange(bp)]
        for read in bam.query(flank,"bam1",strand=strand):
            translated_read=translate_coordinates(flank,read)
            if translated_read.strand=="+" or translated_read.strand==".":
                for start,size in itertools.izip(translated_read.blockStarts,translated_read.blockSizes):
                    for j in xrange(start+offset,start+size+offset):
                        if j>=0 and j<bp:
                            pos[j]+=1
            else:
                for start,size in itertools.izip(translated_read.blockStarts,translated_read.blockSizes):
                    for j in xrange(start+offset,start+size+offset):
                        if j>=0 and j<bp:
                            neg[j]+=1
        for i in xrange(bp):
            pos_sum[i]+=pos[i]
            neg_sum[i]+=neg[i]
            pos_dis[i].append(pos[i])
            neg_dis[i].append(neg[i])
    return pos_sum,pos_dis,neg_sum,neg_dis

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

