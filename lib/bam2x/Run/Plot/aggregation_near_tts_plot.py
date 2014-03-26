from __future__ import print_function
import os
import sys
import logging
import argparse
from bam2x import TableIO,IO,DBI
from bam2x.Tools import compatible_with_transcript,_translate,get_flank_region,translate_coordinates,gini_coefficient
import math
import logging
import matplotlib.pyplot as plt
def help():
    return "aggregation plot figure near tts"
def set_parser(parser):
    parser.add_argument("-b","--bam",type=str,dest="bam")
    parser.add_argument('-I','--input_format',dest="format",choices=("bed3","bed6","bed12"),default="bed12",type=str,help="input file format default=%(default)s")
    parser.add_argument("--up",type=int,dest="up",default=500,help="tts upstream bp")
    parser.add_argument("--down",type=int,dest="down",default=500,help="tts downstream bp")
    #parser.add_argument("--strand",type=str,dest="strand",choices=("read1","read2"),default="read1",help="paired end read strand, default: %(default)s")
  
def run(args):
    logging.basicConfig(level=logging.INFO)
    up=args.up
    down=args.down
    bp_num=up+down
    offset=-up
    bam=DBI.init(args.bam,"bam")
    fin=IO.fopen(args.input,"r")
    out=IO.fopen(args.output,"w")
    bin_sum=[0 for i in xrange(bp_num)]
    bin_e=[0.0 for i in xrange(bp_num)]
    bin_dis=[[] for i in xrange(bp_num)]
    for i0,bed in enumerate(TableIO.parse(fin,args.format)):
        bed_bin=[0 for i in xrange(bp_num)]
        tts=bed.tts()
        tts_flank=get_flank_region(tts,up,down)
        for read in bam.query(tts_flank,"bam1",strand="read1"):
            a=translate_coordinates(tts,read)
            #print(a,file=out)
            for e in a.Exons():
                #print(e,file=out)
                start=e.start-offset
                end=e.stop-offset
                if start < 0: start=0
                if end > bp_num: end=bp_num
                for j in xrange(start,end):
                    bed_bin[j]+=1
        for  i in xrange(bp_num):
            bin_sum[i]+=bed_bin[i]
            bin_dis[i].append(bed_bin[i])
    bed_num=i0+1
    for i in xrange(bp_num):
        bin_e[i]=gini_coefficient(bin_dis[i])
    print("pos_to_tts\taggregation_mean\tgini_coefficient",file=out)
    for i in xrange(bp_num):
        print("{bin}\t{aggregation}\t{E}".format(bin=i+offset,aggregation=float(bin_sum[i])/bed_num,E=bin_e[i]),file=out)
    
    
    ax1=plt.subplot(2,1,1)
    plt.ylabel('gini coeffecient')
    plt.plot(range(-up,down),bin_e)
    ax1.set_ylim(0,1)
    ax1.axes.get_xaxis().set_visible(False)
    plt.axvline(x=0,linewidth=1, color='y')
    ax2=plt.subplot(2,1,2)
    plt.plot(range(-up,down),[float(i)/bed_num for i in bin_sum])
    plt.ylabel('mean coverage')
    plt.xlabel('pos to tts (bp)')
    plt.axvline(x=0,linewidth=1, color='y')
    plt.savefig(args.output+".png")








if __name__=="__main__":
    from bam2x.IO import parser_factory
    p=parser_factory(description=help())
    set_parser(p)
    if len(sys.argv)==1:
        print(p.print_help())
        exit(0)
    run(p.parse_args())

