from __future__ import print_function
import os
import sys
import logging
import argparse
from bam2x import TableIO,IO,DBI
from bam2x.Tools import compatible_with_transcript,_translate,get_flank_region,translate_coordinates,gini_coefficient
import math
import logging

def help():
    return "aggregation plot figure near tss or tts"
def set_parser(parser):
    parser.add_argument("-b","--bam",type=str,dest="bam")
    parser.add_argument('-I','--input_format',dest="format",choices=("bed3","bed6","bed12"),default="bed12",type=str,help="input file format default=%(default)s")
    parser.add_argument("--up",type=int,dest="up",default=500,help="upstream bp")
    parser.add_argument("--down",type=int,dest="down",default=500,help="downstream bp")
    parser.add_argument("--tts",dest="tts",action="store_true",default=False,help="if this option is chosen , aggregate the tts near region")
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
        if args.tts:
            pos=bed.tts()
        else:
            pos=bed.tss()
        pos_flank=get_flank_region(pos,up,down)
        for read in bam.query(pos_flank,"bam1",strand="read1"):
            a=translate_coordinates(pos,read)
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
    if args.tts:
        print("pos_to_tts\taggregation_mean\tgini_coefficient",file=out)
    else:
        print("pos_to_tss\taggregation_mean\tgini_coefficient",file=out)
    for i in xrange(bp_num):
        print("{bin}\t{aggregation}\t{E}".format(bin=i+offset,aggregation=float(bin_sum[i])/bed_num,E=bin_e[i]),file=out)
    
    try:
        import matplotlib.pyplot as plt
        import matplotlib
        matplotlib.rcParams.update({'font.size':9})
        ax1=plt.subplot2grid((7,1),(6,0))
        plt.ylabel('gini coeffecient')
        plt.fill_between(range(-up,down),bin_e,color="r",alpha=0.2,y2=0)
        ax1.set_ylim(0,1)
        ax1.set_xlim(-up,down)
        ax1.axes.get_xaxis().set_visible(False)
        plt.axvline(x=0,linewidth=1, color='y')
        ax2=plt.subplot2grid((7,1),(0,0),rowspan=5)
        ax2.set_xlim(-up,down)
        plt.plot(range(-up,down),[float(i)/bed_num for i in bin_sum])
        plt.ylabel('mean coverage')
        if args.tts:
            plt.xlabel('pos to tts (bp)')
        else:
            plt.xlabel('pos to tss (bp)')
        plt.axvline(x=0,linewidth=1, color='y')
        plt.grid(True)
        plt.savefig(args.output+".png")
    except:
        pass








if __name__=="__main__":
    from bam2x.IO import parser_factory
    p=parser_factory(description=help())
    set_parser(p)
    if len(sys.argv)==1:
        print(p.print_help())
        exit(0)
    run(p.parse_args())

