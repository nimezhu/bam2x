#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 02-13-2014, 00:50:47 EST
from matplotlib.backends.backend_pdf import PdfPages
import numpy
import os,sys,argparse
from bam2x import TableIO
import matplotlib
from pylab import *
import pysam 
def help():
    return "draw reads distribution on chromosomes into pdf file"
def set_parser(p):
    ''' This Function Parse the Argument '''
    p.add_argument('-c','--cex',dest="cex",type=float,default=0.5,help="point size default: %(default)f")
    p.add_argument('-b','--binsize',dest="binsize",type=int,default=200,help="binsize default: %(default)i")
    p.add_argument('-l','--logscale',dest="LogScale",action="store_true",default=False,help="binsize default: %(default)s")
def parseIterChrom(fn):
    sam=pysam.Samfile(fn,"rb")
    lengths=sam.lengths
    chrs=sam.references
    for i,chr in enumerate(chrs):
        bins=[0 for row in range(lengths[i]/binsize+1)]
        for a in sam.fetch(chr):
            bins[a.pos/binsize]+=1
        yield chr,bins
    
def run(local_args):
    global args,out,binsize
    args=local_args
    binsize=args.binsize
    pdf=PdfPages(args.output)
    for chrom,bins in parseIterChrom(args.input):
        fig=figure(figsize=(8,6))
        ax=fig.add_subplot(111)
        print >>sys.stderr,"processing",chrom
        pos=[]
        nz_bins=[]
        for j in range(len(bins)):
            if bins[j]==0: continue
            pos.append(j*binsize)
            nz_bins.append(bins[j])
        if args.LogScale:
            plt.plot(pos,numpy.log(numpy.array(nz_bins)),".",markersize=args.cex,color='b')
        else:
            plt.plot(pos,nz_bins,".",markersize=args.cex,color='b')
        title(chrom)
        if chrom=="chr1 " or chrom=="chr1":
            plt.savefig(pdf,format="pdf")
        else:
            pdf.savefig(fig)
        close()
    pdf.close()
    
if __name__=="__main__":
    run()


