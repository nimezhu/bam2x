#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 07-17-2012, 16:02:11 CDT
from matplotlib.backends.backend_pdf import PdfPages
import numpy
import os,sys,argparse
from xplib import TableIO
import matplotlib
from pylab import *
import pysam 
def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Draw BAM Result Figures', epilog='Library dependency : pysam')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('-i','--input',dest="input",type=str,required=True,help="input file")
    p.add_argument('-o','--output',dest="output",type=str,default="stdout",help="output pdf file",required=True)
    p.add_argument('-c','--cex',dest="cex",type=float,default=0.5,help="point size default: %(default)f")
    p.add_argument('-b','--binsize',dest="binsize",type=int,default=200,help="binsize default: %(default)i")
    p.add_argument('-l','--logscale',dest="LogScale",action="store_true",help="binsize default: %(default)b")
    
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(1)
    return p.parse_args()
def parseIterChrom(fn):
    sam=pysam.Samfile(fn,"rb")
    lengths=sam.lengths
    chrs=sam.references
    for i,chr in enumerate(chrs):
        bins=[0 for row in range(lengths[i]/binsize+1)]
        for a in sam.fetch(chr):
            bins[a.pos/binsize]+=1
        yield chr,bins
    
def Main():
    global args,out,binsize
    args=ParseArg()
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
    Main()


