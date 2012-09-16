#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 16 Jul 2012 15:04:47
from matplotlib.backends.backend_pdf import PdfPages
import numpy
import os,sys,argparse
from xplib import TableIO
import matplotlib
from pylab import *
def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Draw LogRatio Result Figures', epilog='Library dependency : pysam')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('-i','--input',dest="input",type=str,required=True,help="input file")
    p.add_argument('-o','--output',dest="output",type=str,default="stdout",help="output pdf file",required=True)
    p.add_argument('-c','--cex',dest="cex",type=float,default=0.1,help="point size default: %(default)f")
    
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(1)
    return p.parse_args()
def parseIterChrom(fn):
    last_chrom= None
    fin=open(fn)
    positions=[]
    x2s=[]
    coverage=[]
    for x in TableIO.parse(fin,'simple'):
        (chrom,pos,snp,x2,x2_matrix,nt_dist)=x
        b=x2_matrix.replace("( ","")
        b=b.replace(" )","")
        a=b.split(" ")
        x2_matrix=[]
        s=0
        for y in a:
            s+=int(y)
            x2_matrix.append(int(y))

        if (last_chrom==None) or (chrom==last_chrom):
            coverage.append(s)
            positions.append(pos)
            x2s.append(x2)
            last_chrom=chrom
            continue
        yield last_chrom,positions,x2s,coverage
        positions=[]
        x2s=[]
        coverage=[]
        coverage.append(s)
        x2s.append(x2)
        positions.append(pos)
        last_chrom=chrom
    yield last_chrom,positions,x2s,coverage

    
def Main():
    global args,out
    args=ParseArg()
    pdf=PdfPages(args.output)
    for chrom,pos,x2s,coverages in parseIterChrom(args.input):
        highc_pos=[]
        highc_x2s=[]
        lowc_pos=[]
        lowc_x2s=[]
        for i,x in enumerate(coverages):
            if x>4000:
                highc_pos.append(pos[i])
                highc_x2s.append(x2s[i])

            else: 
                lowc_pos.append(pos[i])
                lowc_x2s.append(x2s[i])
        fig=figure(figsize=(8,6))
        ax=fig.add_subplot(211)
        print chrom,len(pos)
        plt.plot(lowc_pos,numpy.log(numpy.array(lowc_x2s)),".",markersize=args.cex,color='b')
        plt.plot(highc_pos,numpy.log(numpy.array(highc_x2s)),".",markersize=args.cex,color='r',alpha=0.5)
        ax=fig.add_subplot(212)
        plt.plot(pos,numpy.log(numpy.array(coverages)),".",markersize=args.cex,color='r')
        #ax.text(0.9,0.2,chrom,transform=ax.transAxes)
        title(chrom)
        if chrom=="chr1 " or chrom=="chr1":
            plt.savefig(pdf,format="pdf")
            
        else:
            pdf.savefig(fig)
        close()
    pdf.close()
    
if __name__=="__main__":
    Main()


