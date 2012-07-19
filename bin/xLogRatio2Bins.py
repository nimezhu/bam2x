#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 18 Jul 2012 14:12:09
from matplotlib.backends.backend_pdf import PdfPages
import numpy
import os,sys,argparse
from xplib import TableIO
import matplotlib
from pylab import *
from scipy.stats import chi2

def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Call Peaks from LogRatio Result', epilog='Library dependency : pysam')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('-i','--input',dest="input",type=str,required=True,help="input file")
    p.add_argument('-o','--output',dest="output",type=str,default="stdout",help="output file")
    p.add_argument('-b','--binsize',dest="binsize",type=int,default=200,help="binsize")
    p.add_argument('-p','--pvalue',dest="pvalue",type=float,default=1e-05,help="pvalue threshold")
    p.add_argument('--reads',dest="reads",type=int,default=5000,help="reads threshold")
    
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(1)
    return p.parse_args()
def parseIterChrom(fn):
    last_chrom= None
    fin=open(fn)
    positions=[]
    x2s=[]
    matrix_x2=[]
    for x in TableIO.parse(fin,'simple'):
        (chrom,pos,snp,x2,x2_matrix,nt_dist)=x
        b=x2_matrix.replace("( ","")
        b=b.replace(" )","")
        a=b.split(" ")
        x2_matrix=[]
        for y in a:
            x2_matrix.append(int(y))
        if (last_chrom==None) or (chrom==last_chrom):
            matrix_x2.append(x2_matrix)
            positions.append(pos)
            x2s.append(x2)
            last_chrom=chrom
            continue
        yield last_chrom,positions,x2s,matrix_x2
        positions=[]
        x2s=[]
        matrix_x2=[]
        matrix_x2.append(x2_matrix)
        x2s.append(x2)
        positions.append(pos)
        last_chrom=chrom
    yield last_chrom,positions,x2s,matrix_x2

    
def Main():
    global args,out
    args=ParseArg()
    binsize=args.binsize
    for chrom,pos,x2s,matrix_x2 in parseIterChrom(args.input):
        max_pos=0
        for i in pos:
            if max_pos<i: max_pos=i
        bins=[[] for row in range(max_pos/binsize+1)]
        for i,x in enumerate(pos):
            # if filter:
            bins[x/binsize].append(i)
        for i,b in enumerate(bins):
            if len(b)==0: continue
            SqOR=0.0
            num=0
            T_CV=0
            entries=[]
            for j in b:
                a=matrix_x2[j]
                MR=float(a[1]+a[3]-2)/(a[1]+a[0]-4+a[2]+a[3])  #minor allele ratio
                CASE_CV=a[1]+a[0]-2
                CONTROL_CV=a[2]+a[3]-2
                CV=CONTROL_CV+CASE_CV

                if CONTROL_CV > args.reads: continue

                if MR<0.05: continue
                num+=1
                T_CV+=CV
                SqOR+=x2s[j]
                entries.append((pos[j],matrix_x2[j]))
            if num==0: continue
            pvalue=1.0-chi2.cdf(SqOR,num)
            if pvalue<=args.pvalue:
                print chrom,i*binsize,(i+1)*binsize,
                print pvalue,"\t(",SqOR,",",num,")","\t",float(T_CV)/num
                for x in entries:
                    print " SNP",x
                print 
            




    
if __name__=="__main__":
    Main()


