#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 17 Sep 2012 21:28:51
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
    p.add_argument('-p','--pvalue',dest="pvalue",type=float,default=1e-05,help="pvalue threshold")
    p.add_argument('--reads',dest="reads",type=int,default=10000,help="reads threshold (if reads number is bigger than this number, don't count. (PCR bias)")
    p.add_argument('--gap',dest="gap",type=int,default=1000000,help="gap threshold")
    p.add_argument('-f','--fast',dest="fast",action="store_true",default=False,help="don't calculate big regions")
    
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(1)
    return p.parse_args()
def parseIterRegion(fn):
    '''
    yield each region ( gap < args.gap)
    
    '''
    last_chrom= None
    last_position = 0
    fin=open(fn)
    positions=[]
    x2s=[]
    snps=[]
    matrix_x2=[]
    for x in TableIO.parse(fin,'simple'):
        (chrom,pos,snp,x2,x2_matrix,nt_dist)=x
        '''
        add filter here
        '''
        
        b=x2_matrix.replace("( ","")
        b=b.replace(" )","")
        a=b.split(" ")
        x2_matrix=[]
        for y in a:
            x2_matrix.append(int(y))
        t=x2_matrix

        CV=int(t[0])+int(t[1])+int(t[2])+int(t[3])
        MR=(float(t[1])+float(t[3]))/CV
        if(CV>args.reads or MR<0.05):continue

        # end of filter
        if (last_chrom==None) or (chrom==last_chrom and pos-last_position < args.gap):
            matrix_x2.append(x2_matrix)
            positions.append(pos)
            x2s.append(x2)
            snps.append(snp)
            last_chrom=chrom
            last_position=pos
            continue
        yield last_chrom,positions,x2s,matrix_x2,snps
        positions=[]
        x2s=[]
        matrix_x2=[]
        snps=[]
        matrix_x2.append(x2_matrix)
        x2s.append(x2)
        positions.append(pos)
        snps.append(snp)
        last_chrom=chrom
        last_position=pos
    yield last_chrom,positions,x2s,matrix_x2,snps


def Segment(p,offset=0,depth=0):
    '''
    recursively report the regions have lowest pvalue

    '''
    mi=0
    mj=0
    n=len(p)
    for i in range(n):
        for j in range(i,n):
            if p[i][j]<p[mi][mj]:
                mi=i
                mj=j
            if p[i][j]==p[mi][mj] and (j-i)>(mj-mi):
                mi=i
                mj=j
 #   print "Segment( depth:",depth,"offset:",offset,") yielding"
 #   print "yield (",mi+offset,mj+offset,")"
 #   print p
    yield (mi+offset,mj+offset)
    up_p=p[0:mi,0:mi]
    if len(up_p > 0):
        for i in Segment(up_p,offset,depth+1):
            yield i
    mj1=mj+1
    down_p=p[mj1:n,mj1:n]
    if len(down_p>0):
        for i in Segment(down_p,offset+mj+1,depth+1):
            yield i

def Main():
    global args,out

    args=ParseArg()
    if args.output=="stdout":
        out=sys.stdout
    else:
        out=open(args.output,"w")
    print >>out,"# INPUT FILE:",args.input
    print >>out,"# PVALUE CUT:",args.pvalue
    print >>out,"# MAX    GAP:",args.gap
    for chrom,pos,x2s,matrix_x2,snps in parseIterRegion(args.input):
        '''
        computer pvalue (from i to j)
        p[i][j] is the pvalue
        we are here to segments the snps into smaller region
        '''
        print >>sys.stderr,"BiGSEGMent",chrom,len(pos),pos[0],pos[-1]
        if args.fast:
            if len(pos)>1000:
                SqOR=0.0
                for k in range(len(pos)):
                    SqOR+=x2s[k]
                p=1.0-chi2.cdf(SqOR,len(pos))
                print >>out,"BIGREGION\t",chrom,"\t",pos[0],"\t",pos[-1]+1,"\t",p

                for i in range(len(pos)):
                    print >>out,"SNP\t",chrom,pos[i],snps[i],x2s[i],matrix_x2[i]
                print >>out,""
                continue




        p=numpy.array([[1.0 for row in range(len(pos))] for col in range(len(pos))])

        for i in range(len(pos)):
            SqOR=0.0
            print >>sys.stderr,i,"\r",
            for j in range(i,len(pos)):
                SqOR+=x2s[j]
                p[i][j]=1.0-chi2.cdf(SqOR,j+1-i)
 #               print "PV",i,j,p[i][j],SqOR,j+1-i
        for (start,end) in Segment(p):
            #print start,end,p[start][end]
            if (p[start][end]<args.pvalue):
                print >>out,"REGION\t",chrom,"\t",pos[start],"\t",pos[end]+1,"\t",p[start][end]

                for i in range(start,end+1):
                    print >>out,"SNP\t",chrom,pos[i],snps[i],x2s[i],matrix_x2[i]
                print >>out,""


            




    
if __name__=="__main__":
    Main()


