#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 17 Sep 2012 19:56:19

import os,sys,argparse
import pysam
import random
from math import log,sqrt
from xplib import TableIO
from xplib.Annotation import Bed
from scipy.stats import chi2


hNtToNum={'a':0,'A':0,
          'c':1,'C':1,
          'g':2,'G':2,
          't':3,'T':3
         }
Nt=['A','C','G','T']

def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency : pysam')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('-E','--elongation',dest="elongation",type=str,help="H3K36me3")
    p.add_argument('-P','--promoter',dest="promoter",type=str,help="H3K4me3 etc.")
    p.add_argument('-C','--control',dest="control",type=str,help="Input or DNA sequencing")
    p.add_argument('-o','--output',dest="output",type=str,default="stdout",help="output file")
    p.add_argument('--trans',dest="trans",type=str,help="trans unit file [output of xbam2tran.py]")
    p.add_argument('-p','--pvalue',dest="pvalue",type=float,default=0.001,help="allele pvalue cutoff")
    
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()


def iterSNPsInBedBAM(bed,bamA,bamB):
    disA=Bam2Dis(bed,bamA)
    disB=Bam2Dis(bed,bamB)

    for i in range(bed.length()):
        x=PosToChi2(disA,disB,i)
        if x:
            y=(bed.chr,bed.start+i,x[2],x[3],x[0],x[1])
            yield y
        
def PosToChi2(disA,disB,pos,min_coverage=10,min_snp=5):
    '''
    detect SNPs
    standard:
        1. coverage    >  10 (optional)
        2. two nucleid acid   
        3. LogRatio
    return :
        1. not a SNP
        2. logRatio
    '''
    a0=disA[pos]
    a1=disB[pos]
    s0=sum(a0)
    s1=sum(a1)
    if s0<min_coverage or s1<min_coverage :
        return None    #  information is not enough for decide.
    s=[0,0,0,0]
    idx=[0,0,0,0]
    for i in range(4):
        s[i]=a0[i]+a1[i]
    for i in range(4):
        for j in range(i+1,4):
            if s[i] < s[j]: idx[i]+=1
            if s[i] >= s[j]: idx[j]+=1
    for i in range(4):
        if idx[i]==0: idx1=i
        if idx[i]==1: idx2=i
    if s[idx2]< min_snp: return None  # no SNP
    (a11,a12,a21,a22)=(a0[idx1],a0[idx2],a1[idx1],a1[idx2])  # zero to one
    if (a11==0 or a12==0 or a21==0 or a22==0):
        ratio=(float(a11+0.5)/float(a12+0.5))/(float(a21+0.5)/float(a22+0.5))
        logratio=log(ratio)
        sigma2=1.0/(a11+1)+1.0/(a12+1)+1.0/(a21+1)+1.0/(a22+1)
    else:
        ratio=(float(a11)/float(a12))/(float(a21)/float(a22))
        logratio=log(ratio)
        sigma2=1.0/a11+1.0/a12+1.0/a21+1.0/a22
    x=logratio*logratio/sigma2
 
    return x,(a11,a12,a21,a22),idx1,idx2

def Bam2Dis(bed,bams):
    segment = [[0,0,0,0] for row in range(bed.length())]
    files=[]
    start=bed.start
    chrom=bed.chr
    end=bed.stop
    if type(bams) != type([1,2,3]):
        bams=[bams]
    for f in bams:
        try: 
            A=f.pileup(bed.chr,bed.start,bed.stop)
        except:
            print "can't pile up",bed
            continue
        for pileupcolumn in A:
            j0=pileupcolumn.pos-start
            if j0<0: continue
            if j0>=end: continue
            for pileupread in pileupcolumn.pileups:
                try:
                    nt=pileupread.alignment.seq[pileupread.qpos]
                    if hNtToNum.has_key(nt):
                        k0=hNtToNum[nt]
                        segment[j0][k0]+=1
                except:
                    pass
    return segment


def Main():
    global args,out
    args=ParseArg()
    if args.output=="stdout":
        out=sys.stdout
    else:
        try:
            out=open(args.output,"w")
        except IOError:
            print >>sys.stderr,"can't open file ",args.output,"to write. Using stdout instead"
            out=sys.stdout
    elongation=pysam.Samfile(args.elongation,'rb')
    control=pysam.Samfile(args.control,'rb')
    promoter=pysam.Samfile(args.promoter,'rb')
    for index,i in enumerate(TableIO.parse(args.trans,"transunit")):
        print >>sys.stderr,index," TransUnit\r",
        SqOR=0.0
        j=-1
        body_pvalue=1.0
        promoter_pvalue=1.0
        OUTSTR=""
        for j,x in enumerate(iterSNPsInBedBAM(i,elongation,control)):
            OUTSTR+="BODY SNP:\t"
            (chr,pos,alleleA,alleleB,x2,dist)=x
            OUTSTR+=chr+"\t"+str(pos)+"\t"+Nt[alleleA]+"/"+Nt[alleleB]+"\t"+str(x2)+"\t"+str(dist)+"\n"
            SqOR+=x2
        if j>=0:
            bp=1.0-chi2.cdf(SqOR,j+1)
            OUTSTR+="BODY PVALUE:\t"+str(bp)+"\n"
            body_pvalue=bp
        for k in i.promoters:
            SqOR=0.0
            j=-1
            PSTR=""
            PSTR+="PRMT:\t"+str(k)+"\n"
            for j,m in enumerate(iterSNPsInBedBAM(k,promoter,control)):
                (chr,pos,alleleA,alleleB,x2,dist)=m
                SqOR+=x2
                PSTR+="\tSNP:\t"
                PSTR+=chr+"\t"+str(pos)+"\t"+Nt[alleleA]+"/"+Nt[alleleB]+"\t"+str(x2)+"\t"+str(dist)+"\n"
            if j>=0:
                p=1.0-chi2.cdf(SqOR,j+1)
                PSTR+="PRMT PVALUE:\t"+str(p)+"\n"
                #print >>sys.stdout,PSTR
                if p<args.pvalue:
                    OUTSTR+=PSTR
                if promoter_pvalue>p:
                    promoter_pvalue=p
        if (promoter_pvalue<args.pvalue or body_pvalue<args.pvalue):
            print >>out,OUTSTR
            print >>out,i

    



    
if __name__=="__main__":
    Main()



