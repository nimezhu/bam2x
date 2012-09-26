#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 26 Sep 2012 09:54:13

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
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -R SNP.bed -A A_Rep1.bam A_Rep2.bam -B B_Rep1.bam B_Rep2.bam -o output.tab --output2 SNP.NoInfo.tab', epilog='Library dependency : pysam')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('-C','--control',dest="control",type=str,help="Input or DNA sequencing")
    p.add_argument('-o','--output',dest="output",type=str,default="stdout",help="output file")
    p.add_argument('--output2',dest="output2",type=str,default="stderr",help="output position that have not enough information")
    p.add_argument('-A','--bamA',dest="bamA",action='append',default=[],required=True,help="bams from Sample A")
    p.add_argument('-B','--bamB',dest="bamB",action='append',default=[],required=True,help="bams from Sample B")
    p.add_argument('-R','--region',dest="region",action='store',type=str,required=True,help="SNP bed file")
    
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()


def iterPos(bed,bamA,bamB):
    disA=Bam2Dis(bed,bamA)
    disB=Bam2Dis(bed,bamB)
    for i in range(bed.length()):
        x=PosToChi2(disA,disB,i)
        y=(bed.chr,bed.start+i,x[2],x[3],x[0],x[1])
        yield y
        
def PosToChi2(disA,disB,pos):
    '''
    return :
        logRatio
    '''
    a0=disA[pos]
    a1=disB[pos]
    s0=sum(a0)
    s1=sum(a1)
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
    (a11,a12,a21,a22)=(a0[idx1],a0[idx2],a1[idx1],a1[idx2])  # zero to one
    if (a11==0 or a12==0 or a21==0 or a22==0):
        ratio=(float(a11+0.5)/float(a12+0.5))/(float(a21+0.5)/float(a22+0.5))
        logratio=log(ratio)
        sigma2=1.0/(a11+1)+1.0/(a12+1)+1.0/(a21+1)+1.0/(a22+1)
    else:
        ratio=(float(a11)/float(a12))/(float(a21)/float(a22))
        logratio=log(ratio)
        sigma2=1.0/a11+1.0/a12+1.0/a21+1.0/a22
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
                    if pileupread.is_del:continue
                    if pileupread.indel!=0: continue
                    nt=pileupread.alignment.seq[pileupread.qpos]
                    if hNtToNum.has_key(nt):
                        k0=hNtToNum[nt]
                        segment[j0][k0]+=1
                except:
                    pass
    return segment
def tab(a):
    S=str(a[0])
    for i in a[1:]:
        S+="\t"+str(i)
    return S
    

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
    if args.output2=="stderr":
        out2=sys.stderr
    else:
        try:
            out2=open(args.output2,"w")
        except IOError:
            print >>sys.stderr,"can't open file",args.output2
            out2=sys.stderr
    bamA=[]
    bamB=[]
    print >>out,"# SNPs Annotation:"
    print >>out,"# ",args.region
    print >>out,"# BAM A:"
    for f in args.bamA:
        bamA.append(pysam.Samfile(f,"rb"))
        print >>out,"#\t",f
    print >>out,"# BAM B:"
    for f in args.bamB:
        bamB.append(pysam.Samfile(f,"rb"))
        print >>out,"#\t",f
    pos1_num=0
    pos2_num=0
    for index,i in enumerate(TableIO.parse(args.region,"simple")):
        print >>sys.stderr,index," bed entry\r",
        b=Bed(i)
        for j,x in enumerate(iterPos(b,bamA,bamB)):
            a=x[5]

            if sum(a)==0:
                pos2_num+=1
                print >>out2,tab(i),"\tNoInfo"
            elif a[1]+a[3]==0:
                pos2_num+=1
                print >>out2,tab(i),"\t"+Nt[x[2]]+" Only"
            else:
                pos1_num+=1
                S=""
                S+=x[0]+"\t"+str(x[1])+"\t"+Nt[x[2]]+"/"+Nt[x[3]]+"\t"+str(x[4])+"\t"+str(x[5])
                print >>out,tab(i),"\t",S
    num=pos1_num+pos2_num
    print >>sys.stderr,pos1_num,"/",num,"SNPs have reads coverage on two type of nucleid acid"
    print >>sys.stderr,pos2_num,"/",num,"SNPs don't have enough information to guess allele preference"

    
if __name__=="__main__":
    Main()



