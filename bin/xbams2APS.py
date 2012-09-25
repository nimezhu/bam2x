#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 25 Sep 2012 15:15:50

import os,sys,argparse
import pysam
import random
from math import log,sqrt

def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency : pysam')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('--bamA',nargs='*',dest="bamA",action="store",default=[],help="add bam to list A")
    p.add_argument('--bamB',nargs='*',dest="bamB",action="store",default=[],help="add bam to list B")
    p.add_argument('--bamlistA',dest="bamlistA",action="store",default=None,type=str,help="add filename in file to list A")
    p.add_argument('--bamlistB',dest="bamlistB",action="store",default=None,type=str,help="add filename in file to list B")
    p.add_argument('--bin',dest="binsize",action="store",default=1000000,help="count binsize bp each time")
    p.add_argument('-o','--out',dest="out",type=str,action="store",default="stdout",help="output file")
    p.add_argument('--min_coverage',dest="min_coverage",type=int,action="store",default=10,help="not count the chi-square if either bam coverage is less than this value [%(default)i]")
    p.add_argument('--min_minor_snp_cov',dest="min_snp",type=int,action="store",default=5,help="not count the chi-square if minor snp coverage (add coverage in two bam) is less than this value [%(default)i]")
    p.add_argument('-r','--region',dest="region",type=str,action="store",default=None,help="only report region example: chr1:1-1000")
   
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()

hNtToNum={'a':0,'A':0,
          'c':1,'C':1,
          'g':2,'G':2,
          't':3,'T':3
         }
Nt=['A','C','G','T']
def parse(region):
    a=region.split(":")
    chrom=a[0]
    (start,stop)=a[1].split("-")
    start=int(start)-1
    stop=int(stop)
    return (chrom,start,stop)
    
def processing():
    a=0
    for j in range(binsize):
        a=diff(j,0,1)


def diff(pos,i0,i1):
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
    a0=segments[i0][pos]
    a1=segments[i1][pos]
    s0=sum(a0)
    s1=sum(a1)
    if s0<args.min_coverage or s1<args.min_coverage :
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
    if s[idx2]<args.min_snp: return None  # no SNP
    (a11,a12,a21,a22)=(a0[idx1],a0[idx2],a1[idx1],a1[idx2])  # zero to one
    print >>out,chrom,"\t",offset+pos,"\t",Nt[idx1]+"/"+Nt[idx2],"\t",
    if (a11==0 or a12==0 or a21==0 or a22==0):
        ratio=(float(a11+0.5)/float(a12+0.5))/(float(a21+0.5)/float(a22+0.5))
        logratio=log(ratio)
        sigma2=1.0/(a11+1)+1.0/(a12+1)+1.0/(a21+1)+1.0/(a22+1)
    else:
        ratio=(float(a11)/float(a12))/(float(a21)/float(a22))
        logratio=log(ratio)
        sigma2=1.0/a11+1.0/a12+1.0/a21+1.0/a22
    x=logratio*logratio/sigma2
    print >>out,x,"\t",
    print >>out,"(",a11,a12,a21,a22,")\t",a0,"vs",a1
    return x


    
def Main():
    global args,chrom,offset,segments,binsize,out
    args=ParseArg()
    if args.out=="stdout":
        out=sys.stdout
    else:
        try:
            out=open(args.out,"w")
        except IOError:
            print >>sys.stderr,"can't open file",args.out,"to write, using stdout instead"
            out=sys.stdout
    binsize=args.binsize
    segments=[]
    if args.region!=None:
	(chrom,start,end)=parse(args.region)
        binsize=end-start
    for i in range(2):
        segments.append([[0,0,0,0] for row in range(binsize)])
    Afiles=[]
    Bfiles=[]
    for bamB in args.bamB:
        Bfiles.append(pysam.Samfile(bamB,"rb"))
    for bamA in args.bamA:
        Afiles.append(pysam.Samfile(bamA,"rb"))
    if args.bamlistA:
        f=open(args.bamlistA,"r")
        for line in f:
            line =line.strip()
            Afiles.append(pysam.Samfile(line,"rb"))
        f.close()
    if args.bamlistB:
        f=open(args.bamlistB,"r")
        for line in f:
            line =line.strip()
            Bfiles.append(pysam.Samfile(line,"rb"))
        f.close()
    #=====================
    chrs=()
    lengths=()
    for f in Afiles:
        new_chrs=f.references
        new_lengths=f.lengths
        if len(chrs)<len(new_chrs):  # some bam don't have chrY
            chrs=new_chrs
            lengths=new_lengths
    for f in Bfiles:
        new_chrs=f.references
        new_lengths=f.lengths
        if len(chrs)<len(new_chrs):  # some bam don't have chrY
            chrs=new_chrs
            lengths=new_lengths



    #=====================
    hChrToTid={}
    hChrToLength={}
    for i,x in enumerate(chrs):
        hChrToTid[x]=i
        hChrToLength[x]=lengths[i]

    print >>out,"# Compare The SNP preferences between files:"
    print >>out,"# Group A:"
    for f in args.bamA:
        print >>out,"#",f
    if args.bamlistA:
        f=open(args.bamlistA)
        for line in f:
            line=line.strip()
            print >>out,"#",line
        f.close()
    print >>out,"# Group B:"
    for f in args.bamB:
        print >>out,"#",f
    if args.bamlistB:
        f=open(args.bamlistB)
        for line in f:
            line=line.strip()
            print >>out,"#",line
        f.close()
    print >>out,"# Parameters:"
    print >>out,"# Minimal Coverage:          ",args.min_coverage
    print >>out,"# Minimal Minor SNP Coverage:",args.min_snp
    if args.region==None:
        for i,chrom in enumerate(chrs):
           for j in range(lengths[i]/binsize+1):
               offset=j*binsize
               start=offset
               end=offset+binsize
               for i0 in range(len(segments)):
                   for j0 in range(binsize):
                       for k0 in range(4):
                           segments[i0][j0][k0]=0
               k=0
               for f in Afiles:
                   try: 
                       A=f.pileup(chrom,start,end)
                   except:
                       print >>sys.stderr,"can't pile up",chrom,start,end
                       continue
                   for pileupcolumn in A:
                       i0=k
                       j0=pileupcolumn.pos-offset
                       if j0<0: continue
                       if j0>=binsize: break
                       for pileupread in pileupcolumn.pileups:
                           try:
                               if pileupread.is_del : continue
                               if pileupread.indel : continue
                               nt=pileupread.alignment.seq[pileupread.qpos]
                               if hNtToNum.has_key(nt):
                                   k0=hNtToNum[nt]
                                   segments[i0][j0][k0]+=1
                           except:
                               pass
               k=1 
               for f in Bfiles:
                   try: 
                       B=f.pileup(chrom,start,end)
                   except:
                       print >>sys.stderr,"can't pile up",chrom,start,end
                       continue
                   for pileupcolumn in B:
                       i0=k
                       j0=pileupcolumn.pos-offset
                       if j0<0: continue
                       if j0>=binsize: break
                       for pileupread in pileupcolumn.pileups:
                           try:
                               if pileupread.is_del : continue
                               if pileupread.indel : continue
                               nt=pileupread.alignment.seq[pileupread.qpos]
                               if hNtToNum.has_key(nt):
                                   k0=hNtToNum[nt]
                                   segments[i0][j0][k0]+=1
                           except:
                               pass
               a=processing()
               if start%100000==0:
                   print >>sys.stderr,chrom,"\t",start,"bp now                   \r",
        
    else:
                   k=0
                   offset=start
	           for f in Afiles:
	               try: 
	                   A=f.pileup(chrom,start,end)
	               except:
	                   continue
	               for pileupcolumn in A:
	                   i0=k
	                   j0=pileupcolumn.pos-offset
	                   if j0<0: continue
	                   if j0>=binsize: break
	                   for pileupread in pileupcolumn.pileups:
	                       try:
                                   if pileupread.is_del : continue
                                   if pileupread.indel : continue
	                           nt=pileupread.alignment.seq[pileupread.qpos]
	                           if hNtToNum.has_key(nt):
	                               k0=hNtToNum[nt]
	                               segments[i0][j0][k0]+=1
	                       except:
	                           pass
	           k=1
	           for f in Bfiles:
	               try: 
	                   B=f.pileup(chrom,start,end)
	               except:
	                   continue
	               for pileupcolumn in B:
	                   i0=k
	                   j0=pileupcolumn.pos-offset
	                   if j0<0: continue
	                   if j0>=binsize: continue
	                   for pileupread in pileupcolumn.pileups:
	                       try:
                                   if pileupread.is_del : continue
                                   if pileupread.indel : continue
	                           nt=pileupread.alignment.seq[pileupread.qpos]
	                           if hNtToNum.has_key(nt):
	                               k0=hNtToNum[nt]
	                               segments[i0][j0][k0]+=1
	                       except:
	                           pass
	           a=processing()
    
if __name__=="__main__":
    Main()


