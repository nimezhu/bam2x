#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 07-09-2012, 13:56:34 CDT

import os,sys,argparse
import pysam

def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency : pysam')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('--bam',nargs='*',dest="bams",action="store",default=[],help="add bam to list")
    p.add_argument('--bin',dest="binsize",action="store",default=1000,help="binsize")
   
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()

hNtToNum={'a':0,'A':0,
          'c':1,'C':1,
          'g':2,'G':2,
          't':3,'T':3
         }
def processing():
    a=0
    for j in range(binsize):
        for i0 in range(len(args.bams)):          
            for i1 in range(i0,len(args.bams)):
                if diff(j,i0,i1):
                    a+=1
    return a    
def diff(pos,i0,i1):
    a0=segments[i0][pos]
    a1=segments[i1][pos]
    s0=sum(a0)
    s1=sum(a1)
    if s0==0 or s1==0:
        return 0
    max0=0
    max1=0
    for i in range(4):
        if a0[max0]< a0[i]: max0=i
        if a1[max1]< a1[i]: max1=i
    if max0==max1: return 0

    #### need revised to t test ######### 
    #### for simplicity
    return 1



    
def Main():
    global args,chrom,offset,segments,binsize
    
    args=ParseArg()
    binsize=args.binsize
    segments=[]
    for i in range(len(args.bams)):
        segments.append([[0,0,0,0] for row in range(binsize)])
    files=[]
    for bam in args.bams:
        print bam
        files.append(pysam.Samfile(bam,"rb"))
    #=====================
    chrs=()
    lengths=()
    for f in files:
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
    for i,chrom in enumerate(chrs):
        for j in range(lengths[i]/binsize):
            offset=j*binsize
            start=offset
            end=offset+binsize
            for i0 in range(len(segments)):
                for j0 in range(binsize):
                    for k0 in range(4):
                        segments[i0][j0][k0]=0
            for k,f in enumerate(files):
                for pileupcolumn in f.pileup(chrom,start,end):
                    i0=k
                    j0=pileupcolumn.pos-offset
                    if j0<0: continue
                    if j0>=binsize: continue
                    for pileupread in pileupcolumn.pileups:
                        nt=pileupread.alignment.seq[pileupread.qpos]
                        if hNtToNum.has_key(nt):
                            k0=hNtToNum[nt]
                            segments[i0][j0][k0]+=1
            a=processing()
            if a>0: print chrom,"\t",start,"\t",end,"\t",a
            if start%1000000==0:
                print >>sys.stderr,chrom,"\t",start,"bp now                   \r",
    
    
if __name__=="__main__":
    Main()


