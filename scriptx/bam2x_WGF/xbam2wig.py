#!/usr/bin/python
# Wei Guifeng, <guifengwei@gmail.com>
import os,re,sys,argparse
import pysam
from collections import defaultdict

'''
This script was used to make the wiggle file from ChIP-seq Bam
And you can generate the bigwig with the Kent utility: wigToBigWig

'''

## Compute the read intensity
## binsize bp resolution 
## signal intensity = 1,000,000.0 * reads_counts/ total_mapped

class BamBins:
    '''
        count the reads number , not the coverage.
    '''
    def __init__(self,bamfilename,binsize=200,ReadFile=True):
        self.bins=defaultdict(list)
        self.bamfilename=bamfilename
        self.binsize=binsize
        self.samfile=pysam.Samfile(bamfilename,"rb")
        self.chrs=self.samfile.references
        self.lengths=self.samfile.lengths
        self.chromsize={}
        self.chr2tid={}
        for i,chrs in enumerate(self.chrs): 
            self.chromsize[chrs]=self.lengths[i]
            self.chr2tid[chrs]=i
        '''
        initialize binss
        '''
        for chrs,length in self.chromsize.iteritems():
            n_bin=length/self.binsize
            if length%self.binsize!=0: n_bin+=1
            self.bins[self.chr2tid[chrs]]=[0 for row in range(n_bin)]
        if ReadFile:
            self.mapped=0
            self.unmapped=0
            self.readBam()
    def __str__(self):
        s="# Bamfile: "+ os.path.basename(self.bamfilename)+"\n"
        s+="# Binsize: "+ str(self.binsize)+"\n"
        s+="# Mapped Reads Number: "+ str(self.mapped)+"\n"
        s+="# Unmapped Reads Number: "+ str(self.unmapped)+"\n"
        return s
    def readBam(self):
        '''
        read all alignment reads in bamfile into bins
        '''
        print >>sys.stderr,"# reading ",self.bamfilename
        i=0
        for s in self.samfile:
            i+=1
            if i%100000==0: print >>sys.stderr,i,"reads\r",
            if s.tid==-1:
                self.unmapped+=1
                continue
            self.mapped+=1
            start,stop=0,0
            if s.is_reverse:
                start,stop=s.pos+s.qlen-extend, s.pos+s.qlen
                if start<=0:
                    start=1
            else:
                start,stop=s.pos,s.pos+extend
                if stop>self.lengths[s.tid]:
                    stop= self.lengths[s.tid]
            startBin,stopBin=start/self.binsize,stop/self.binsize
            if (start-startBin*self.binsize)>self.binsize/2:
                startBin=startBin+1
            if (stop-stopBin*self.binsize)<self.binsize/2:
                stopBin=stopBin-1
            for k in range(startBin,stopBin+1):
                self.bins[s.tid][k]+=1
        print >>sys.stderr,"\r# reading ",os.path.basename(self.bamfilename)," done\t\t\t"
        if norm:
            for chrs in self.chrs:
                print >>sys.stderr,"# Normalizing ... ... \r",
                for bin_id,value in enumerate(self.bins[self.chr2tid[chrs]]):
                    # read intensity every i million mapped reads
                    self.bins[self.chr2tid[chrs]][bin_id] *= 1000000.0/float(self.mapped)
            print >>sys.stderr,"# Normalizing ... ... done"
    def correct(self,other):
        '''
        correct the ChIP-seq Signal from the control file
        '''
        for chrs in self.chrs:
            for bin_id,signal in enumerate(self.bins[self.chr2tid[chrs]]):
                self.bins[self.chr2tid[chrs]][bin_id] -= other.bins[other.chr2tid[chrs]][bin_id]
                if self.bins[self.chr2tid[chrs]][bin_id]<0:
                    self.bins[self.chr2tid[chrs]][bin_id]=0
        print >>sys.stderr,"# Correct ... ... OK"

def ArrayToWig(a):
    ''' 
    bambins to wiggles file,resolution=Binsize
    '''
    state=0
    start=0
    stop=0
    segment=[]
    for i,x in enumerate(a):
        if x>0 and state==0:
            start=i
            state=1
            stop=i
        if x>0 and state==1:
            stop=i
        if x==0 and state==0:
            continue
        if x==0 and state==1:
            state=0
            segment.append((start,stop))
    if state==1:
        segment.append((start,stop))
    return segment

def parse_argument():
    '''
        argument parser
    '''
    p=argparse.ArgumentParser(description='example: %(prog)s -b binsize -e extend -C input.bam --bam chip.bam -w output.wig',epilog='dependency pysam')
    # Binsize
    p.add_argument("-b","--binsize",dest="binsize",type=int,metavar="Binsize",default=100,
            help="The genome will be divided into binsize bp bins. default=%(default)d")
    # extend the reads in the 3'-direction.
    # the estimated median size of ChIP fragment is 200bp
    p.add_argument('-e','--extend',dest='extend',type=int,metavar="Extend",default=200,help="extend the reads in the 3'-direction.default=%(default)d")
    # the Control experimnt
    p.add_argument("-c",'--control',dest="control",type=str,metavar="Control",help="The Input bam")
    # the chip experiment
    p.add_argument("--bam",dest="bam",type=str,metavar='Bam',help="The ChIP experiment bam file")
    # show help if no argument
    p.add_argument("-w",'--wiggle',dest='wiggle',metavar="Wiggle",required=True,type=str,help="The wiggles output")
    if len(sys.argv)==1:
        sys.exit(p.print_help())
    args = p.parse_args()
    return args

def main():
    global args,norm,extend
    # obtain the args
    args=parse_argument()
    control=args.control
    bam=args.bam
    binsize=args.binsize
    norm = True
    extend=args.extend

    out=open(args.wiggle,'w')
    print >>sys.stderr,"# xbam2wig.py 0.1.4 "
    print >>sys.stderr,"# ChIP: ",os.path.basename(bam)
    if args.control:
        print >>sys.stderr,"# Control: ",os.path.basename(control)
    print >>sys.stderr,"# Wiggle: ",args.wiggle
    print >>sys.stderr,"# Binsize: ",binsize
    print >>sys.stderr,"# Extend: ",extend

    # bam to bins
    chip_bin=BamBins(bam,binsize=binsize)
    if args.control:
        control_bin=BamBins(control,binsize=binsize)
    # normlized input signal
        chip_bin.correct(control_bin)
    print >>out,"track type=wiggle_0 name="+os.path.basename(bam)+" description="+os.path.basename(bam)
    for i,chrs in enumerate(chip_bin.bins):
        print >>sys.stderr,"# Making the wiggles ... ...", chip_bin.chrs[i],"\r",
        s=ArrayToWig(chip_bin.bins[i])
        for each in s:
            (start,stop)=each
            print >>out, "fixedStep chrom="+chip_bin.chrs[i],"start="+str(start*binsize+1),"step="+str(binsize),"span="+str(binsize)
            for k in range(start,stop+1):
                print >>out, "%.3f" %(chip_bin.bins[i][k])
    out.close()
    print >>sys.stderr, "# Making the wiggles ... ... done!"


if __name__=="__main__":
    main()

