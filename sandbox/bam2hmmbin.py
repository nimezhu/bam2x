#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 01 Jul 2012 00:19:07
import os,sys,argparse
import pysam
import prob
from ghmm import *

class BamBins:
    '''
        count the reads number , not the coverage.
    '''
    sigma=IntegerRange(0,2)
    def __init__(self,bamfilename,binsize=200,ReadFile=True):
        self.bins=[]
        self.bamfilename=bamfilename
        self.binsize=binsize
        self.samfile=pysam.Samfile(bamfilename,"rb")
        self.chrs=self.samfile.references
        self.lengths=self.samfile.lengths
        '''
        initialize binss
        '''
        self.total_length=0
        for l in self.lengths:
            self.total_length+=l
            n_bin=l/self.binsize
            if l%self.binsize!=0: n_bin+=1
            self.bins.append([0 for row in range(n_bin)])
        if ReadFile:
            self.mapped=0
            self.unmapped=0
            self.readBam()
            self.lam=float(self.mapped)/(self.total_length/self.binsize)
        
    def __str__(self):
        s="# Bamfile: "+self.bamfilename+"\n"
        s+="# Binsize: "+str(self.binsize)+"\n"
        s+="# Mapped Reads Number: "+ str(self.mapped)+"\n"
        s+="# Unmapped Reads Number: "+ str(self.unmapped)+"\n"
        s+="# Lambda: "+ "%.4f"%self.lam +"\n"
        return s
        
    def readBam(self):
        '''
        read all alignment reads in bamfile  into bins
        '''
        print >>sys.stderr,"reading ",self.bamfilename
        i=0
        for s in self.samfile:
            i+=1
            if i%100000==0: print >>sys.stderr,i,"reads\r",
            if s.tid==-1:
                self.unmapped+=1
                continue
            self.mapped+=1
            bin_id=(s.pos+s.qlen/2)/self.binsize
            self.bins[s.tid][bin_id]+=1
        print >>sys.stderr,"\rreading ",self.bamfilename," done\t\t\t"
    
    def callBinPeak(self,pvalue=1e-05):
        self.pvalue=pvalue
        self.binpeak=[]
        threshold=1
        while 1:
            if prob.poisson_cdf(threshold,self.lam,False) < pvalue: break
            threshold+=1
        for i in range(len(self.bins)):
            print >>sys.stderr,"Call Pvalue of ",self.chrs[i],"\r",
            self.binpeak.append([False for row in range(len(self.bins[i]))])
            for j in range(len(self.bins[i])):
                if self.bins[i][j]>=threshold:
               # if prob.poisson_cdf(self.bins[i][j],self.lam,False)<pvalue:
               # if prob.poisson_cdf(self.bins[i][j],self.lam,False)<pvalue:
                    self.binpeak[i][j]=True
    def PrintBinPeak(self,output="stdout"):
         if output=="stdout":
             out=sys.stdout
         else:
             try:
                out=open(output,"w")
             except IOerror:
                print >>sys.stderr,"Warning: can't write to file",output,", using stdout instead."
                out=sys.stdout
         ind=0
         print >>out,str(self)
         for i in range(len(self.bins)):
            for j in range(len(self.bins[i])):
                if self.binpeak[i][j]:
                    ind+=1
                    s=self.chrs[i]+"\t"+str(j*self.binsize)+"\t"+str((j+1)*self.binsize)
                    s+="\t"+"BinPeak_"+str(ind)
                    s+="\t"+str(self.bins[i][j])
                    s+="\t."
                print >>out,s
    def trainHMM(self,seq=None):
        T=[[0.9,0.1],[0.1,0.9]]
        e1=[0.1,0.9]
        e0=[0.9,0.1]
        E=[e0,e1]
        pi=[0.9,0.1] # initial 10% are peak?
        sigma=BamBins.sigma
       # if not seq:
       #     seq=self.binpeak[0]
       # print seq
        m = HMMFromMatrices(sigma,DiscreteDistribution(sigma),T,E,pi)
        m.baumWelch(EmissionSequence(sigma,self.binpeak[0]))
        self.m=m
    def decodeHMM(self):
        self.hmm_states=[]
        for i in self.binpeak:
            a,score=self.m.viterbi(EmissionSequence(BamBins.sigma,i))
            self.hmm_states.append(a)
    def write(self,out="stdout"):
        for i,chr in enumerate(self.chrs):
            if out=="stdout": out=sys.stdout
            if type(out)==type("s"): 
                try:
                    out=open(out,"w")
                except IOerror:
                    print >>sys.stderr,"can't open file [",out,"] to write"
            print >>out,"chrom=",chr
            for j in range(len(self.bins[i])):
                print >>out,j*self.binsize,self.bins[i][j],self.binpeak[i][j],self.hmm_states[i][j]

    def process(self,pvalue=1e-05):
        self.callBinPeak(pvalue)
        self.trainHMM()
        self.decodeHMM()
    def parse_segment(self):
        '''
         Usage :  for i in bambin.parse_segment()
        '''
        for i,x in enumerate(self.hmm_states):
            state=0
            start=0
            stop=0
            for j,y in enumerate(x):
                if state==0:
                    if y==1:
                        state=1
                        start=j
                        stop=j
                    if y==0:
                        continue
                if state==1:
                    if y==1:
                        stop=j
                    if y==0:
                        yield (self.chrs[i],start*self.binsize,(stop+1)*self.binsize)
                        state=0
            if state==1:
                 yield (self.chrs[i],start*self.binsize,(stop+1)*self.binsize)






def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -i file.bam -o file.bed', epilog='Library dependency : pysam')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('--input','-i',type=str,dest='Bamfile',help="the input alignment bamfile")
    p.add_argument('--output','-o',type=str,dest='Output',default="stdout", help="the output bed file  default: %(default)s")
    p.add_argument('--pvalue','-p',type=float,dest='Pvalue',default=1e-05, help="p-value threshold  default: %(default)f")
    p.add_argument('--bin','-b',type=int,dest='Binsize',default=200, help="p-value threshold  default: %(default)f")
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)

    return p.parse_args()
def Main():
    global args,Table
    args=ParseArg()
    bambin=BamBins(args.Bamfile,args.Binsize)
    bambin.process(args.Pvalue)
    print bambin.m
    for i in bambin.parse_segment():
        print i
  
        
   

    
if __name__=="__main__":
    Main()


