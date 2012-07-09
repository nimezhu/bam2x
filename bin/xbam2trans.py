#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 09 Jul 2012 15:21:08
import os,sys,argparse
import pysam
import xplib.Stats.prob as prob
from xplib.Annotation.Utils import *
from ghmm import *
from xplib import TableIO
'''

'''
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
        self.chr2tid={}
        for i,chr in enumerate(self.chrs):
            self.chr2tid[chr]=i
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
            print >>sys.stderr,"Call Bins Pvalue in",self.chrs[i],"\r",
            self.binpeak.append([False for row in range(len(self.bins[i]))])
            for j in range(len(self.bins[i])):
                if self.bins[i][j]>=threshold:
               # if prob.poisson_cdf(self.bins[i][j],self.lam,False)<pvalue:
               # if prob.poisson_cdf(self.bins[i][j],self.lam,False)<pvalue:
                    self.binpeak[i][j]=True
        print >>sys.stderr,"Call Bins Pvalue Done                         "
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
    def iterPeaks(self):
        self.peaks=[]
        print >>sys.stderr,"Reading ",self.bamfilename," HMM Segments"
        k=0
        for i in self.parse_segment():
            k+=1
            if k%10000==0 :
                print >>sys.stderr,"parsed ",k," segments\r",
            (chr,start,stop)=i
        
            a=self.samfile.fetch(chr,start,stop)
            refine_start=None
            refine_stop=None
            reads_num=0
            for j in a:
                reads_num+=1
                if refine_start==None:
                    refine_start=j.pos
                if refine_stop==None:
                    refine_stop=j.pos+j.qlen
                if j.pos < refine_start: refine_start=j.pos
                if j.pos+j.qlen > refine_stop: refine_stop=j.pos+j.qlen
            peak_intensity=0
            peak_pos=0
            coverage=0.0
            for pileupcolumn in self.samfile.pileup(chr,refine_start,refine_stop):
                
                if pileupcolumn.n > peak_intensity:
                    peak_intensity=pileupcolumn.n
                    peak_pos=pileupcolumn.pos
                coverage+=pileupcolumn.n
            coverage/=(refine_stop-refine_start)
            lam=float(self.mapped)/(self.total_length/(refine_stop-refine_start))

            pvalue=prob.poisson_cdf(reads_num,lam,False)
            yield (chr,refine_start,refine_stop,reads_num,pvalue,coverage,peak_pos,peak_intensity)



        print >>sys.stderr,"Reading HMM Segments Done!                   " 

def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -E H3K36me3.bam -P H3K4me3.bam -o output.tab', epilog='Library dependency : pysam ghmm')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('--elongation','-E',type=str,dest='eBamfile',required=True,help="the input alignment bamfile of H3K36me3")
    p.add_argument('--promoter','-P',type=str,dest='pBamfile',required=True,help="the input alignment bamfile of H3K4me3 (or H3K27ac etc.)")
    p.add_argument('--output','-o',type=str,dest='Output',default="stdout", help="the output bed file  default: %(default)s")
    p.add_argument('--pvalue','-p',type=float,dest='Pvalue',default=1e-05, help="p-value threshold  default: %(default)f")
    p.add_argument('--bin','-b',type=int,dest='Binsize',default=200, help="binsize default: %(default)d")
    p.add_argument('--near','-n',type=int,dest='near',default=5000, help="Promoter less than  %(default)d are counted as Nearby Promoter")
    p.add_argument('--gene','-g',type=str,dest="genetab",help="Known Gene Tab (Download From UCSC genome browser )")
    p.add_argument('--filters','-f',dest='filters',nargs="+",action='store',help="filter other genome feature annotations(bed format). Example: Known_ncRNA.bed")
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()




    
    
    
def iterNearPromoter(chrom,start,stop):
    if not pTable.has_key(chrom):
        raise StopIteration
    near_start=start-args.near
    near_stop=stop+args.near
    if near_start<0:
        near_start=0
    for bin in iterRangeOverlapBins(near_start,near_stop):
        for p in pTable[chrom][bin]:
                if p[1] < near_stop and p[2] > near_start: 
                    yield p

def tab(a):
    s=str(a[0])
    for i in a[1:]:
        s+="\t"+str(i)
    return s
def Main():
    global args,Table,pTable
    near=5000
    args=ParseArg()
    if args.Output=="stdout":
        out=sys.stdout
    else:
        try:
          out=open(args.Output,"w")
        except IOerror:
          print >>sys.stderr,"can't write to %s",args.Output
          out=sys.stdout
    
    if args.genetab:
        genes=readIntoBinIndex(TableIO.parse(args.genetab,'genebed'))
    other_filters=[]
    if args.filters:
        for i in args.filters:
            other_filters.append(readIntoBinIndex(TableIO.parse(i,'bed')))
    
    ebambin=BamBins(args.eBamfile,args.Binsize)
    ebambin.process(args.Pvalue)
    pbambin=BamBins(args.pBamfile,args.Binsize)
    pbambin.process(args.Pvalue)
    pTable=readIntoBinIndex(pbambin.iterPeaks())
    print >>out,"# this file was generated by %(prog)s"
    print >>out,"# the columns are  id,(chrom,start,stop,reads_number,pvalue,coverage,peak_pos,peak_intesity),guess_strand"
    print >>out,"# NV is Novel, NP is No Promoter, KN is Known"
    print >>out,"# Guess Strand Choices : None ( No promoter nearby)"
    print >>out,"#                        .    ( both ends have promoters)"
    print >>out,"#                        +    ( only 5' end has promoter)"
    print >>out,"#                        -    ( only 3' end has promoter)"
    print >>out,"# Options:"
    print >>out,"# Pvalue cutoff:",args.Pvalue
    print >>out,"# Neayby Promoter:",args.near,"bp"
    print >>out,"# H3K36me3 bam file"
    print >>out,ebambin
    print >>out,"# Promoter bam file"
    print >>out,pbambin
    k=0
    for i in ebambin.iterPeaks():
        (chrom,start,stop,reads_num,pvalue,coverage,peak_pos,peak_i)=i
        k+=1
        t=0
        plus=0
        plus_reads=0
        minus=0
        minus_reads=0  
        thres=500 # core promoter down 500bp? 
        if (stop-start)<1000:
            thres=(stop-start)/2
        for p in iterNearPromoter(chrom,start,stop):
                if p[1]<start+thres and p[2] < stop-thres: 
                    plus+=1
                    plus_reads+=p[3]
                if p[2]>stop-thres and p[1] > start+thres: 
                    minus+=1
                    minus_reads+=p[3]
                t+=1
        gene=[]
        if args.genetab:
            gene=getOverlapFeatures(Bed([chrom,start,stop]),genes)
        fs=[]
        if len(other_filters)>0:
            for f in other_filters:
                for o in getOverlapFeatures(Bed([chrom,start,stop]),f):
                    fs.append(o)
        if len(gene)==0 and len(fs)==0 and t>0:
            KorN="NV"
        elif len(gene) > 0 or len(fs) > 0:
            KorN="KN"
        elif t==0:
            KorN="NP"
        print >>out,KorN,"\t",
        print >>out,"Transcript_"+str(k),"\t",tab(i),
        if t==0:
            print >>out,"\t","None"
        elif minus>0 and plus==0:
            print >>out,"\t","-"
        elif plus>0 and minus==0:
            print >>out,"\t","+"
        else:
            print >>out,"\t","."
        print >>out,KorN,"\tPromoter Info:\t[ Total: "+str(t)+",",str(plus)+" >>>> "+str(minus)+" ]"
        for g in gene:
            print >>out,KorN,"\tOverlapGene:\t",g
        if len(other_filters)>0:
            for f in fs:
                print >>out,KorN,"\tOverlapFeat:\t",f
        for p in iterNearPromoter(chrom,start,stop):
            print >>out,KorN,"\tNearbyPromoter:\t",p,
            if p[1] < start+thres and p[2] < stop-thres:
                print >>out,"+"
            elif p[1] > start+thres and p[2] > stop-thres:
                print >>out,"-"
            else:
                print >>out,"."
            
        
       
        print >>out,"//\n"


    
if __name__=="__main__":
    Main()


