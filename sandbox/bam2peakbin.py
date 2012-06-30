#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 29 Jun 2012 16:23:45
import os,sys,argparse
import pysam
import prob


class BamBins:
    '''
        count the reads number , not the coverage.
    '''
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
            n_bin=l/binsize
            if l%binsize!=0: n_bin+=1
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
        self.binpeak=[]
        for i in range(len(self.bins)):
            print >>sys.stderr,"Call Pvalue of ",self.chrs[i],"\r",
            self.binpeak.append([False for row in range(len(self.bins[i]))])
            for j in range(len(self.bins[i])):
                if self.bins[i][j]==0: continue
                if prob.poisson_cdf(self.bins[i][j],self.lam,False)<pvalue:
                    self.binpeak[i][j]=True
    def PrintBinPeak(self,output="stdout"):
         if output=="stdout":
             out=sys.stdout
         else:
             try:
                out=open(output,"w")
             except:
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
        
        

        
        
        

    

def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -i file.bam -o file.bed', epilog='Library dependency : pysam')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('--input','-i',type=str,dest='Bamfile',help="the input alignment bamfile")
    p.add_argument('--output','-o',type=str,dest='Output',default="stdout", help="the output bed file  default: %(default)s")
    p.add_argument('--pvalue','-p',type=float,dest='Pvalue',default=1e-05, help="p-value threshold  default: %(default)f")
    p.add_argument('--bin','-b',type=int,dest='Binsize',default=1e-05, help="p-value threshold  default: %(default)f")
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)

    return p.parse_args()
def Main():
    global args,Table
    args=ParseArg()
    bambin=BamBins(args.Bamfile,args.Binsize)
    bambin.callBinPeak(args.Pvalue)
    bambin.PrintBinPeak(args.Output)
   

    
if __name__=="__main__":
    Main()


