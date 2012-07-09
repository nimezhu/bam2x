#!/usr/bin/python
# Wei Guifeng, <guifengwei@gmail.com>
'''
    Calculate the reproducibility, like as:
    
    H3K4me3rep1.bam    H3K4me3rep2.bam 

    Calculate the coefficient between similar chromatin state marker, like as:

    H3K4me1.bam        H3K27ac.bam
'''

import os,sys,argparse
import pysam
import scipy.stats



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
        print >>sys.stderr,"\rreading ",os.path.basename(self.bamfilename)," done\t\t\t"
        if norm:
            for chrs in self.chrs:
                print >>sys.stderr,"Normalizing ... ... \r",
                for each_bin in self.bins[self.chr2tid[chrs]]:
                    each_bin = 100000000.0*each_bin/float(self.binsize*self.mapped)
            print >>sys.stderr,"Normalizing ... ... done"
    def makevector(self):
        ''' 
        combine all bins into one vector
        '''
        vector=[]
        for chrs in self.chrs:
            vector+= self.bins[self.chr2tid[chrs]]
        return vector


def pearson(x=None,y=None):
    '''
        Calculate the pearson correlation coefficient
    '''
    if isinstance(x,list) and isinstance(y,list):
        if len(x) == len(y):
            coeff,pvaule=scipy.stats.pearsonr(x,y)
            return coeff
        else:
            print >>sys.stderr,"The length between x and y is not equal."
    else:
        print >>sys.stderr,"The type of x and y must be list!"
def parse_argument():
    '''
        argument parser
    '''
    p=argparse.ArgumentParser(description='Example: %(prog)s -b Binsize -N -P H3K4me1.bam -Q H3K27ac.bam',
        epilog="This script is used to calculate the pearson correlation of coefficient between bam")
    p.add_argument('-V',"--version",action='version',version='%(prog)s 0.1.1')
    p.add_argument('--bin','-b',type=int,dest='Binsize',default=200, help="binsize default: %(default)d")
    p.add_argument('-N','--norm',action='store_true',dest="norm",
                    help="Normalization. RPKM = 10^9*C/(N*L)")
    p.add_argument('-P','-p',type=str,dest="pbam",help="The 1st bam file")
    p.add_argument('-Q','-q',type=str,dest="qbam",help="The 2nd bam file")

    if len(sys.argv)==1:
        sys.exit(p.print_help())
    args = p.parse_args()
    return args


def test():
    global norm

    # process the bam file
    args = parse_argument()
    norm = args.norm
    print >>sys.stderr,"# xbam_correlation.py"
    print >>sys.stderr,"# pBamfile: ",os.path.basename(args.pbam)
    print >>sys.stderr,"# qBamfile: ",os.path.basename(args.qbam)
    print >>sys.stderr,"# Binsize: ", args.Binsize
    print >>sys.stderr,"# Normalization: ", norm

    pbambins= BamBins(args.pbam,binsize=args.Binsize)
    qbambins= BamBins(args.qbam,binsize=args.Binsize)

	#if (pbambins.mapped<10000000) or (qbambins.mapped<10000000) or abs( (pbambins-qbambins)>5000000 ):
	#	if not norm:
	#		print >>sys.stderr,"# warning: recommend -N "

    p_vector = pbambins.makevector()
    q_vector = qbambins.makevector()

    coeff = pearson(x=p_vector,y=q_vector)
    print >>sys.stderr, "# The Pearson Correlation of Coefficient is",coeff

        
if __name__ == '__main__':
    test()
