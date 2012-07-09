#!/usr/bin/python
# Wei Guifeng, <guifengwei@gmail.com>
'''
    This script is used to calculate the pearson correlation of coefficient between bam

    OUTPUT: cor_matrix.txt
'''
import os,re,sys,argparse
import pysam
import scipy.stats
import numpy as np
from collections import defaultdict

#### hg19_coordinats
_references=('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM')

_lengths=(249250621L, 243199373L, 198022430L, 191154276L, 180915260L, 171115067L, 159138663L, 146364022L, 141213431L, 135534747L, 135006516L, 133851895L, 115169878L, 107349540L, 102531392L, 90354753L, 81195210L, 78077248L, 59128983L, 63025520L, 48129895L, 51304566L, 155270560L, 59373566L, 16571L)
####

class BamBins:
    '''
        count the reads number , not the coverage.
    '''
    def __init__(self,bamfilename,binsize=200,ReadFile=True):
        self.bins=defaultdict(list)
        self.bamfilename=bamfilename
        self.binsize=binsize
        self.samfile=pysam.Samfile(bamfilename,"rb")
        self.chrs=_references
        self.lengths=_lengths
        self.chromsize={}
        for i in range(len(self.chrs)): 
            self.chromsize[self.chrs[i]]=self.lengths[i]
        self.chr2tid={} 
        for i,chrs in enumerate(self.chrs):
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
                for bin_id,value in enumerate(self.bins[self.chr2tid[chrs]]):
                    self.bins[self.chr2tid[chrs]][bin_id] *= 100000000.0/float(self.binsize*self.mapped)
            print >>sys.stderr,"Normalizing ... ... done"
    def makevector(self):
        ''' 
        combine all bins into one vector
        '''
        vector=[]
        for chrs in self.chrs:
            vector+= self.bins[self.chr2tid[chrs]]
        return vector
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
    p=argparse.ArgumentParser(description='Example: %(prog)s -b Binsize -C control.bam --bam <H3K4me1.bam> <H3K4me2.bam> <H3k4me3.bam> <H3K36me3.bam> bamfile -O outputfile')
    p.add_argument('-V',"--version",action='version',version='%(prog)s 0.1.3')
    p.add_argument('--bin','-b',type=int,dest='Binsize',default=200, help="binsize default: %(default)d")
    p.add_argument('-C','--control',type=str,dest='control',metavar="input",help="control files [Option]")
    p.add_argument('--bam',type=str,dest='bams',metavar="Bamfile",nargs='+',help="treatment file.")
    p.add_argument('-O','-o','--output',type=str,dest="output",metavar='Output',default='cor_matrix.txt',help='output the correlation matrix into this file if you set. default: %(default)s')

    if len(sys.argv)==1:
        sys.exit(p.print_help())
    args = p.parse_args()
    return args

def process_name(path):
    '''
        eliminate the unneeded name
    '''
    names=path
    pattern='(BroadHistone|HaibTfbs|OpenChromChip|SydhTfbs)(H1hesc|K562)(\w+?)Aln'
    ss=re.findall(pattern,names)
    lab,cell,ChIP=ss[0]
    return cell+ChIP

def process_bamfile(bams):
	'''
		Check whether it is a bamfile or not.
	'''
	bamlist=[]
	for i in bams:
		if os.path.getsize(i) < 1000000:
			f=open(i,'r')
			for line in f:
				if not line.startswith('#'):
					line=line.strip().split()
					bamlist.append(line[0])
		else:
			bamlist.append(i)
	return bamlist

def out_matrix(matrix,names):
    '''
        output the matrix to the cor_matrix.txt file
    '''
    cor_matrix= matrix
    names     = names
    f=open(output,'w')
    f.write("title")
    for each in names:
        f.write( "\t" + process_name(os.path.basename(each)) )
    f.write("\n")
    for row,name in enumerate(names):
		f.write( process_name( os.path.basename(name) ) )
		for column in range(len(names)):
			f.write("\t")
			f.write(str(cor_matrix[row,column]))
		f.write("\n")
    f.close()

def test():
    global args,norm,output
    # parse the argument line
    args = parse_argument()
    norm = True
    bams = args.bams
    output = args.output

    print >>sys.stderr,"# xbam_cor_matrix.py 0.1.3 "
    if args.control:
		print >>sys.stderr,"# control file: ",os.path.basename(args.control)
    else:
		print >>sys.stderr,"# control file: None"

    bams=process_bamfile(bams)
    
    for i in bams:
        print >>sys.stderr, "# Bamfile: ", os.path.basename(i)
    print >>sys.stderr,"# Binsize: ", args.Binsize
    print >>sys.stderr,"# Normalization: ", norm
	
    if args.control:
		control_bins= BamBins(args.control,binsize=args.Binsize)
    bambins=defaultdict(list)
    bambins_vector=defaultdict(list)
    for x in range(len(bams)):
        # process the bamfile
        bambins[x]= BamBins(bams[x],binsize=args.Binsize)
        # ChIP-seq signal = RPKM(treatment)-RPKM(control)
        if args.control:
			bambins[x].correct(control_bins)
        bambins_vector[x]=bambins[x].makevector()
    # Computing the correlation matrix
    # empty zeros matrix
    print >>sys.stderr,"# Computing the correlation matrix "
    cor_matrix = np.zeros((len(bams),len(bams)))
    for i in range(len(bams)):
        cor_matrix[i,i]=1
    for i in range(len(bams)):
        for j in range(len(bams)):
            cor_matrix[i,j]= pearson(bambins_vector[i],bambins_vector[j])
            cor_matrix[j,i]= cor_matrix[i,j]
    # output
    print >>sys.stderr, "# Writting the",output
    out_matrix(cor_matrix,bams)
    print >>sys.stderr, "# All Done "



if __name__ == '__main__':
    test()

