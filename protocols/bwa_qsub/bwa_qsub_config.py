#!/usr/bin/env python
# Programmer : zhuxp
# Date: 
# Last-modified: 02-08-2013, 00:39:28 EST
VERSION="0.1"
import os,sys,argparse
from xplib.Annotation import Bed
from xplib import TableIO
import signal
signal.signal(signal.SIGPIPE,signal.SIG_DFL)
import gzip
import time
def Template():
    pass

'''    
# working dir  e.g. $HOME/scratch/job....
WORKING_DIR=.

# where you put your samples
# sample name (fastq file name without the suffix )
# if your fastq file is  SampleA.fastq  
# then  SAMPLE=SampleA
#       FASTQ_SUFFIX=fastq
#       PAIRED_END=0
# if you hava paired end data
#     your have two files ;  e.g.  SampleA_1.fastq SampleA_2.fastq
# then  SAMPLE=SampleA
#       FASTQ_SUFFIX=fastq
#       PAIRED_END=1
SAMPLE_DIR=.
SAMPLE=small_test     
FASTQ_SUFFIX=fastq
PAIRED_END=1          #1 or 0


# bwa index directory
# usally should put genome bwa_index file in nearline
#  the genome index should in the bwa_index directory
BWA_INDEX=$HOME/nearline/bwa_index
GENOME=mm10      

UNIQ_MAP=1   # 1 or 0  if it is 1 , it will report another bam which only have uniq map. 

# where you put the bam output file
OUTPUT_DIR=output

# where you submit the job to hpcc or just generate the scripts

SUBMIT=0   # 0 or 1 

MAIL=nimezhu@gmail.com


#=========================== USING THE DEFAULT  IF YOU WANT ===============================
# where you put samtools and bwa
BIN=/home/wangj2/bin                      
BWA=$BIN/bwa
SAMTOOLS=$BIN/samtools
# bwa aln options 
CPU=8
BWA_ALN_OPTION="-Y"
'''


def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency : xplib')
    p.add_argument('-v','--version',action='version',version='%(prog)s '+VERSION)
    p.add_argument('-p','--paired_end',dest="b_paired_end",default=False,action="store_true",help="input file format")
    p.add_argument('-o','--output',dest="output",type=str,default="./",help="output dir")
    p.add_argument('-g','--genome',dest="genome",type=str,default="mm10",help="genome file")
    p.add_argument('-w','--working_dir',dest="working",type=str,default="./",help="working dir")
    
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()
def Main():
    '''
    IO TEMPLATE
    '''
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
    if args.input=="stdin":
        fin=sys.stdin
    else:
        try:
            x=args.input.split(".")
            if x[-1]=="gz":
                fin=gzip.open(args.input,"r")
            else:
                fin=open(args.input,"r")
        except IOError:
            print >>sys.stderr,"can't read file",args.input
            fin=sys.stdin
    '''
    END OF IO TEMPLATE 
    '''




    
if __name__=="__main__":
    Main()





