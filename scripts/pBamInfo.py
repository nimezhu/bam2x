#!/usr/bin/env python
# programmer : zhuxp
# usage:
import sys
from getopt import getopt
import pysam
from multiprocessing import Process,Pool
def show_help():
    print >>sys.stderr,"pBamInfo.py:  show basic information on bam file"
    print >>sys.stderr,"Version:"+Version+"\n"
    print >>sys.stderr,"Usage: pBamInfo.py <options> file.bam > file.bamlist"
    print >>sys.stderr,"Example:  pBamInfo.py file1.bam file2.bam > file.bamlist"
    print >>sys.stderr,"          pBamInfo.py --bam old.bamlist > new.bamlist"
    print >>sys.stderr,"                       bamlist example 1:"
    print >>sys.stderr,"                             /data/user/1.bam"
    print >>sys.stderr,"                             /data/user/2.bam"
    print >>sys.stderr,"Library dependency: pysam\n\n"
    print >>sys.stderr,"Options:"
    print >>sys.stderr,"   -h,--help          show this help message"
    print >>sys.stderr,"   --bam bamlistfile"
    print >>sys.stderr,"   -n,                processor number"




def Main():
    global Version
    Version="0.1"
    if len(sys.argv)<2:
        show_help()
        exit
    opts,restlist = getopt(sys.argv[1:],"hb:n:",\
                        ["help","bam=","number="])
    File_reads_number={}
    ifBamFileList=0
    n=4
    for o, a in opts:                       
        if o in ("-h","--help"): 
            show_help()
            exit
        if o in ("--bam"):
            ifBamFileList=1
            BamFileList=a
        if o in ("-n","--number"):
            n=int(a)
    processes=[]
    bamlist=[]
    if ifBamFileList:
        f=open(BamFileList)
        for i in f:
            i=i.strip()
            a=i.split("\t")
            bamlist.append(a[0])
    for i in restlist:
        i=i.strip()
        a=i.split('.')
        if a[-1]=="bam": 
            bamlist.append(i)
    pool = Pool(processes=n)
    result=pool.map(count_read_number,bamlist)
    for i,x in zip(bamlist,result):
            print i+"\t"+str(x)
    

def count_read_number(filename):
    print >>sys.stderr,"counting ",filename
    samfile=pysam.Samfile(filename,"rb")
    ss=0
    for chr in samfile.references:
        ss+=samfile.count(chr)
    samfile.close()
    return ss
if __name__=="__main__":
    Main()

