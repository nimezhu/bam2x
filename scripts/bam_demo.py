#!/usr/bin/python
# programmer : zhuxp
# usage:
import sys
from getopt import getopt
import pysam
def show_help():
    print >>sys.stderr,"bam2wig.py -b binsize file.bam > file.wig"
    print >>sys.stderr,"   convert bamfile into wigfile"
    print >>sys.stderr,"   library dependency: pysam"



def Main():
    opts,restlist = getopt(sys.argv[1:],"hb:r:",\
                        ["help","binsize=","region=","bam="])
    global Binsize,ifRegion,Region
    Binsize=200
    ifRegion=0
    for o, a in opts:
        if o in ("-h","--help"): show_help()
        if o in ("-b","--binsize"):
                Binsize=int(a)
        if o in ("-r","--region"):
                ifRegion=1
                Region=a
        if o in ("--bam"):
            ifBamFileList=1
            BamFileList=a
    bamlist=[]
    if ifBamFileList:
        f=open(BamFileList)
        for i in f:
            i=i.strip()
            bamlist.append(i)
    for i in restlist:
        a=i.split("\.")
        if a[-1]=="bam": 
            bamlist.append(i)

    for i in bamlist:
        print i
        s=pysam.Samfile(i,"rb")
        ss=0
        for chr in s.references:
            ss+=s.count(chr)
        print ss
    


    
if __name__=="__main__":
    Main()

