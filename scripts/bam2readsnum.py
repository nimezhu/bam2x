#!/usr/bin/python
# programmer : zhuxp
# usage:
import sys
from getopt import getopt
import pysam
def show_help():
    print >>sys.stderr,"bam2readsnum.py:  count bam reads number"
    print >>sys.stderr,"Version:"+Version+"\n"
    print >>sys.stderr,"Usage: bam2readsnum.py <options> file.bam > file.bamlist"
    print >>sys.stderr,"Example:  bam2readsnum.py file1.bam file2.bam > file.bamlist"
    print >>sys.stderr,"          bam2readsnum.py --bam old.bamlist > new.bamlist"
    print >>sys.stderr,"                       new bamlist example 1:"
    print >>sys.stderr,"                             /data/user/1.bam  10000000  # filename<tab>reads number<enter>"
    print >>sys.stderr,"                             /data/user/2.bam  11001100"
    print >>sys.stderr,"                       old bamlist example 2:"
    print >>sys.stderr,"                             /data/user/1.bam"
    print >>sys.stderr,"                             /data/user/2.bam"
    print >>sys.stderr,"Library dependency: pysam\n\n"
    print >>sys.stderr,"Options:"
    print >>sys.stderr,"   -h,--help          show this help message"
    print >>sys.stderr,"   --bam bamlistfile"
    print >>sys.stderr,"                      a file contains bam filenames and reads number<optional> in each line"
    print >>sys.stderr,"   "




def Main():
    global Version
    Version="0.1"
    if len(sys.argv)<2:
        show_help()
        exit
    opts,restlist = getopt(sys.argv[1:],"hb:",\
                        ["help","bam="])
    File_reads_number={}
    ifBamFileList=0


    for o, a in opts:                       
        if o in ("-h","--help"): 
            show_help()
            exit
        if o in ("--bam"):
            ifBamFileList=1
            BamFileList=a

    bamlist=[]
    if ifBamFileList:
        f=open(BamFileList)
        for i in f:
            i=i.strip()
            a=i.split("\t")
            bamlist.append(a[0])
            if len(a)>1:
                File_reads_number[a[0]]=long(a[1])   # using total reads number pre-computed in bamlistfile
            else:
                print >>sys.stderr,"counting ",a[0]  # counting reads number in bam , might be slow. 
                samfile=pysam.Samfile(a[0],"rb")
                ss=0
                for chr in samfile.references:
                    ss+=samfile.count(chr)
                File_reads_number[a[0]]=ss
                samfile.close()
    for i in restlist:
        i=i.strip()
        a=i.split('.')
        if a[-1]=="bam": 
            bamlist.append(i)
            samfile=pysam.Samfile(i,"rb")
            ss=0
            print >>sys.stderr,"counting ",i
            for chr in samfile.references:
                ss+=samfile.count(chr)
            File_reads_number[i]=ss
            samfile.close()
    for i in bamlist:
            print i+"\t"+str(File_reads_number[i])
    


    
if __name__=="__main__":
    Main()

