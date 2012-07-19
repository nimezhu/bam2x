#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 12 Jul 2012 14:34:20

import os,sys,argparse
from xplib import TableIO

def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Extract possible allele specific SNPs candidate from LogRatio Output', epilog='Library dependency : pysam')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('-X','--chi2',dest="X2",type=float,default=6.0,help="chi-square score threshold, default %(default)f")
    p.add_argument('-i','--input',dest="input",type=str,required=True,help="input file")
    p.add_argument('-o','--output',dest="output",type=str,default="stdout",help="output file")
    p.add_argument('-c',dest="cat",default=False,action="store_true",help="mark the SNP number")
    
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(1)
    return p.parse_args()
def Main():
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
    try:
        infile=open(args.input,"r")
    except IOError:
        print >>sys.stderr,"can't open input file",args.input
    state=0
    last_chr=""
    last_pos=0
    idx=0
    for line in infile:
        line=line.strip()
        if line[0]=="#": 
            print >>out,line
            continue
        idx+=1

        x=line.split("\t")
        b=x[4].replace("( ","")
        b=b.replace(" )","")
        a=b.split(" ")
        for i,y in enumerate(a):
            a[i]=int(y)
        if (a[0]-a[1])*(a[2]-a[3])<0 and float(x[3])>args.X2:
            chr=x[0]
            pos=int(x[1])
            if chr != last_chr:
                print 
            if pos-last_pos > 100000:
                print
            last_chr=chr
            last_pos=pos
            if args.cat:
                print >>out,idx,"\t",line
            else:
                print >>out,line
    
if __name__=="__main__":
    Main()


