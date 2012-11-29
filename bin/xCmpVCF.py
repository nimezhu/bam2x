#!/usr/bin/env python
# Programmer : zhuxp
# Date: 
# Last-modified: 29 Nov 2012 16:28:38

import os,sys,argparse
from xplib import TableIO
from xplib import DBI
import signal
signal.signal(signal.SIGPIPE,signal.SIG_DFL)
import gzip

def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency : xplib')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('-o','--output',dest="output",type=str,default="stdout",help="output file DEFAULT: STDOUT")
    #p.add_argument('-A','--BAM_A',dest="BAM_A",nargs="*",required=True,help="BAM A")
    p.add_argument('-a','--VCF_A',dest="VCF_A",type=str,required=True,help="VCF A")
    #p.add_argument('-B','--BAM_B',dest="BAM_B",nargs="*",required=True,help="BAM B")
    p.add_argument('-b','--VCF_B',dest="VCF_B",type=str,required=True,help="VCF B")
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
    '''
    END OF IO TEMPLATE 
    '''
    print >>out,"# QUERY (VCF A):",args.VCF_A
    print >>out,"# DATA  (VCF B):",args.VCF_B
    print >>out,"# A11 VCF in A and B, and alt nt is the same : A VCF entry" 
    print >>out,"# B11 VCF in A and B, and alt nt is the same : B VCF entry" 
    print >>out,"# A12 VCF position in A and B, and alt nt is not the same : A VCF entry" 
    print >>out,"# B12 VCF position in A and B, and alt nt is not the same : B VCF entry" 
    print >>out,"# A10 VCF, only exists in A"
    print >>out,"# B01 VCF, only exists in B"
    print >>sys.stderr,"Initialize data: reading ",args.VCF_A
    VCF_A_DBI=DBI.init(args.VCF_A,"vcf")
    print >>sys.stderr,"Initialize data: reading ",args.VCF_B
    VCF_B_DBI=DBI.init(args.VCF_B,"vcf")
    
    A11=0
    A12=0
    A10=0
    B01=0
    
    i0=0
    print >>sys.stderr,"Query ",args.VCF_A
    for (x,i) in enumerate(VCF_A_DBI):
        if x%1000==0: print >>sys.stderr,x," entries\r",
        flag=0
        hit=None
        for j in VCF_B_DBI.query(i):
            if i==j: 
                hit=j
                flag=1
                continue
            else:
                hit=j
                flag=2
        if flag==1:
            print >>out,"A11_%08d\t"%A11,i
            print >>out,"B11_%08d\t"%A11,hit
            print >>out,""
            A11+=1
        elif flag==2:
            print >>out,"A12_%08d\t"%A12,i
            print >>out,"B12_%08d\t"%A12,hit
            print >>out,""
            A12+=1
        else:
            print >>out,"A10_%08d\t"%A10,i
            print >>out,""
            print >>out,""
            A10+=1
    print >>sys.stderr,"Query ",args.VCF_B
    for (x,i) in enumerate(VCF_B_DBI):
        if x%1000==0: print >>sys.stderr,x," entries\r",
        flag=0
        for j in VCF_A_DBI.query(i):
            flag=1
        if flag==0:
            print >>out,"B01_%08d\t"%B01,i
            print >>out,""
            print >>out,""
            B01+=1
    print >>out,"# [AB]11 number:",A11
    print >>out,"# [AB]12 number:",A12
    print >>out,"# A10 number:",A10
    print >>out,"# B01 number:",B01





    
if __name__=="__main__":
    Main()





