#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 07-17-2012, 13:28:40 CDT

import os,sys,argparse
import pysam
def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency : pysam')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('-i','--input',dest="input",type=str,help="bamfile")
    p.add_argument('-o','--output',dest="output",type=str,default="stdout",help="output bed file")
    
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)
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
    bam=pysam.Samfile(args.input,"rb")
    for x,i in enumerate(bam):
        if x%100000==0: print >>sys.stderr,"converted",x,"reads\r",
        if i.tid<0: continue
        print >>out,bam.references[i.tid]+"\t"+str(i.pos)+"\t"+str(i.pos+i.alen)




    
if __name__=="__main__":
    Main()



