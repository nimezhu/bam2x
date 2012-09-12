#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 10 Aug 2012 14:03:24

import os,sys,argparse
from xplib.Annotation import Bed
from xplib.Annotation import Utils
from xplib import TableIO
def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency : xplib')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('-i','--input',dest="input",type=str,help="input file")
    p.add_argument('-o','--output',dest="output",type=str,default="stdout",help="output file")
    p.add_argument('-b','--beds',dest="beds",action="store",default=[],help="feature annotation files",nargs="+")
    p.add_argument('-m',dest="m",action="store_true",default=False,help="print in table format")

    
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()
def Main():
    global args,out
    h={} #count the hits number
    hcode={} #count the hits number
    args=ParseArg()
    if args.output=="stdout":
        out=sys.stdout
    else:
        try:
            out=open(args.output,"w")
        except IOError:
            print >>sys.stderr,"can't open file ",args.output,"to write. Using stdout instead"
            out=sys.stdout
    datas=[]
    input_handle=TableIO.parse(args.input,"bed")
    beds=[]
    h[args.input]=0
    for i in input_handle:
        beds.append(i)
        h[args.input]+=1
    
    for f in args.beds:
        print "reading",f
        h[f]=0
        fi=TableIO.parse(f,"bed")
        '''
        data tuple structure
        (filename,data_binindex)
        usage example:
            bed=Bed([chr,start,end])
            for (file,data) in data:
                for feature in Utils.iterOverlapFeature(bed,data):
                    print feature
        '''
        datas.append((f,Utils.readIntoBinIndex(fi)))
    print >>out,"# Input:",args.input
    for i,f in enumerate(args.beds):
        print >>out,"# Feature No."+str(i+1),":",f

    for i in beds:

        print >>out,i,
        code=""
        if not args.m:
            print >>out,""
            print >>out,""
        for (f,d) in datas:
            s=""
            mark=0
            for feature in Utils.iterOverlapFeature(i,d):
                mark+=1
                s+=str(feature)+"\n"
            if mark>0:
                h[f]+=1
                code+="1"
                if not args.m:
                    print >>out,"OVERLAP FEATURES IN",f
                    print >>out,s
                else:
                    print >>out,"\t",1,
            else:
                code+="0"
                if not args.m:
                    print >>out,"NO OVERLAP FEATURES IN",f
                else:
                    print >>out,"\t",0,
        if hcode.has_key(code):
            hcode[code]+=1
        else:
            hcode[code]=1
        if args.m:
            print >>out,""
        else:
            print >>out,"//"
    for i,f in enumerate(args.beds):
        print >>out,"# ",h[f],"/",h[args.input],"overlap with No."+str(i+1),f
    for code in sorted(hcode.iterkeys()):
        print >>out,"# code:",code,"\t",hcode[code]

    
if __name__=="__main__":
    Main()




