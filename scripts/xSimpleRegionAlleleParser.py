#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 09 Aug 2012 13:57:50

import os,sys,argparse
import types
import math
class SimpleRegionAllele:
    def __init__(self):
        self.string=""
        self.body_pvalue=1.0
        self.promoter_pvalue=1.0
        self.score=0.0
    def __str__(self):
        return self.string+"SCORE: "+str(self.score)+"\n//"
    def __cmp__(self,other):
        return  cmp(other.score,self.score)
    def update_score(self,MAX_BODY_SCORE=140,MAX_PRMT_SCORE=60):
        self.score=0.0
        if self.body_pvalue!=1:
            if self.body_pvalue==0:
                a=MAX_BODY_SCORE
            else:
                a=-7*math.log(self.body_pvalue,10)
                if a>MAX_BODY_SCORE:
                    a=MAX_BODY_SCORE
            self.score+=a
        if self.promoter_pvalue!=1:
            if self.promoter_pvalue==0:
                a=MAX_PRMT_SCORE
            else:
                a=-3*math.log(self.promoter_pvalue,10)
                if a>MAX_PRMT_SCORE:
                    a=MAX_PRMT_SCORE
            self.score+=a
def SimpleRegionAlleleIterator(handle):
    if type(handle)==type("s"):
        try:
            handle=open(handle,"r")
        except:
            raise ValueError("Can't open file %s"%handle)
    SRA=SimpleRegionAllele()
    for line in handle:
        line=line.strip()
        if len(line)==0: continue
        if line[0]=="#": continue
        if line=="//" or line=="// ": 
            SRA.update_score()
            yield SRA
            SRA=SimpleRegionAllele() #Reset
            continue
        SRA.string+=line+"\n"
        x=line.split(":")
        if x[0]=="BODY PVALUE":
            a=float(x[1])
            SRA.body_pvalue=a
        if x[0]=="PRMT PVALUE":
            a=float(x[1])
            if SRA.promoter_pvalue > a:
                SRA.promoter_pvalue=a
def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency : pysam')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('-i','--input',dest="input",type=str,help="input file")
    p.add_argument('-o','--output',dest="output",type=str,default="stdout",help="output file")
    
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
    f=SimpleRegionAlleleIterator(args.input)
    a=[]
    for r in f:
        a.append(r)
    a.sort()
    for r in a:
        print >>out,r


    
if __name__=="__main__":
    Main()



