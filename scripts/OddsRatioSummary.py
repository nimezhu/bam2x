#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 12 Sep 2012 20:15:35

import os,sys,argparse
from xplib.Annotation import Bed
from xplib import TableIO
from xplib.Annotation import Utils
def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Automaticly Summary OddsRatioOutput. Summary the coverage distribution and cetromere affect. Output the regions that might be useful.  Example: %(prog)s -i file.OddsRatio.out', epilog='Library dependency : xplib')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('-i','--input',dest="input",type=str,help="input file")
    p.add_argument('-o','--output',dest="output",type=str,default="stdout",help="output file")
    p.add_argument('-c','--centromere',dest="centromere",type=str,help="centromere bed file",required=True)
    p.add_argument('-t','--threshold',dest="t",type=int,default=1000,help="coverage threshold")
    p.add_argument('-X','--chi2',dest="chi2",type=float,default=10.83,help="chi square threshold, default: %(default)f   equals to pvalue 0.001 when freedom is 1")
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()
def sum(x):
    s=0
    for i in x:
        s+=i
    return s
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
    h=[[0,0],[0,0]]
    l=[[0,0],[0,0]]
    d=Utils.readIntoBinIndex(TableIO.parse(args.centromere,"bed"))
    print >>out,"# Coverage Threshold: ",args.t
    print >>out,"# Chi2 Threshold:",args.chi2
    for i in TableIO.parse(args.input,"oddsratiosnp"):
        mark=0
        for j in Utils.iterOverlapFeature(i,d):
            mark+=1
        if mark>1: mark=1
        if sum(i.A_nt_dis) > args.t and sum(i.B_nt_dis) > args.t:
            if i.odds_ratio > args.chi2:
                h[mark][1]+=1
            else:
                h[mark][0]+=1
            print >>out,i,"\tHigh\t",mark
        else:
            print >>out,i,"\tLow\t",mark
            if i.odds_ratio > args.chi2:
                l[mark][1]+=1
            else:
                l[mark][0]+=1
    print >>out,"# HighOddsRatio:",h[0][1]+l[0][1]+h[1][1]+l[1][1]
    print >>out,"# LowOddsRatio:",h[0][0]+l[0][0]+h[1][0]+l[1][0]
    print >>out,"#"
    print >>out,"# HighCoverage:",sum(h[1])+sum(h[0])
    print >>out,"# LowCoverage :",sum(l[1])+sum(l[0])
    print >>out,"#"
    print >>out,"# HighCoverage, HighOddsRatio",h[1][1]+h[0][1]
    print >>out,"# HighCoverage, LowOddsRatio",h[1][0]+h[0][0]
    print >>out,"#"
    print >>out,"# HighCoverage, InCentromere",sum(h[1])
    print >>out,"# HighCoverage, NotInCentromere",sum(l[1])
    print >>out,"#"

    print >>out,"# HighCoverage, HighOddsRatio, InCentromere",h[1][1]
    print >>out,"# HighCoverage, HighOddsRatio, NotInCentromere",h[0][1]
    print >>out,"# HighCoverage, LowOddsRatio,  InCentromere",h[1][0]
    print >>out,"# HighCoverage, LowOddsRatio, NotInCentromere",h[0][0]
    print >>out,"# LowCoverage, HighOddsRatio, InCentromere",l[1][1]
    print >>out,"# LowCoverage, HighOddsRatio, NotInCentromere",l[0][1]
    print >>out,"# LowCoverage, LowOddsRatio, InCentromere",l[1][0]
    print >>out,"# LowCoverage, LowOddsRatio, NotInCentromere",l[0][0]
    

if __name__=="__main__":
    Main()




