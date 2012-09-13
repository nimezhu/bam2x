#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 11 Sep 2012 09:15:55
'''
convert OddsRatio  Out file  to wig
OddsRatio Out file is 0-index
Wig is 1-index

'''
import os,sys,argparse
from xplib import TableIO
def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -i file.OddsRatio.out -o file.wig', epilog='Library dependency : xplib')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('-i','--input',dest="input",type=str,help="input file")
    p.add_argument('-o','--output',dest="output",type=str,default="stdout",help="output file")
    p.add_argument('-p','--prefix',dest="prefix",type=str,default="",help="prefix of names")
    
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
    lastchr=""
    x=args.input.split("/")
    x[-1]=x[-1].replace(".OddsRatio.out","")
    x[-1]=x[-1].replace(".LogR.out","")
    name=x[-1]
    color='50,150,255'

    if name.find('K9me3')>-1 or name.find('k9me3')>-1:
        color='255,0,0'
    elif name.find('K4me3')>-1 or name.find('k4me3')>-1:
        color='0,255,255'
    elif name.find('K4me1')>-1 or name.find('k4me1')>-1:
        color='0,255,127'
    elif name.find('K36me3')>-1 or name.find('k36me3')>-1:
        color='150,255,0'
    elif name.find('K27me3')>-1 or name.find('k27me3')>-1:
        color='255,170,200'
    elif name.find('K27ac')>-1 or name.find('k27ac')>-1:
        color='50,100,255'
    
    print >>out,"browser position chr1:16890000-16990000"
    print >>out,"browser hide all"
    print >>out,"track type=wiggle_0 name=\"OddsRatio "+args.prefix+" "+x[-1]+"\" description=\"OddsRatio of "+args.prefix+" "+x[-1]+"\" visibility=full autoScale=off viewLimits=0.0:15.0 color="+color+" yLineMark=5.00 yLineOnOff=on priority=10"
    for a in TableIO.parse(args.input,"simple"):
        a[0]=a[0].strip()
        if lastchr!=a[0]:
            print >>out,"variableStep chrom="+a[0]
            lastchr=a[0]
        print >>out,a[1]+1,a[3]



    
if __name__=="__main__":
    Main()




