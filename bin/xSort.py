#!/usr/bin/env python
# Programmer : zhuxp
# Date: 
# Last-modified: 08 Oct 2012 17:52:58

import os,sys,argparse
from xplib.Annotation import Bed
from xplib import TableIO
import signal
signal.signal(signal.SIGPIPE,signal.SIG_DFL)
import gzip
'''
sort genome annotations
and print out
'''
def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -i file.bed -I bed -o file.sorted.bed', epilog='Library dependency : xplib')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('-i','--input',dest="input",default="stdin",type=str,help="input file DEFAULT: STDIN")
    p.add_argument('-I','--input_format',dest="format",default="bed",type=str,help="input file format")
    p.add_argument('-o','--output',dest="output",type=str,default="stdout",help="output file DEFAULT: STDOUT")
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
    if args.input=="stdin":
        fin=sys.stdin
    else:
        try:
            x=args.input.split(".")
            if x[-1]=="gz":
                fin=gzip.open(args.input,"r")
            else:
                fin=open(args.input,"r")
        except IOError:
            print >>sys.stderr,"can't read file",args.input
            fin=sys.stdin
    '''
    END OF IO TEMPLATE 
    '''
    a={}
    for i in TableIO.parse(fin,args.format):
        if not a.has_key(i.chr):
            a[i.chr]=[]
        a[i.chr].append(i)
    for chr in sorted(a.keys()):
        for i in a[chr]:
            print >>out,i
        




    
if __name__=="__main__":
    Main()





