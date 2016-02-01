#!/usr/bin/env python
# Programmer : zhuxp
# Date: 
# Last-modified: 01-31-2016, 10:36:11 EST
from __future__ import print_function
VERSION="0.1"
import os,sys,argparse
from bam2x import TableIO,Tools
from bam2x import IO,DBI
import time
'''
report bed in selected regions ( bed6 tabix format ).

'''

def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency : bam2x')
    p.add_argument('-v','--version',action='version',version='%(prog)s '+VERSION)
    p.add_argument('-i','--input',dest="input",default="stdin",type=str,help="input file DEFAULT: STDIN")
    p.add_argument('-I','--input_format',dest="format",default="bed6",type=str,help="input file format")
    p.add_argument('-o','--output',dest="output",type=str,default="stdout",help="output file DEFAULT: STDOUT")
    p.add_argument('-b','--bed',dest="bed",type=str,help="bed tabix region")
    #p.add_argument('-n','--num_cpus',dest="num_cpus",type=int,default=4,help="number of cpus DEFAULT: %(default)i")
    if len(sys.argv)==1:
        print(p.print_help(),file=sys.stderr)
        exit(0)
    return p.parse_args()
def Main():
    '''
    IO TEMPLATE
    '''
    global args,out
    args=ParseArg()
    fin=IO.fopen(args.input,"r")
    out=IO.fopen(args.output,"w")
    '''
    END OF IO TEMPLATE 
    '''
    db=DBI.init(args.bed,"tabix",cls="bed3")
    for i in TableIO.parse(fin,args.format):
        signal=0;
        for j in db.query(i):
            signal=1;
            break
        if signal==1:
            print(i,file=out)




    
if __name__=="__main__":
    Main()





