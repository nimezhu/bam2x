#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 26 Jul 2012 11:12:32

import os,sys,argparse

def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -i file.tab -s file.namelist -o output.tab', epilog='extract a subset from table.')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('-i','--input',dest="input",type=str,help="input table")
    p.add_argument('-s','--subset',dest="subset",type=str,help="file with IDs that you want to extract from input file")
    p.add_argument('-o','--output',dest="output",type=str,default="stdout",help="output file, default: stdout")
    p.add_argument('-k','--ik',dest="ik",type=int,default=1,help="the column of ID in input file [default: %(default)i") 
    p.add_argument('--sk',dest="sk",type=int,default=1,help="the column of ID in subset namelist file [default %(default)i")
    
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
    h={}
    subset=open(args.subset,"r")
    for line in subset:
        line=line.strip()
        if len(line)==0: continue
        if line[0]=="#": 
            continue
        x=line.split("\t")
        h[x[args.sk-1].strip()]=1
    
    table=open(args.input,"r")
    for line in table:
        line=line.strip()
        if len(line)==0: continue
        if line[0]=="#": 
            print >>out,line
            continue
        x=line.split("\t")
        if h.has_key(x[args.ik-1].strip()):
            print >>out,line
            h[x[args.ik-1].strip()]=2
    for i in h.keys():
        if h[i]==1:
            print >>sys.stderr,"can't find ",i,"in file"

        
         


    
if __name__=="__main__":
    Main()



