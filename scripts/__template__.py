#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 29 Jun 2012 13:54:48
import os,sys,argparse
import pysam
def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -i file.bam -o file.bed', epilog='Library dependency : pysam')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('--input','-i',type=str,dest='Bamfile',help="the input alignment bamfile")
    p.add_argument('--output','-o',type=str,dest='Output',default="stdout", help="the output bed file  default: %(default)s")
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)

    return p.parse_args()









    

def Main():
    global args,Table
    args=ParseArg()
    if args.Output=="stdout":
        out=sys.stdout
    else:
        out=open(args.Output,"w")
    samfile=pysam.Samfile(args.Bamfile,"rb")
    chrs=samfile.references
    lengths=samfile.lengths
    Table=[]


    
if __name__=="__main__":
    Main()


