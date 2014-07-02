#!/usr/bin/env python
# Programmer : zhuxp
# Date: 
# Last-modified: 07-02-2014, 11:38:05 EDT
from bam2x import __version__ as VERSION
import pysam
import sys
import gzip
import argparse
import os.path
import logging
suffixToFormat={
    'fa':'fasta',
    'fq':'fastq',
    'genetab':'genebed',
    'bw':'bigwig',
    'tab':'genebed',
    '2bit':'genome'

}

def suffix(string):
    x=string.split('.')
    return x[-1]
from string import lower

def guess_format(string):
    
    x=string.split(".")
    suffix=lower(x[-1])
    if x[-1]=="gz":
        suffix=lower(x[-2])
    if suffixToFormat.has_key(suffix):
        return suffixToFormat[suffix]
    else:
        return suffix

def get_col_number(fn,sep="\t"):
    f=fopen(fn,"r")
    a=f.next().split(sep)
    retv=len(a)
    f.close()
    return retv

def open_output(output):
    out=None
    if output=="stdout":
        out=sys.stdout
    else:
        try:
            if os.path.isfile(output):
                i=1;
                newname=None
                while(True):
                    name,ext=os.path.splitext(output)
                    newname=name+"("+str(i)+")"+ext
                    if not os.path.isfile(newname):
                        break
                    i+=1
                output=newname
                logging.warn("output file exists, automatically rename the output file to "+output)
            out=open(output,"w")
        except IOError:
            print >>sys.stderr,"can't open file ",output,"to write. Using stdout instead"
            out=sys.stdout
    return out

def open_input(input):
    fin=None
    if input=="stdin":
        fin=sys.stdin
    else:
        try:
            x=input.split(".")
            if x[-1]=="gz":
                fin=gzip.open(input,"r")
            else:
                fin=open(input,"r")
        except IOError:
            print >>sys.stderr,"can't read file",input
            fin=sys.stdin
    return fin
    
def fopen(file,mode="r",**kwargs):
    '''
    '''
    if guess_format(file)=="bam" and mode=="r":
        return pysam.Samfile(file,"rb")
    if mode=="w":
        return open_output(file)
    if mode=="r":
        return open_input(file)
    return None
    
def get_col_num(filehandle,sep="\t"):
    try:
        filehandle.seek(0)
        for i in filehandle:
            if i[0]=="#": continue
            l=i.split(sep)
            filehandle.seek(0)
            return len(l)
    except:
        raise

def read_config(file,sep=":"):
    config={}
    for i in TableIO.parse(file,sep=sep):
        config[i[0]]=i[1]
    return config













def parser_factory(**dict):
    ''' This Function Parse the Argument '''
    group=argparse.ArgumentParser( **dict)
    group.add_argument('-v','--version',action='version',version='%(prog)s '+VERSION)
    group.add_argument('-i','--input',dest='input',default='stdin',type=str,help="input file Default: %(default)s")
    group.add_argument('-o','--output',dest='output',default='stdout',type=str,help="output file Default: %(default)s")
    return group
def mparser_factory(**dict):
    ''' This Function Parse the Argument '''
    group=argparse.ArgumentParser( **dict)
    group.add_argument('-v','--version',action='version',version='%(prog)s '+VERSION)
    group.add_argument('-i','--input',dest='input',default='stdin',type=str,help="input file Default: %(default)s")
    group.add_argument('-o','--output',dest='output',default='stdout',type=str,help="output file Default: %(default)s")
    group.add_argument('-n','--num_cpus',dest='num_cpus',default=4,type=int,help="number of cpus, Default: %(default)i")
    return group

    
if __name__=="__main__":
    pass








