#!/usr/bin/env python
# Programmer : zhuxp
# Date: 
# Last-modified: 02-12-2014, 15:40:13 EST
VERSION="0.1"
import pysam
import sys
import gzip

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

def open_output(output):
    out=None
    if output=="stdout":
        out=sys.stdout
    else:
        try:
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
    

def read_config(file,sep=":"):
    config={}
    for i in TableIO.parse(file,sep=sep):
        config[i[0]]=i[1]
    return config

