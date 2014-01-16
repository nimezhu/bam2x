#!/usr/bin/env python
# Programmer : zhuxp
# Date: 
# Last-modified: 01-16-2014, 13:25:06 EST
VERSION="0.1"
import os,sys,argparse
import gzip
import time
import pysam
import xplib.Tools as Tools
from xplib import TableIO
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
    if Tools.guess_format(file)=="bam" and mode=="r":
        return pysam.Samfile(file,"rb")
    if mode=="w":
        return open_output(file)
    if mode=="r":
        return open_input(file)
    return None
    
def formatter(metabed,attr):
    '''
    format output of MetaBed3
    example:
        print formatter(metabed,["chr","start","stop"])
    '''
    s=""
    for i in attr:
        s+=str(metabed.__dict__[i])+"\t"
    s.strip('\t')
    return s

def read_config(file,sep=":"):
    config={}
    for i in TableIO.parse(file,sep=sep):
        config[i[0]]=i[1]
    return config

