#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 30 Jun 2012 00:16:38
from xplib.Annotation import Bed
import types
def BedIterator(handle):
    if type(handle)==type("s"):
        try:
            handle=open(handle,"r")
        except:
            raise ValueError("Can't open file %s"%handle)
    for line in handle:
        line=line.strip()
        if line[0]=="#": continue
        x=line.split("\t")
        b=Bed(x)
        yield b
    return
    
def BedWriter(handle,bed):
     if type(handle)==type("s"):
        try:
            handle=open(handle,"w")
        except:
            raise ValueError("Can't open file %s"%handle)
        handle.write(bed+"\n")


    

