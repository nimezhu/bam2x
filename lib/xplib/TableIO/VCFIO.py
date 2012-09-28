#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 28 Sep 2012 01:22:59
from xplib.Annotation import VCF
import types
def VCFIterator(handle):
    if type(handle)==type("s"):
        try:
            handle=open(handle,"r")
        except:
            raise ValueError("Can't open file %s"%handle)
    for line in handle:
        line=line.strip()
        if line[0]=="#": continue
        x=line.split("\t")
        b=VCF(x)
        yield b
    return
    
