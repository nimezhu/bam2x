#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 30 Jun 2012 01:15:18
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
    
