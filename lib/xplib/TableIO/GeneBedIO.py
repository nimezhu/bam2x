#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 07-02-2012, 23:20:31 CDT
from xplib.Annotation import GeneBed
import types
def GeneBedIterator(handle):
    '''
    '''
    if type(handle)==type("s"):
        try:
            handle=open(handle,"r")
        except:
            raise ValueError("Can't open file %s"%handle)
    for line in handle:
        line=line.strip()
        if line[0]=="#": continue
        x=line.split("\t")
        b=GeneBed(x)
        yield b
    return
    
