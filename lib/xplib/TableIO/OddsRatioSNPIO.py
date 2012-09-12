#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 12 Sep 2012 10:31:21
import types
from xplib.Annotation import OddsRatioSNP

def OddsRatioSNPIterator(handle):
    if type(handle)==type("s"):
        try:
            handle=open(handle,"r")
        except:
            raise ValueError("Can't open file %s"%handle)
    for line in handle:
        line=line.strip()
        if len(line)==0: continue
        if line[0]=="#": continue
        x=line.split("\t")
        x=OddsRatioSNP(x)
        yield x
    return
    
