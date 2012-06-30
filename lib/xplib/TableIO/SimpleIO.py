#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 06-30-2012, 13:47:38 CDT
import types

def SimpleIterator(handle):
    if type(handle)==type("s"):
        try:
            handle=open(handle,"r")
        except:
            raise ValueError("Can't open file %s"%handle)
    for line in handle:
        line=line.strip()
        if line[0]=="#": continue
        x=line.split("\t")
        for i,y in enumerate(x):
            try:
              x[i]=float(x[i])
              if x[i]==int(x[i]):
                  x[i]=int(x[i])
            except:
              pass
        yield x
    return
    
