#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 07-12-2012, 14:00:51 CDT
import types

def SimpleIterator(handle):
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
        for i,y in enumerate(x):
            try:
              x[i]=float(x[i])
              if x[i]==int(x[i]):
                  x[i]=int(x[i])
            except:
              pass
        yield x
    return
    
