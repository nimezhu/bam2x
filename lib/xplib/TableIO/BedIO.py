#!/usr/bin/python
# Programmer : zhuxp
# Date:
# Last-modified: 02 Oct 2012 10:47:20
from xplib.Annotation import Bed
import types
def BedIterator(handle,**dict):
    '''
    Bed file iterator
    Usage Example:
        from xplib.TableIO.BedIO import BedIterator
        for bed in BedIterator(file or filename):
            print bed
    Wrapper in TableIO.parse(file or filename,"bed")
    Usage:
        for i in TableIO.parse(file or filename,"bed"):
            print bed
        
    Bed is xplib.Annotation.Bed object.
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
        b=Bed(x)
        yield b
    return
    
