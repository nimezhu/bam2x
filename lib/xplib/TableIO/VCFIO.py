#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 02 Oct 2012 10:44:54
from xplib.Annotation import VCF
import types
def VCFIterator(handle,**kwargs):
    '''
    Wrapper in TableIO.parse(file or filename,"vcf")
    Usage:
        for i in TableIO.parse(file or filename,"vcf"):
            print i
        i is a Annotation.VCF class object.
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
        b=VCF(x)
        yield b
    return
    
