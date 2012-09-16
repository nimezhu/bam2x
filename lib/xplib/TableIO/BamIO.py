#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 16 Sep 2012 14:07:50
import types
import pysam
from xplib.Annotation import Bed

def BamIterator(filename):
    f=pysam.Samfile(filename,"rb")
    for i in f:
        yield i
def SamIterator(filename):
    f=pysam.Samfile(filename,"r")
    for i in f:
        yield i
def BamToBedIterator(filename):
    f=pysam.Samfile(filename,"rb")
    for i in f:
        if i.tid<0:continue
        strand="+"
        if i.is_reverse:
            strand="-"
        score=i.mapq
        bed=Bed([f.references[i.tid],i.pos,i.aend,i.qname,score,strand])
        yield bed
