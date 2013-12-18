#!/usr/bin/env python
# Programmer : zhuxp
# Date: 
# Last-modified: 12-17-2013, 23:53:03 EST

from xplib.Annotation import *
from bitarray import bitarray
def method():
    print "aha"


def twobitarray_and(bitarrayA,bitarrayB):
    '''
    '11' code unknown
    '10' code block on
    '01' code block off

    usage:
    twobitarray_and(a,b)

    return true if they are compatible  
    which means there no '00' in each twobit and operation
    '''
    l1=len(bitarrayA)
    l2=len(bitarrayB)
    l=l1
    if l>l2: l=l2
    for i in range(0,l,2):
        a0=bitarrayA[i]
        a1=bitarrayA[i+1]
        b0=bitarrayB[i]
        b1=bitarrayB[i+1]
        if (a0!=b0 and a1!=b1): return False
    return True

def bitarray_and(bitarrayA,bitarrayB):
    l1=len(bitarrayA)
    l2=len(bitarrayB)
    if l1!=l2: return None
    C=bitarray(l1)
    for i in range(l1):
        C[i]=bitarrayA[i] and bitarrayB[i]
    return C

def bitarray_to_intron_number(bitarray):
    l=len(bitarray)
    s=0
    state=0
    for i in range(0,l,2):
        if bitarray[i] and (not bitarray[i+1]):
            state=1
        if (not bitarray[i]) and bitarray[i+1]:
            if state!=-1:
                s+=1
            state=-1
    return s    

def TuringFactory(bed12):
    import xplib.Turing.TuringCodeBook as cb
    from xplib.Turing import TuringCode, TuringGraph
    g=[]
    g.append(TuringCode(bed12.start,cb.ON))
    g.append(TuringCode(bed12.stop,cb.OFF))
    if isinstance(bed12,Bed12) or isinstance(bed12,GeneBed):
        for i in bed12.Exons():
            g.append(TuringCode(i.start,cb.BLOCKON))
            g.append(TuringCode(i.stop,cb.BLOCKOFF))
    g.sort()
    return TuringGraph(g)
