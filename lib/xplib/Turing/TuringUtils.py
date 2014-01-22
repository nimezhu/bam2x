#!/usr/bin/env python
# Programmer : zhuxp
# Date: 
# Last-modified: 01-22-2014, 15:04:30 EST

from xplib.Annotation import *
from bitarray import bitarray
import TuringCodeBook as cb
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

def bitarray_to_rep(bitstr):
    s=""
    for i in range(0,len(bitstr),2):
        if bitstr[i] and (not bitstr[i+1]):
            s+="E"
        elif (not bitstr[i]) and bitstr[i+1]:
            s+="i"
        else:
            s+="_"
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
def TupleTuringFactory(bed12):
    import xplib.Turing.TuringCodeBook as cb
    from xplib.Turing import TuringCode, TuringGraph
    from xplib.Tuple.Bed12Tuple import CHROMSTART,CHROMEND,get_exons
    g=[]
    g.append((bed12[CHROMSTART],cb.ON))
    g.append((bed12[CHROMEND],cb.OFF))
    if len(bed12)==12:
        for i in get_exons(bed12):
            g.append((i[CHROMSTART],cb.BLOCKON))
            g.append((i[CHROMEND],cb.BLOCKOFF))
    g.sort()
    return g

def count_xor(path):
    s=0
    for i in range(0,len(path),2):
        if path[i]^path[i+1]: s+=1 # ^ is xor 
    return s        

def isCompatible(pathA,pathB):  ##pathA and B shoulc be same length
    #print "debug pathA",pathA
    #print "debug pathB",pathB
    for i in range(0,len(pathA),2):
        if not (pathA[i]&pathB[i] | pathA[i+1]&pathB[i+1]):
            return False
    return True
def isSelfIncompatible(pathA):
    #print "debug",pathA
    for i in range(0,len(pathA),2):
        if not pathA[i]:
            if not pathA[i+1]:
                return True
    return False

def isOverlap(pathA,pathB):   ## Fir paired end reads overlap is not start > stop, overlap definition is not the same as overlap reades
    for i in range(0,len(path),2):
        if (pathA[i]^pathA[i+1]) and (pathB[i]^pathB[i+1]):
            return True
    return False


