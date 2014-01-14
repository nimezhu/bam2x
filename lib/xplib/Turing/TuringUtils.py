#!/usr/bin/env python
# Programmer : zhuxp
# Date: 
# Last-modified: 01-14-2014, 14:28:59 EST

from xplib.Annotation import *
from bitarray import bitarray
import TuringCodeBook as cb
from xplib.Turing import  TuringCode
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

def count_xor(path):
    s=0
    for i in range(0,len(path),2):
        if path[i]^path[i+1]: s+=1 # ^ is xor 
    return s        

def isCompatible(pathA,pathB):  ##pathA and B shoulc be same length
    for i in range(0,len(pathA),2):
        if not (pathA[i]&pathB[i] | pathA[i+1]&pathB[i+1]):
            return False
    return True
def isSelfIncompatible(pathA):
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


import struct
import array
import collections
def write_int_array(ints,f):
    sign=0
    if type(f)==type("str"):
        try:
            f=open(f,"wb")
            sign=1
        except:
            print >>sys.stderr,"can't open file to write"
            exit(0)
    for i in ints:
        f.write(struct.pack('i',i))
    if sign==1:
        f.close()


def read_int_array(f):
    '''
    read any iterable int array or read from file
    '''
    sign=0
    if type(f)==type("str"):
        try:
            f=open(f,"rb")
            sign=1
        except:
            print >>sys.stderr,"can't open file to read"
            exit(0)
    if isinstance(f,file):
        while True:
            a = array.array('i')
            a.fromstring(f.read(4096))
            if not a:
                break
            for i in a:
                yield i
        if sign==1:
            f.close()
        raise StopIteration
    elif isinstance(f,collections.Iterable):
        for i in f:
            yield i
        raise StopIteration

def write_turing_array(turings,f,mode="n"):
    '''
    HEADERCODE: normal 3721
    HEADERCODE: compress 9981
    '''
    if mode=="n":
        header=cb.NORMAL_CODES_ARRAY_HEADER
    elif mode=="c":
        header=cb.COMPRESS_CODES_ARRAY_HEADER

    sign=0
    if type(f)==type("str"):
        try:
            f=open(f,"wb")
            sign=1
        except:
            print >>sys.stderr,"can't open file to write"
            exit(0)
    f.write(struct.pack("i",header))
    if mode=="n":
        for i in turings:
            f.write(struct.pack('i',i.pos))
            f.write(struct.pack('i',i.code))
    elif mode=="c":
        last=turings.next()
        number=1
        for i in turings:
            if i.pos==last.pos and i.code==last.code:
                number+=1
            else:
                f.write(struct.pack('i',last.pos))
                f.write(struct.pack('i',last.code))
                f.write(struct.pack('i',number))
                last=i
                number=1
        f.write(struct.pack('i',last.pos))
        f.write(struct.pack('i',last.code))
        f.write(struct.pack('i',number))
    if sign==1:
        f.close()



def read_turing_array(f,mode="n"):
    '''
    read any turing encode array 
    parse it into NORMAL MODE or COMPRESS MODE
    '''
    r=read_int_array(f)
    header=r.next()
    if mode=="n":
        if header==cb.NORMAL_CODES_ARRAY_HEADER:
            while True:
                try:
                    pos=r.next()
                    code=r.next()
                    yield TuringCode(pos,code)
                except StopIteration:
                    raise StopIteration
        elif header==cb.COMPRESS_CODES_ARRAY_HEADER:
            while True:
                try:
                    pos=r.next()
                    code=r.next()
                    number=r.next()
                    for i in range(number):
                        yield TuringCode(pos,code)
                except StopIteration:
                    raise StopIteration
        raise StopIteration
    elif mode=="c":
        last_pos=None
        last_code=None
        number=None
        if header==cb.COMPRESS_CODES_ARRAY_HEADER:
            while True:
                try:
                    pos=r.next()
                    code=r.next()
                    number=r.next()
                    yield TuringCode(pos,code),number
                except:
                    raise StopIteration
        elif header==cb.NORMAL_CODES_ARRAY_HEADER:
            try:
                last_pos=r.next()
                last_code=r.next()
                number=1
            except StopIteration:
                raise StopIteration
            while True:
                try:
                    pos=r.next()
                    code=r.next()
                    if pos==last_pos and code==last_code:
                        number+=1
                    else:
                        yield TuringCode(last_pos,last_code),number
                        number=1
                        last_pos=pos
                        last_code=code
                except StopIteration:
                    break
            yield TuringCode(last_pos,last_code),number
            raise StopIteration
        raise StopIteration



