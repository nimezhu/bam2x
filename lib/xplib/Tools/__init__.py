# Programmer : zhuxp
# Date:  Sep 2012
# Last-modified: 06-12-2013, 15:01:47 EDT

hNtToNum={'a':0,'A':0,
          'c':1,'C':1,
          'g':2,'G':2,
          't':3,'T':3
         }
Nt=['A','C','G','T']

    
def rc(seq):
   comps = {'A':"T", 'C':"G", 'G':"C", 'T':"A",
           'B':"V", 'D':"H", 'H':"D", 'K':"M",
           'M':"K", 'R':"Y", 'V':"B", 'Y':"R",
           'W':'W', 'N':'N', 'S':'S'}
   return ''.join([comps[x] for x in seq.upper()[::-1]])
def shuffle(seq):
   import random
   a=list(seq)
   random.shuffle(a)
   return "".join(a)

def seq_wrapper(seq,width=60):
    s=""
    seqlen=len(seq)
    for i in range(0,seqlen,width):
        stop=i+width
        if stop>seqlen:stop=seqlen
        s+=seq[i:stop]+"\n"
    return s
def distance(A,B):
    if A.chr!=B.chr: return None
    if overlap(A,B): return 0
    m=abs(A.start-B.start)
    if( m > abs(A.start-B.stop)): m = abs(A.start-B.stop)
    if( m > abs(A.stop-B.stop)): m = abs(A.stop-B.stop)
    if( m > abs(A.stop-B.start)): m = abs(A.stop-B.start)
    return m
def translate_coordinate(A,B):
    '''
    translate B's coordiante based on A
    '''
    if A.chr!=B.chr: return None
    if A.strand=="+" or A.strand==".":
        return (B.start-A.start,B.stop-A.start,B.strand)
    if A.strand=="-":
        strand="."
        if B.strand=="-":strand="+"
        if B.strand=="+":strand="-"
        return (A.stop-B.stop,A.stop-B.start,strand)

def overlap(A,B):
    '''
    if A is overlapping with B.
    A and B are ? extends Bed class.
    '''
    if(A.chr != B.chr) : return 0
    if (A.stop < B.start) : return 0
    if (B.stop < A.start) : return 0
    return 1
