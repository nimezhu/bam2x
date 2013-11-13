#!/usr/bin/env python
# Programmer : zhuxp
# Date: 
# Last-modified: 11-12-2013, 15:58:01 EST
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
