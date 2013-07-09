#!/usr/bin/python
# Wei Guifeng, <guifengwei@gmail.com>

import sys, math, re, os, string
import scipy.stats
import numpy as np
from itertools import izip

def pearson(x=None, y=None):
    ''' 
        Calculate the pearson correlation coefficient 
    '''
    if isinstance(x,list) and isinstance(y,list):
        if len(x) == len(y):
            coeff,pvaule=scipy.stats.pearsonr(x,y)
            return coeff
        else:
            print >>sys.stderr,"The length between",x,"and",y," is not equal."
    else:
        print >>sys.stderr,"The type of",x,"and",y,"must be list!"

def spearman(x=None, y=None):
    ''' 
        Calculate the spearman correlation coefficient 
    '''
    if isinstance(x,list) and isinstance(y,list):
        if len(x) == len(y):
            coeff,pvaule=scipy.stats.spearmanr(x,y)
            return coeff
        else:
            print >>sys.stderr,"The length between",x,"and",y," is not equal"
    else:
        print >>sys.stderr,"The type of",x,"and",y,"must be list!"

def median(numbers):
    '''Return the median of the list numbers.'''
    # sort the list of numbers and take the middle element
    n=len(numbers)
    copy=numbers[:]
    copy.sort()
    if n & 1: # There is an odd number of elements
        return copy[n / 2]
    else:
        return (copy[n / 2 - 1] + copy[n / 2]) / 2

def mode(numbers):
    '''Return the mode of the list numbers '''
    #Find the value that occurs the most frequently in a data set
    freq={}
    for i in range(len(numbers)):
        try:
            freq[number[i]] += 1
        except KeyError:
            freq[numbers[i]] = 1
        max = 0
        mode = None
        for k, v in freq.iteritems():
            if v > max:
                max = v
                mode = k
        return mode

##################################################
#  Tissue Specificity measurement
#   Usage: a=[1,2,3,4,5,6]
#          socre, whichone = SpecificityScore(a)
#
def ShannonEntropy(e):
    ''' Return the entropy of the list numbers '''
    # entropy_score = (1..N) -Pi*log2(Pi)
    if isinstance(e,list):
        if sum(e)>0 and min(e)>=0:
            probs = [float(i)/sum(e) for i in e]
            entropy = -sum([ i*math.log(i,10) for i in probs if i > 0 ])
            return entropy
        else:
            return 0
    else:
        print >>sys.stderr, 'ShannonEntropy: the input is not a list'

def JSdist(p,q):
    '''JSdist is the square root of the JS divergence'''
    if isinstance(p,list) and isinstance(q,list):
        JS_div = ShannonEntropy([0.5*(i+j) for i,j in izip(p,q)])-0.5*(ShannonEntropy(p) + ShannonEntropy(q))
        if JS_div >=0:
            return math.sqrt(JS_div)
        else:
            return Inf
    else:
        print >>sys.stderr, 'JSdist: the input is not a list'

def SpecificityScore(e):
    ''' return the tissue-specific score (see Moran N. Cabili et al. 2011 Genes & Dev for details)
        as well as the which one in the list
        usage:
            Score, which_one = SpecificityScore(e)
    '''
    if isinstance(e, list):
        size = len(e)
        # normalize the expression vector to density vector
        # add pseudo-count 1 to raw FPKM and normalized density vector
        densityV =  [ math.log(i+1, 2) for i in e ]
        ss = sum(densityV)
        if ss == 0:
            # if the list are full of zero, the specificity score is 0 and the output of which_one is 0.
            return (0, 0)
        newVector = [ 1.0 * i /ss for i in densityV ]
        # Expr is numpy.array and diagonal
        Expr = np.diag([1] * size)
        js_distance_list = []
        for i in Expr:
            js_distance_list.append( 1 - float(JSdist(newVector,list(i))) )
        score = max(js_distance_list)
        which = js_distance_list.index(score)
        return score, which+1
    else:
        print >>sys.stderr, "SpecificityScore: the input is not a list" 

##################################################

if __name__ == '__main__':
    header = ["Gene1", "Gene2", "Gene3", "Gene4", "Gene5", "Gene6", "Gene7"]
    a=[0,0,0,0,0,0,0]
    print a
    (score, which_one) = SpecificityScore(a)
    if which_one == 0:
        print score, "NONE"
        pass
    else:
        print score, header[which_one-1]

