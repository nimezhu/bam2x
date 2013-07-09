#!/usr/bin/python
# Wei Guifeng, <guifengwei@gmail.com>
# -- coding: utf-8 --
''' Usage: classify_CGI(bed) == 'H'
    This script can be used to assign the promoter as HCP,ICP,LCP.
'''
import sys, os, string, re
from zbed import Bed

_ratio_Cutoff_high = 0.6
_ratio_Cutoff_low  = 0.4
_gc_Cutoff_high    = 0.55

def getSequence(bed):
    ''' get the sequence from the hg19 genomic coordinates '''
    seq=bed.getSeq()
    return seq

def GC_content(Sequence):
    ''' Calculate the GC percent in the DNA sequence 
        Number of C or G was divied by the length of this sequence
    '''
    seq = Sequence.upper()
    # GC number
    GC = seq.count('G') + seq.count('C')
    return float(GC)/len(seq)

def Obs_Exp_CpG_ratio(Sequence):
    '''
        Obs/Exp CpG = Number of CpG * N / (Number of C * Number of G)
    '''
    seq = Sequence.upper()
    CpG_number = len( re.findall(r'CG',seq) )
    C_number = seq.count('C')
    G_number = seq.count('G')
    try:
        ratio = float(CpG_number) * len(seq) / (C_number * G_number)
    except ZeroDivisionError:
        ratio = 0.0
    return ratio

def classify_CGI(bed):
    ''' input: bed region
        output: high-CpG ==> 'H', Low-CpG ==> 'L', intermediate-CpG ==> 'I'
        sliding window = 500bp; offset = 5bp
        'H': GC_content > 0.55   +   Obs/Exp_CpG_ratio > _ratio_Cutoff_high (0.6)
        'I':
        'L': Obs/Exp_CpG_ratio < _ratio_Cutoff_low (default = 0.4)
    '''
    bed = Bed(bed)
    seq = bed.getSeq()
    window = 500

    # Obs/Exp_CpG_ratio list
    ratio=[]

    for i in range(0, len(seq) - window, 5):
        DNA = seq[i:i+500]
        r = Obs_Exp_CpG_ratio(DNA)
        gc = GC_content(DNA)
        if r > _ratio_Cutoff_high and gc > _gc_Cutoff_high:
            return 'H'
        ratio.append(r)
    #
    if not ratio or max(ratio) < _ratio_Cutoff_low:
        return 'L'
    # otherwise
    return 'I'
    
if __name__== "__main__":
    main()
##
##
"""
References:

     Mikkelsen, T.S. et al. Genome-wide maps of chromatin state in pluripotent and lineage-committed cells.
     Nature 448, 553-60 (2007)
"""

