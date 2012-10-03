# Programmer : zhuxp
# Date: 
# Last-modified: 03 Oct 2012 16:50:07

import os,sys,argparse
import pysam
import random
from math import log,sqrt
from xplib import TableIO
from xplib.Annotation import *
from scipy.stats import chi2


hNtToNum={'a':0,'A':0,
          'c':1,'C':1,
          'g':2,'G':2,
          't':3,'T':3
         }
Nt=['A','C','G','T']

        
def ntDisToChi2(disA,disB):
    '''
    return :
       APS 
    '''
    x=OddsRatioSNP(A=disA,B=disB)
    return x.APS
