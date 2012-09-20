#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 19 Sep 2012 23:36:53

import os,sys,argparse
from xplib.Annotation import *
from xplib import TableIO
import pysam
'''
BASIC QUERY FUNCTIONS
'''
class DBI:
    def __init__(self,data,**dict):
        self.data=data
    def query(self,x):
        pass
class BinIndexI(DBI):
    def __init__(self,file,**dict):
        format=dict['format']
        self.data=Utils.readIntoBinIndex(file,format)
    def query(self,x):
        for i in Utils.iterOverlapFeature(x,self.data):
            yield i
        
class TabixI(DBI):
    def __init__(self, tabix_file_name,**dict):
        self.tabix_file_name=tabix_file_name
        try:
            self.data=pysam.Tabixfile(tabix_file_name)
        except:
            print >>sys.stderr,"WARNING: Can't init the tabix file",tabix_file_name

    def query(self,x):
        try:
            for item in self.data.fetch(x.chr,x.start,x.stop):
                yield item
        except:
           StopIteration

class BamI(DBI):
    hNtToNum={'a':0,'A':0,
          'c':1,'C':1,
          'g':2,'G':2,
          't':3,'T':3
         }
    Nt=['A','C','G','T']
    
    
    #def __init__(self,bamfile,flat=None):
    def __init__(self,bamfile,**dict):
        if type(bamfile)==type("str"):
            try:
                self.bamfile=pysam.Samfile(bamfile,"rb")
            except:
                print >>sys.stderr,"WARNING: Can't init the bam file",tabix_file_name
    def query(self,x):
        try:
            A=self.bamfile.pileup(x.chr,x.start,x.stop)
        except:
            print >>sys.stderr,"Can't pile up",x.chr,x.start,x.stop
            StopIteration
            return 
        s=[[0,0,0,0] for row in range(x.stop-x.start)]
        for pileupcolumn in A:
            j=pileupcolumn.pos-x.start
            if j<0: continue
            for pileupread in pileupcolumn.pileups:
                try:
                    nt=pileupread.alignment.seq[pileupread.qpos]
                    if BamI.hNtToNum.has_key(nt):
                        k=BamI.hNtToNum[nt]
                        s[j][k]+=1
                except:
                    pass
        for i in s:
            yield i


class BamsI(DBI):
    pass



