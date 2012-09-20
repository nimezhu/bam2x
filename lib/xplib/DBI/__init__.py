#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 20 Sep 2012 16:33:43
__all__=['BamI','TabixI','DBI',"BinIndexI"]
from DBI import *

FormatToDBI = { 
             "bed":BinIndexI,
             "genebed":BinIndexI,
             "bam":BamI,
             "bams":BamlistI,
             "bamlist":BamlistI,
             "tabix":TabixI,
             "oddsratiosnp":BinIndexI,
            }
def query(x,dbi):
    for i in dbi.query(x):
        yield i
        
def init(handle,dbformat,**dict):
        if dbformat in FormatToDBI:
            dbi=FormatToDBI[dbformat]
            return dbi(handle,format=dbformat,**dict)

