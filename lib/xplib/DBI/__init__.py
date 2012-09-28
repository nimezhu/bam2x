#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 28 Sep 2012 10:07:22
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
             "vcf":BinIndexI,
            }
def query(x,dbi):
    for i in dbi.query(x):
        yield i
        
def init(handle,dbformat,**dict):
        if dbformat in FormatToDBI:
            dbi=FormatToDBI[dbformat]
            return dbi(handle,format=dbformat,**dict)

