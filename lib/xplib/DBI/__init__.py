#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 10-02-2012, 02:12:22 CDT
__all__=['BamI','TabixI','DBI',"BinIndexI","BamlistI"]
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

