#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 19 Sep 2012 23:34:21
__all__=['BamI','TabixI','DBI',"BinIndexI"]
from DBI import *

FormatToDBI = { 
             "binindex":BinIndexI,
             "bam":BamI,
             "tabix":TabixI

            }
#FormatToWrite    = {  
#                     "bed":BedIO.BedWriter,
#                     "genebed":GeneBedIO.GeneBedWriter,
#                     "simple":SimpleIO.SimpleWriter
#                   }
def query(x,dbi):
    print "IN Query"
    for i in dbi.query(x):
        yield i
        
def init(handle,dbformat="binindex",**dict):
        if dbformat in FormatToDBI:
            dbi=FormatToDBI[dbformat]
            return dbi(handle,**dict)

