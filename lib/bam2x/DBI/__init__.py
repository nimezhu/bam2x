# Programmer : zhuxp
# Date: 
# Last-modified: 02-12-2014, 12:43:09 EST
__all__=['BamI','TabixI','MetaDBI',"BinIndexI","BamlistI","TwoBitI","GenomeI"]
from DB import *

FormatToDBI = { 
             "flatfile":BinIndexI,
             "binindex":BinIndexI,
             "bam":BamI,
             "bams":BamlistI,
             "bamlist":BamlistI,
             "tabix":TabixI,
             "genome":GenomeI,
             "bigwig":BigWigI
            }
def query(x,dbi,**dict):
    '''
    DBI.query is wrapper for query all kinds of data|file|database. 
    '''
    for i in dbi.query(x,**dict):
        yield i
        
def init(handle,dbformat,**dict):
    '''
    A Wrapper for initialize the DBI
    '''
    if dbformat in FormatToDBI:
        dbi=FormatToDBI[dbformat]
        return dbi(handle,**dict)
    else:
        dbi=BinIndexI
        return dbi(handle,**dict)

