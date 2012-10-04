# Programmer : zhuxp
# Date: 
# Last-modified: 03 Oct 2012 21:39:31
__all__=['BamI','TabixI','MetaDBI',"BinIndexI","BamlistI"]
from DB import *

FormatToDBI = { 
             "bed":BinIndexI,
             "genebed":BinIndexI,
             "bam":BamI,
             "bams":BamlistI,
             "bamlist":BamlistI,
             "tabix":TabixI,
             "oddsratiosnp":BinIndexI,
             "vcf":BinIndexI,
             "2bit":TwoBitI,
             "genome":TwoBitI,
            }
def query(x,dbi):
    '''
    DBI.query is wrapper for query all kinds of data|file|database. 
    
    Usage Example 1:
        
        dbi=DBI.init(file or filname,"bed")
        bed=Bed(["chr1",1,1000,".",".","."]) 
        for i in DBI.query(bed,dbi):
            print i
        
        the program above will return all overlap features.
    Usage Example 2:
        
        dbi=DBI.init(filname,"bam")
        bed=Bed(["chr1",1,10,".",".","."]) 
        for i in DBI.query(bed,dbi):
            print i
        
        the yeild variable i is a list with length 4.
        [A number,C number, G number, T number] in each postion of bed region.
    
    The Query Result is in different format for different db format.
    Most flat annotation files are stored into BinIndex data structure for query.

    This Wrapper is not neccesary.
    the code below is the same.
    
    Usage Example 1:
        dbi=DBI.init(file or filname,"bed")
        bed=Bed(["chr1",1,1000,".",".","."]) 
        for i in dbi.query(bed):
            print i
    '''
    for i in dbi.query(x):
        yield i
        
def init(handle,dbformat,**dict):
    '''
    A Wrapper for initialize the DBI
    Usage:
        dbi=DBI.init(data or filename or list,dbformat)

    Examples:
        dbi=DBI.init(file,"bed")

        dbi=DBI.init(list,"bed")
        
        dbi=DBI.init("filename.bam","bam")
        
        dbi=DBI.init(["file1.bam","file2.bam","bamlist")


    '''
    if dbformat in FormatToDBI:
        dbi=FormatToDBI[dbformat]
        return dbi(handle,format=dbformat,**dict)

