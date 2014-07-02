# Programmer : zhuxp
# Date: 
# Last-modified: 07-02-2014, 12:02:42 EDT
__all__=['BamI','TabixI','MetaDBI',"BinIndexI","BamlistI","TwoBitI","GenomeI"]
from DB import *
import os
from bam2x import IO,Annotation,TableIO

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

def smart_init(handle,**dict):
    '''
    test version
    '''
    if isinstance(handle,str):
        fn,ext=os.path.splitext(handle)
        if ext==".bam":
            dbi=FormatToDBI["bam"]
            return dbi(handle,**dict)
        elif ext==".gz":
            if os.path.isfile(handle+".tbi"):
                fn1,ext1=os.path.splitext(fn)
                if ext1==".bed":
                    col_num=IO.get_col_number(handle)
                    t=ext1+str(col_num)
                else:
                    t=ext1
                t=t[1:]

                if TableIO.hclass.has_key(t):
                    return TabixI(handle,cls=t,**dict)
                else:
                    return TabixI(handle,**dict)
            else:
                fn1,ext1=os.path.splitext(fn)
                if ext1==".bed":
                    col_num=IO.guess(handle)
                    t=ext1+col_num
                else:
                    t=ext1
                t=t[1:]
                if TableIO.hclass.has_key(t):
                    return BinIndexI(handle,cls=t,**dict)
                else:
                    return BinIndexI(handle,**dict)
        elif ext==".bed":
            col_num=IO.get_col_number(handle)
            t=ext+str(col_num)
            if TableIO.hclass.has_key(t):
                return BinIndexI(handle,cls=t,**dict)
            else:
                if col_num >=12 :
                    return BinIndexI(handle,cls="bed12",**dict)
                elif col_num>=6:
                    return BinIndexI(handle,cls="bed6",**dict)
                else:
                    return BinIndexI(handle,cls="bed3",**dict)

        else:
            if TableIO.hclass.has_key(t):
                return BinIndexI(handle,cls=t,**dict)
    elif isinstance(handle,file):
        return smart_init(handle.name,**dict)
    else:
        return BinIndex(handle,**dict)
 


        

