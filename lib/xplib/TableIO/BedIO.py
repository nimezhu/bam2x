# Programmer : zhuxp
# Date:
# Last-modified: 01-22-2014, 12:17:23 EST
from xplib.Annotation import Bed,Bed12
import types
import gzip
import SimpleIO 
def BedIterator(handle,**dict):
    '''
    Bed file iterator
    Usage Example:
        from xplib.TableIO.BedIO import BedIterator
        for bed in BedIterator(file or filename):
            print bed
    Wrapper in TableIO.parse(file or filename,"bed")
    Usage:
        for i in TableIO.parse(file or filename,"bed"):
            print bed
        
    Bed is xplib.Annotation.Bed object.
    '''
    for i in SimpleIO.SimpleIterator(handle,**dict):
        if(len(i)>=12):
            yield Bed12(i)
        else:
            yield Bed(i)

    
from xplib.Tuple.Bed12Tuple import *
def BedTupleIterator(handle,**dict):
    '''
    Bed file Tuple iterator
    Usage Example:
        from xplib.TableIO.BedIO import BedTupleIterator
        for bed in BedIterator(file or filename):
            print bed
    Usage:
        for i in TableIO.parse(file or filename,"bed2tuple"):
            print bed
        
    '''
    for x in SimpleIO.SimpleIterator(handle,**dict):
        l=[]
        if len(x)==12:
            x[BLOCKSIZES]=x[BLOCKSIZES].strip(",").split(",")
            x[BLOCKSTARTS]=x[BLOCKSTARTS].strip(",").split(",")
            for i in range(x[BLOCKCOUNT]):
                x[BLOCKSIZES][i]=int(x[BLOCKSIZES][i])
                x[BLOCKSTARTS][i]=int(x[BLOCKSTARTS][i])
            x[BLOCKSIZES]=tuple(x[BLOCKSIZES])
            x[BLOCKSTARTS]=tuple(x[BLOCKSTARTS])
        yield tuple(x)
    

    
