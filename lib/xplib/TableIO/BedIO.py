# Programmer : zhuxp
# Date:
# Last-modified: 06-12-2013, 13:13:51 EDT
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

    
