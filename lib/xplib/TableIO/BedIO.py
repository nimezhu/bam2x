# Programmer : zhuxp
# Date:
# Last-modified: 03 Oct 2012 16:49:20
from xplib.Annotation import Bed
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
        yield Bed(i)

    
