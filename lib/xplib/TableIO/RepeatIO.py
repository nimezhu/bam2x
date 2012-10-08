# Programmer : zhuxp
# Date:
# Last-modified: 08 Oct 2012 17:51:53
from xplib.Annotation import Repeat
import types
import gzip
import SimpleIO 
def RepeatIterator(handle,**dict):
    '''
    Repeat Table file iterator
    Usage Example:
        from xplib.TableIO.BedIO import RepeatIterator
        for i in RepeatIterator(file or filename):
            print i
    Wrapper in TableIO.parse(file or filename,"repeat")
    Usage:
        for i in TableIO.parse(file or filename,"repeat"):
            print i
        
    Bed is xplib.Annotation.Bed object.
    '''
    for i in SimpleIO.SimpleIterator(handle,**dict):
        yield Repeat(i)

    
