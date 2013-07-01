# Programmer : zhuxp
# Date:
# Last-modified: 06-19-2013, 16:40:52 EDT
from xplib.Annotation import Fimo
import types
import gzip
import SimpleIO 
def FimoIterator(handle,**dict):
    '''
    FIMO file iterator
    Wrapper in TableIO.parse(file or filename,"fimo")
    Usage:
        for i in TableIO.parse(file or filename,"fimo"):
            print i
        
    Bed is xplib.Annotation.Bed object.
    '''
    for i in SimpleIO.SimpleIterator(handle,**dict):
            yield Fimo(i)

