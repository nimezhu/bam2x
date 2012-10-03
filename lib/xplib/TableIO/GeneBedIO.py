# Programmer : zhuxp
# Date: 
# Last-modified: 03 Oct 2012 16:49:30
from xplib.Annotation import GeneBed
import types
import SimpleIO
def GeneBedIterator(handle,**kwargs):
    '''
    '''
    for x in SimpleIO.SimpleIterator(handle,**kwargs):
        b=GeneBed(x)
        yield b
    
