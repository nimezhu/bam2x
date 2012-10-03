# Programmer : zhuxp
# Date: 
# Last-modified: 03 Oct 2012 16:49:52
from xplib.Annotation import VCF
import types
import SimpleIO
def VCFIterator(handle,**kwargs):
    '''
    Wrapper in TableIO.parse(file or filename,"vcf")
    Usage:
        for i in TableIO.parse(file or filename,"vcf"):
            print i
        i is a Annotation.VCF class object.
    '''
    for x in SimpleIO.SimpleIterator(handle,**kwargs):
        b=VCF(x)
        yield b
    
