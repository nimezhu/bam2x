# Programmer : zhuxp
# Date:
# Last-modified: 12-11-2012, 14:02:03 CST
from xplib.Annotation import MetaBed
import types
import gzip
import SimpleIO 
def MetaBedIterator(handle,**dict):
    '''
    MetaBed file iterator
    Usage Example:
        from xplib.TableIO.MetaBedIO import BedIterator
        for bed in BedIterator(file or filename):
            print bed
    Wrapper in TableIO.parse(file or filename,"metabed",attr=header)
    Usage:
        for i in TableIO.parse(file or filename,"metabed",header=True):
            print i
        
    MetaBed is xplib.Annotation.MetaBed object.
    '''
    f=SimpleIO.SimpleIterator(handle,**dict)
    if dict.has_key("header") and dict["header"]==True:
        attr=f.next()
        dict["attr"]=attr
    for i in f:
        yield MetaBed(value=i,attr=dict["attr"]);

    
