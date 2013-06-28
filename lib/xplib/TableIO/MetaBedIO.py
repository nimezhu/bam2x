# Programmer : zhuxp
# Date:
# Last-modified: 06-28-2013, 11:16:08 EDT
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
    
    header could be a list or filename 

    MetaBed is xplib.Annotation.MetaBed object.
    '''
    f=SimpleIO.SimpleIterator(handle,**dict)
    
    if dict.has_key("header") and dict["header"]==True:
        attr=f.next()
    elif dict.has_key("header") and isinstance(dict["header"],list):
        attr=dict["header"]
    elif dict.has_key("header") and isinstance(dict["header"],str):
        fh=TableIO.parse(dict["header"])
        attr=fh.next()
    if dict.has_key("header"):
        for i in f:
            yield MetaBed(value=i,attr=attr);
    else:
        for i in f:
            yield i

    
