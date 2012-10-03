# Programmer : zhuxp
# Date: 
# Last-modified: 03 Oct 2012 16:49:40
import types
import gzip

def SimpleIterator(handle,**kwargs):
    '''
    Wrapper in TableIO.parse(filename or file,"simple") or TableIO.parse(filename or file). 
    The default iterator for TableIO.parse.
    Usage:
        for i in TableIO(f):
            print i

        i is a list contains str,int and float variables.
    the default separator is "\t".
    The sparator can be changed in this way.
    Example:
        for i in TableIO(f,sep=","):
            print i
    '''
    if type(handle)==type("s"):
        try:
            handle=handle.strip()
            x=handle.split(".")
            if x[-1]=="gz":
                handle=gzip.open(handle,"r")
            else:
                handle=open(handle,"r")
            
        except:
            raise ValueError("Can't open file %s"%handle)
    sep="\t"
    if kwargs.has_key("sep"):
        sep=kwargs["sep"]
    for line in handle:
        line=line.strip()
        if len(line)==0: continue
        if line[0]=="#": continue
        x=line.split(sep)
        for i,y in enumerate(x):
            try:
              x[i]=float(x[i])
              if x[i]==int(x[i]):
                  x[i]=int(x[i])
            except:
              pass
        yield x
    return
    
