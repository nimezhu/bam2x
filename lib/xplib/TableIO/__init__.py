#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 24 Jul 2012 14:56:22
import BedIO
import GeneBedIO
import SimpleIO
import TransIO
FormatToIterator = { "bed":BedIO.BedIterator,
                     "genebed":GeneBedIO.GeneBedIterator,
                     "simple":SimpleIO.SimpleIterator,
                     "transunit":TransIO.TransUnitIterator
                     
                   }
#FormatToWrite    = {  
#                     "bed":BedIO.BedWriter,
#                     "genebed":GeneBedIO.GeneBedWriter,
#                     "simple":SimpleIO.SimpleWriter
#                   }
def parse(handle,format):
    """
    - handle  - handle to the file, or the filename
    - format  - lower case string describing the file format
                example 'bed' 'genebed' 'gtf'
    Example:
        from xplib import TableIO
        for i in TableIO.parse(file or filename,"bed"):
            print i

    """
    mode='rU'
    if format in FormatToIterator:
        iterator_generator=FormatToIterator[format]
        i=iterator_generator(handle)
    for r in i:
        yield r

def write(handle,format):
    pass



