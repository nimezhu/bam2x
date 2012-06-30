#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 30 Jun 2012 00:49:29
import BedIO
import GeneBedIO
import SimpleIO
FormatToIterator = { "bed":BedIO.BedIterator,
                     "genebed":GeneBedIO.GeneBedIterator,
                     "simple":SimpleIO.SimpleIterator
                     
                   }
FormatToWrite    = {  
                     "bed":BedIO.BedWriter,
                     "genebed":GeneBedIO.GeneBedWriter,
                     "simple":SimpleIO.SimpleWriter
                   }
def parse(handle,format):
    """
    - handle  - handle to the file, or the filename
    - format  - lower case string describing the file format
                example 'bed' 'genebed' 'gtf'
    """
    mode='rU'
    if format in FormatToIterator:
        iterator_generator=FormatToIterator[format]
        i=iterator_generator(handle)
    for r in i:
        yield r

def write(handle,format):
    pass



