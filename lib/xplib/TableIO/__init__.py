#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 12 Sep 2012 10:38:10
import BedIO
import GeneBedIO
import SimpleIO
import TransIO
import OddsRatioSNPIO
FormatToIterator = { "bed":BedIO.BedIterator,
                     "genebed":GeneBedIO.GeneBedIterator,
                     "simple":SimpleIO.SimpleIterator,
                     "transunit":TransIO.TransUnitIterator,
                     "oddsratiosnp":OddsRatioSNPIO.OddsRatioSNPIterator
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



