#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 03 Oct 2012 13:18:15
import BedIO
import GeneBedIO
import SimpleIO
import TransIO
import OddsRatioSNPIO
import BamIO
import VCFIO
FormatToIterator = { "bed":BedIO.BedIterator,
                     "genebed":GeneBedIO.GeneBedIterator,
                     "simple":SimpleIO.SimpleIterator,
                     "transunit":TransIO.TransUnitIterator,
                     "oddsratiosnp":OddsRatioSNPIO.OddsRatioSNPIterator,
                     "aps":OddsRatioSNPIO.OddsRatioSNPIterator,
                     "bam":BamIO.BamIterator,
                     "sam":BamIO.SamIterator,
                     "bam2bed":BamIO.BamToBedIterator,
                     "vcf":VCFIO.VCFIterator
                   }
def parse(handle,format="simple",**dict):
    """
    - handle  - handle to the file, or the filename
    - format  - lower case string describing the file format
                example 'bed' 'genebed' 'bam' 'sam' 'vcf'
    Example:
        from xplib import TableIO
        for i in TableIO.parse(file or filename,"bed"):
            print i

    """
    mode='rU'
    if format in FormatToIterator:
        iterator_generator=FormatToIterator[format]
        i=iterator_generator(handle,**dict)
    for r in i:
        yield r

