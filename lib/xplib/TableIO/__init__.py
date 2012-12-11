# Programmer : zhuxp
# Date: 
# Last-modified: 12-11-2012, 11:50:26 CST
import BedIO
import GeneBedIO
import SimpleIO
import TransIO
import OddsRatioSNPIO
import BamIO
import VCFIO
import RepeatIO
import MetaBedIO
FormatToIterator = { "bed":BedIO.BedIterator,
                     "genebed":GeneBedIO.GeneBedIterator,
                     "simple":SimpleIO.SimpleIterator,
                     "transunit":TransIO.TransUnitIterator,
                     "oddsratiosnp":OddsRatioSNPIO.OddsRatioSNPIterator,
                     "aps":OddsRatioSNPIO.OddsRatioSNPIterator,
                     "bam":BamIO.BamIterator,
                     "sam":BamIO.SamIterator,
                     "bam2bed":BamIO.BamToBedIterator,
                     "vcf":VCFIO.VCFIterator,
                     "repeat":RepeatIO.RepeatIterator,
                     "metabed":MetaBedIO.MetaBedIterator
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

