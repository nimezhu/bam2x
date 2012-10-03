# Programmer : zhuxp
# Date: 
# Last-modified: 03 Oct 2012 16:49:10
import types
import pysam
from xplib.Annotation import Bed

def BamIterator(filename,**kwargs):
    '''
    iterator for reading a bam file. 
    Usage:
        from xplib.TableIO.BamIO import BamIterator
        for read in BamIterator(filename):
            print read
    read is a alignment in pysam.AlignedRead format.
    Wrapper In TableIO.parse(filename,"bam")
    Usage:
        for i in TableIO.parse(filename,"bam"):
            print i
    '''
    f=pysam.Samfile(filename,"rb")
    for i in f:
        yield i
def SamIterator(filename,**kwargs):
    '''
    iterator for reading a sam file. 
    Usage:
        from xplib.TableIO.BamIO import SamIterator
        for read in SamIterator(filename):
            print read
    read is a alignment in pysam.AlignedRead format.
    Wrapper In TableIO.parse(filename,"sam")
    Usage:
        for i in TableIO.parse(filename,"sam"):
            print i
    '''
    f=pysam.Samfile(filename,"r")
    for i in f:
        yield i
def BamToBedIterator(filename,**kwargs):
    '''
    iterator for reading a bam file, yield Bed Object instead of pysam.AlignedRead Object 
    Usage:
        from xplib.TableIO.BamIO import BamToBedIterator
        for read in BamToBedIterator(filename):
            print read
    read is a alignment in pysam.AlignedRead format.
    Wrapper In TableIO.parse(filename,"bam2bed")
    Usage:
        for i in TableIO.parse(filename,"bam2bed"):
            print i

    A simple bam2bed.py which will read bam file and print aligned read in bed format:
        import sys
        from xplib import TableIO
        filename=sys.args[1]
        for i in TableIO.parse(filename,"bam2bed"):
            print i
    '''

    f=pysam.Samfile(filename,"rb")
    for i in f:
        if i.tid<0:continue
        strand="+"
        if i.is_reverse:
            strand="-"
        score=i.mapq
        bed=Bed([f.references[i.tid],i.pos,i.aend,i.qname,score,strand])
        yield bed
