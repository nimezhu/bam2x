# Programmer : zhuxp
# Date: 
# Last-modified: 08-28-2013, 14:37:27 EDT
import types
import pysam
from xplib.Annotation import Bed,Bed12
from xplib import Tools

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
def BamToBed12Iterator(handle,**kwargs):
    '''
    handle is an bam iterator
    need references hash if handle is not filename.
    '''
    if type(handle)==type("string"):
        handle=pysam.Samfile(handle,"rb");
    for i in handle:
        if i.tid<0: continue
        strand="+"
        if i.is_reverse:
            strand="-"
        score=i.mapq
        
        '''
        test
        '''
        if kwargs.has_key("references"):
            chr=kwargs["references"][i.tid];
        else:
            try:
                 chr=handle.references[i.tid];
            except:
                 chr="chr"
        
        start=i.pos
        end=i.aend
        name=i.qname
        cds_start=start
        cds_end=start
        itemRgb="0,0,0"
        (block_starts,block_sizes)=Tools.cigar_to_coordinates(i.cigar);
        bed=Bed12([chr,start,end,name,score,strand,cds_start,cds_end,itemRgb,len(block_sizes),block_sizes,block_starts])
        yield bed
