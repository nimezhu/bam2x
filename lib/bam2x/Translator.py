# Programmer : zhuxp
# Date: 
# Last-modified: 02-12-2014, 00:21:15 EST
import types
import pysam
from bam2x.Annotation import BED12 as Bed12
from bam2x import Tools

def BamToBed12(handle,**kwargs):
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
            if isinstance(kwargs["references"],str):
                chr=kwargs["references"]  
            else:
                chr=kwargs["references"][i.tid];
        else:
            try:
                 chr=handle.references[i.tid];
            except:
                 chr="chr"
        if kwargs.has_key("strand"):
            if kwargs["strand"]=="read1" or kwargs["strand"]=="firstMate":
                read1=True
            else:
                read1=False
        else:
            read1=True   
        start=i.pos
        end=i.aend
        name=i.qname
        cds_start=start
        cds_end=start
        itemRgb="0,0,0"
        '''
        debug
        import sys
        if i.cigar is None:
            print >>sys.stderr,"why cigar is Nonetype?"
            print >>sys.stderr,i
            exit(0)
        end of debug
        '''
        if i.cigar==None: continue # IGNORE THIS READS?
        (block_starts,block_sizes)=Tools.cigar_to_coordinates(i.cigar);
        if i.is_read1 and not read1:
            strand=Tools.reverse_strand(strand)
        elif i.is_read2 and read1:
            strand=Tools.reverse_strand(strand)
        bed=Bed12(chr,start,end,name,score,strand,cds_start,cds_end,itemRgb,len(block_sizes),block_sizes,block_starts)
        yield bed

