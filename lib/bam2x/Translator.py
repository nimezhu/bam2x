# Programmer : zhuxp
# Date: 
# Last-modified: 07-31-2014, 14:42:36 EDT
import types
import pysam
from bam2x.Annotation import BED12 as Bed12
from bam2x.Annotation import Fragment
from bam2x import Tools

def BamToBed12(handle,uniq=False,**kwargs):
    '''
    handle is an bam iterator
    need references hash if handle is not filename.
    '''
    if type(handle)==type("string"):
        handle=pysam.Samfile(handle,"rb");
    use_chrs=False
    if kwargs.has_key("references"):
        if isinstance(kwargs["references"],str):
            chr=kwargs["references"]  
        else:
            use_chrs=True
            chrs=kwargs["references"];
    else:
        try:
            chrs=handle.references;
            use_chrs=True
        except:
            chr="chr"
 
    if kwargs.has_key("strand"):
        if kwargs["strand"]=="read1" or kwargs["strand"]=="firstMate":
            read1=True
        else:
            read1=False
    else:
        read1=True   
    for i in handle:
        if i.tid<0: continue
        if uniq:
            nh=_get_tag_score(i,"NH")
            if nh:
                if nh > 1: 
                    continue
            else:
                raise "no NH tag in your bam file"
        if use_chrs:
            chr=chrs[i.tid]
        strand="+"
        if i.is_reverse:
            strand="-"
        score=i.mapq
        
        '''
        test
        '''
        start=i.pos
        end=i.aend
        name=i.qname
        cds_start=start
        cds_end=start
        itemRgb="0,0,0"
        if not uniq:
            '''
            try to put NH score in itemRgb
            '''
            try:
                nh=_get_tag_score(i,"NH")
                if nh:
                    itemRgb=str(nh)+",0,0"
            except:
                pass
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
        #bed=Bed12(chr,i.pos,i.aend,i.qname,score,strand,i.pos,i.pos,itemRgb,len(block_sizes),block_sizes,block_starts)
        if nh:
            setattr(bed,"nh",nh)
        yield bed
def _get_tag_score(read,tag):
    for itag,iscore in read.tags:
        if tag==itag:
            return iscore
    return None


def BamToFragmentIterator(handle,**kwargs):
    '''
    handle is an bamfile or an iterator

    if it is an iterator , 
    a "bam" option could add to find the mate read which are not in the iterator.
    
    The strategy is:
        find all the pairs in iterator first and yield them 
        then
        try to find the rest reads' mate in bamfile.

    One problem is that :
        if too much read are clustered together
        it might be very cost memory if we read too much first read and the iterator still doesn't find their mates.
        TODO: FIX THIS PROBLEM!

    '''
    paired_reads={}
    if type(handle)==type("string"):
        handle=pysam.Samfile(handle,"rb");
    for read in handle:
        fragment_name=strip_mate_id(read.qname)
        if paired_reads.has_key(fragment_name):
            yield Fragment(paired_reads[fragment_name],read) #Paired End fragment
            del paired_reads[fragment_name]
        else:
            if read.is_qcfail or read.is_unmapped:
                continue 
            if  read.mate_is_unmapped or (not read.is_paired):
                yield Fragment(read) # Single end Fragment
            else:
                paired_reads[fragment_name]=read

    
    '''
    the rest of paired end read which haven't find mate yet 
    '''
    db=None
    if kwargs.has_key("bam"): 
        db=kwargs["bam"]
        if type(db)==type("string"):
            db=pysam.Samfile(db,"rb")
    elif isinstance(handle,pysam.Samfile): 
        db=handle
    if db is not None:
        pos=db.tell()
        for read in paired_reads.values():
            try: 
                mate=db.mate(read)
            except ValueError:
                mate=None
                continue
            finally: 
                db.seek(pos)
            yield Fragment(read,mate)
    else:
        for read in paired_reads.values():
            yield Fragment(read)
    del paired_reads


def strip_mate_id(read_name):
    '''
    this function was copied from http://pydoc.net/Python/misopy/0.4.7/misopy.sam_utils/ 
    
    Strip canonical mate IDs for paired end reads, e.g.
    #1, #2
    or:
    /1, /2
    '''
    if read_name.endswith("/1") or read_name.endswith("/2") or read_name.endswith("#1") or read_name.endswith("#2"):
        read_name = read_name[0:-3]
    return read_name


