# Programmer : zhuxp
# Date:  Sep 2012
# Last-modified: 09-24-2014, 17:49:44 EDT
from string import upper,lower
from bam2x.Annotation import BED6 as Bed
from bam2x.Annotation import BED12 as Bed12
from bam2x.Annotation import BED12 
from bam2x.Annotation import BED3
from bam2x.Annotation import Fragment
import bam2x
import copy
import logging
import itertools
# __all__=["IO","codon"]

hNtToNum={'a':0,'A':0,
          'c':1,'C':1,
          'g':2,'G':2,
          't':3,'T':3
         }
Nt=['A','C','G','T']

suffixToFormat={
    'fa':'fasta',
    'fq':'fastq',
    'genetab':'genebed',
    'bw':'bigwig',
    'tab':'genebed',
    '2bit':'genome'

}

def get_flank_region(bed,up,down,chr_lengths=None):
    down_bp=down
    up_bp=up
    if bed.strand=="-":
        start=bed.start-down
        stop=bed.stop+up
        if start < 0:
            start=0
            down_bp=bed.start
        if chr_lengths is not None and chr_lengths.has_key(bed.chr):
            if stop > chr_lengths[bed.chr]:
                stop=chr_lengths[bed.chr]
                up_bp=stop-bed.stop
    elif bed.strand=="+" or bed.strand==".":
        start=bed.start-up
        stop=bed.stop+down
        if start < 0:
            start=0
            up_bp=bed.start
        if chr_lengths is not None and chr_lengths.has_key(bed.chr):
            if stop > chr_lengths[bed.chr]:
                stop=chr_lengths[bed.chr]
                down_bp=stop-bed.stop
    return bed._replace(id="{id}_up{up}_down{down}".format(id=bed.id,up=up_bp,down=down_bp),start=start,stop=stop)
    
def rc(seq):
   comps = {'A':"T", 'C':"G", 'G':"C", 'T':"A",
           'B':"V", 'D':"H", 'H':"D", 'K':"M",
           'M':"K", 'R':"Y", 'V':"B", 'Y':"R",
           'W':'W', 'N':'N', 'S':'S'}
   return ''.join([comps[x] for x in seq.upper()[::-1]])
def shuffle(seq):
   import random
   a=list(seq)
   random.shuffle(a)
   return "".join(a)

def seq_wrapper(seq,width=60):
    s=""
    seqlen=len(seq)
    for i in range(0,seqlen,width):
        stop=i+width
        if stop>seqlen:stop=seqlen
        s+=seq[i:stop]+"\n"
    return s
def distance(A,B):
    if A.chr!=B.chr: return None
    if overlap(A,B): return 0
    return min(abs(A.start-B.stop),abs(A.stop-B.start))
def translate_coordinate(coord,bed,reverse=False):
    '''
    translate bed's coordiante based on coord
    if reverse is True:
        bed is in coord's coordinates and translate it back to the chromosome coordinates
    coord is a simple BED6 or BED3, didn't consider the splicing in this function
    '''
    if reverse:
        #TO TEST
        chr=coord.chr
        id=bed.id
        if coord.strand=="+" or coord.strand==".":
            return (coord.start+bed.start,coord.start+bed.stop,bed.strand)
        else:
            strand="."
            if bed.strand=="-":strand="+"
            if bed.strand=="+":strand="-"
            return (coord.stop-bed.stop,coord.stop-bed.start,strand)

    else:
        if coord.chr!=bed.chr: return None
        if coord.strand=="+" or coord.strand==".":
            return (bed.start-coord.start,bed.stop-coord.start,bed.strand)
        if coord.strand=="-":
            strand="."
            if bed.strand=="-":strand="+"
            if bed.strand=="+":strand="-"
            return (coord.stop-bed.stop,coord.stop-bed.start,strand)
def translate_coordinates(coord,bed,reverse=False): # bed is Bed12 format
    '''
    Translate Bed12 Object bed to BED3 or BED6 coordinates ( no splicing )
    if reverse is True
        bed is in coord's coordinates and translate it back to the chromosome coordinates.
        TO TEST
    '''
    if reverse:
        chr=coord.chr
    else:
        chr=coord.id
    id=bed.id
    (start,stop,strand)=translate_coordinate(coord,bed,reverse)
    score=bed.score
    if isinstance(bed,Bed12):
        (cds_start,cds_stop,cds_strand)=translate_coordinate(coord,Bed(bed.chr,bed.cds_start,bed.cds_stop,bed.id+"_cds",bed.score,bed.strand),reverse)
        itemRgb=bed.itemRgb
        blockCount=bed.blockCount
        blockSizes=list(bed.blockSizes)
        blockStarts=list(bed.blockStarts)
        
        if coord.strand=="+" or coord.strand==".":
            for i,x in enumerate(blockStarts):
                blockStarts[i]=blockStarts[i]
        elif coord.strand=="-":
            for i,x in enumerate(blockStarts):
                #print "debug",i,x,blockStarts,blockSizes
                blockStarts[i]+=blockSizes[i]
            blockStarts=blockStarts[::-1]
            for i,x in enumerate(blockStarts):
                #print blockStarts[i]
                blockStarts[i]=bed.stop-(blockStarts[i]+bed.start)
            blockSizes=blockSizes[::-1]
        return Bed12(chr,start,stop,id,score,strand,cds_start,cds_stop,itemRgb,blockCount,tuple(blockSizes),tuple(blockStarts))
    else:
        C=copy.copy(bed)
        C.chr=chr
        C.start=start
        C.stop=stop
        C.end=stop
        C.strand=strand
        if bed.__dict__.has_key('pos'):
            C.pos=C.start+1  # position is 1  index
        return C
        
def merge_bed(bedA,bedB,id="noname"):
    '''
    merge two bed into a bed , ratain all exon region.
    '''
    if bedA.chr!=bedB.chr: return None
    l=[]
    for i in bedA.Exons():
        l.append(i)
    for i in bedB.Exons():
        l.append(i)
    return _merge_bed6(l,id=id)
    #TODO 

def merge_beds(beds,id="noname"):
    '''
    merge a bed list
    '''
    l=[]
    chr=beds[0].chr
    for bed in beds:
        if bed.chr!=chr: return None
        for e in bed.Exons():
            l.append(e)
    return _merge_bed6(l,id=id)


def _merge_bed6(beds,id="noname"):
    '''
    a simple turing state change
    '''
    l=[]
    strand=beds[0].strand
    for i in beds:
        l.append((i.start,-1))
        l.append((i.stop,1))
        if strand!=i.strand: strand="."
    chr=beds[0].chr
    l.sort()
    switch=0
    
    start=l[0][0]
    end=l[-1][0]
    
    state=0
    assert i[0][1] > 1
    last_pos=l[0][0]
    state=1
    switch=1
    blockStarts=[]
    blockSizes=[]
    for i in l[1:]:
        state-=i[1]
        if switch==1 and state==0:
            blockStarts.append(last_pos-start)
            blockSizes.append(i[0]-last_pos)
            switch=0
        if switch==0 and state>0:
            last_pos=i[0]
            switch=1
    
    assert state==0
    blockCount=len(blockSizes)
    return BED12(chr,start,end,id,0.0,strand,start,start,"0,0,0",blockCount,blockSizes,blockStarts) 

def translate(bedA,bedB):
    #TODO
    '''
    INPUT: two BED12 bedA and bedB 
    return the location for bedA in merge bed, the location for bedB in merge bed, merged bed
    '''
    meta=merge_bed(bedA,bedB,bedA.id+"_"+bedB.id+"_merged")
    new_bedA=_translate_to_meta(meta,bedA)
    new_bedB=_translate_to_meta(meta,bedB)
    return new_bedA,new_bedB,meta
def _translate_to_meta(meta,bed):
    '''
    '''
    l=[]
    '''
    for i in meta.Exons():
        l.append((i.start,-1,1))
        l.append((i.stop,1,1))
    '''
    
    for start,size in itertools.izip(meta.blockStarts,meta.blockSizes):
        l.append((meta.start+start,-1,1))
        l.append((meta.start+start+size,1,1))
    '''
    for i in bed.Exons():
        l.append((i.start,-1,2))
        l.append((i.stop,1,2))
    '''
    if hasattr(bed,"blockStarts"):
        for start,size in itertools.izip(bed.blockStarts,bed.blockSizes):
            l.append((bed.start+start,-1,2))
            l.append((bed.start+start+size,1,2))
    else:
        l.append((bed.start,-1,2))
        l.append((bed.stop,1,2))
    if hasattr(bed,"cds_start"):
        l.append((bed.cds_start,0,3))
        l.append((bed.cds_stop,0,4))
    else:
        '''
        no cds 
        '''
        l.append((bed.start,0,3))
        l.append((bed.start,0,4))

    l.sort()
    meta_state=0
    bed_state=0
    meta_last_pos=0
    meta_coordinate=0
    tl=[]
    for i in l:
        if i[2]==1: #META
            if i[1]==-1 and meta_state==0:
                meta_state=1
                meta_last_pos=i[0]
            elif i[1]==1 and meta_state==1:
                meta_state=0
                meta_coordinate+=i[0]-meta_last_pos 
                meta_last_pos=i[0]
        elif i[2]==2:  #BED
            tl.append((i[0]-meta_last_pos+meta_coordinate,i[1]))
        elif i[2]==3:
            cds_start=i[0]-meta_last_pos+meta_coordinate
        elif i[2]==4:
            cds_stop=i[0]-meta_last_pos+meta_coordinate
    tl.sort()
    state=0
    switch=0
    blockStarts=[]
    blockSizes=[]
    for i in tl:
        state-=i[1]
        if switch==0:
            if state > 0:
                b_start=i[0]
                switch=1
                blockStarts.append(b_start)
        elif switch==1:
            if state == 0:
                b_end=i[0]
                switch=0
                blockSizes.append(b_end-b_start)
    blockCount=len(blockSizes)
    itemRgb="0,0,0"
    try:
        itemRgb=bed.itemRgb
    except:
        pass
    
    if meta.strand=="." or meta.strand=="+":
        return BED12(meta.id,blockStarts[0],blockStarts[-1]+blockSizes[-1],bed.id,bed.score,bed.strand,cds_start,cds_stop,itemRgb,blockCount,tuple(blockSizes),tuple([blockStart-blockStarts[0] for blockStart in blockStarts]))
    else:
        len_meta=meta.cdna_length()
        strand=reverse_strand(bed.strand)
        return BED12(meta.id,len_meta-blockStarts[-1]-blockSizes[-1],len_meta-blockStarts[0],bed.id,bed.score,strand,len_meta-cds_stop,len_meta-cds_start,itemRgb,blockCount,tuple(blockSizes[::-1]),tuple([len_meta-i0-j0-(len_meta-blockStarts[-1]-blockSizes[-1]) for i0,j0 in itertools.izip(blockStarts[::-1],blockSizes[::-1])]))
    
    
_translate=_translate_to_meta

def reverse_translate(meta,bed):
    '''
    reverse translate bed in meta coordinates ( a BED12 object ) to chromosome coordinates
    '''
    assert meta.id==bed.chr
    start=bed.start
    stop=bed.stop
    start_sign=True
    stop_sign=True

    if meta.strand=="-":
        for blockSize,blockStart in itertools.izip(meta.blockSizes[::-1],[ i0+j0 for i0,j0 in itertools.izip(meta.blockStarts[::-1],meta.blockSizes[::-1])]):
            start-=blockSize
            if start <= 0 and start_sign:
                new_stop=blockStart-blockSize-start+meta.start
                start_sign=False
            stop-=blockSize
            if stop <= 0 and stop_sign:
                new_start=blockStart-blockSize-stop+meta.start
                stop_sign=False
                break
    else:
        for blockSize,blockStart in itertools.izip(meta.blockSizes,meta.blockStarts):
            start-=blockSize
            stop-=blockSize
            if start <= 0 and start_sign:
                new_start=blockStart + start + blockSize + meta.start
                start_sign=False
            if stop  <= 0 and stop_sign:
                new_stop=blockStart + stop + blockSize + meta.start
                stop_sign=False
                break
    return meta._slice(new_start,new_stop,bed.id)






def overlap(A,B):
    '''
    if A is overlapping with B.
    A and B are ? extends Bed class.
    '''
    if(A.chr != B.chr) : return 0
    if (A.stop < B.start) : return 0
    if (B.stop < A.start) : return 0
    return 1
def find_nearest(bed,dbi,extends=50000,**dict):
    start=bed.start-extends
    stop=bed.stop+extends
    chr=bed.chr
    if start<0: start=0
    new_bed=BED3(chr,start,stop)

    results=dbi.query(new_bed,**dict)
    d=2*extends
    flag=0
    
    for result in results:
        d0=distance(bed,result)
        if d0<d:
            d=d0
            nearest=result
            if  result.strand=="." or bed.strand==".":
                strand="."
            elif result.strand==bed.strand:
                strand="+"
            else:
                strand="-"
            flag=1
    if flag==0:
        return (None,None,None)
    else:
        return (d,nearest,strand)

def extend_slice(gene,start,end):
    if gene.start > start:
        gene=extend_start(gene,start)
    if gene.end < end:
        gene=extend_end(gene,end)
    gene=gene._slice(start,end)
    return gene

def extend_start(bed12,start):
    '''
        extend the start for easy gene
    '''
    if bed12.strand=="+" : exon=bed12.Exons()[0]
    elif bed12.strand=="-" : exon=bed12.Exons()[bed12.blockCount-1]
    else: exon=bed12.Exons()[0]
    if start < exon.stop: 
        blockStarts=[]
        offset=start-bed12.start
        for i in bed12.blockStarts:
            blockStarts.append(bed12.start+i-start)
        blockStarts[0]=0
        blockSizes=[ i for i in bed12.blockSizes]
        blockSizes[0]=bed12.blockSizes[0] - offset
        return bed12._replace(blockStarts=blockStarts,blockSizes=blockSizes,start=start)
    else:
        logging.WARN("start site might be in intron , don't extend the start",bed12)
        return bed12
def extend_end(bed12,stop):
    '''
        extend the end for easy gene
    '''
    if bed12.strand=="-" : exon=bed12.Exons()[0]
    elif bed12.strand=="+" : exon=bed12.Exons()[bed12.blockCount-1]
    else: exon=bed12.Exons()[0]
    if stop > exon.start: 
        blockSizes=[ i for i in bed12.blockSizes]
        blockSizes[bed12.blockCount-1]=bed12.blockSizes[bed12.blockCount-1] - ( bed12.stop - stop)
        return bed12._replace(blockSizes=blockSizes,stop=stop)
    else:
        logging.WARN("stop site might be in intron, don't fix the stop",bed12)
        return bed12  




def compatible(a,b,**kwargs):
    '''
    VERSION: TEST
    if a and b are compatible return true
    a and b are BED12 class 
    definition of compatible
       the overlap region should be same transcript structure.
    '''
    if (not overlap(a,b)): return True;
    
    if (a.strand!=b.strand): 
        if kwargs.has_key("unstranded") and kwargs["unstranded"]:
            pass
        else:    
            return False;
    '''
    if two bed is not overlap, they are compatible.
    '''
    start=a.start;
    if (start < b.start): start=b.start 
    '''
    start is the max start of a.start and b.start 
    '''
    stop=a.stop
    if (stop > b.stop): stop=b.stop
    '''
    stop is the min stop of a.stop and b.stop
    find the overlap region from [start,stop)
    '''
    '''
    sliced_a=a._slice(start,stop)
    sliced_b=b._slice(start,stop)
    if sliced_a is None: return False
    if sliced_b is None: return False
    if sliced_a.blockCount!=sliced_b.blockCount : 
        return False
    else:
        for i in range(sliced_a.blockCount):
            if sliced_a.blockStarts[i]!=sliced_b.blockStarts[i]:
                return False
            if sliced_a.blockSizes[i]!=sliced_b.blockSizes[i]:
                return False
    return True
    ''' 
    l=[] 
    for bstart,bsize in itertools.izip(a.blockStarts,a.blockSizes):
        l.append((a.start+bstart,0))
        l.append((a.start+bstart+bsize,1))
    for bstart,bsize in itertools.izip(b.blockStarts,b.blockSizes):
        l.append((b.start+bstart,2))
        l.append((b.start+bstart+bsize,3))
    l.append((a.start,5))
    l.append((a.stop,6))
    l.append((b.start,7))
    l.append((b.stop,8))
    l.sort()
    last_pos=0
    a_on=False
    b_on=False
    a_block_on=False
    b_block_on=False
    for i in l:
        if i[0]!=last_pos:
            if(a_on and b_on):
                if(a_block_on != b_block_on):
                    return False
            last_pos=i[0]
        if i[1]==0:
            a_block_on=True
        elif i[1]==1:
            a_block_on=False
        elif i[1]==2:
            b_block_on=True
        elif i[1]==3:
            b_block_on=False
        elif i[1]==5:
            a_on=True
        elif i[1]==6:
            a_on=False
            break
        elif i[1]==7:
            b_on=True
        elif i[1]==8:
            b_on=False
            break
    return True


    


    



def reverse_strand(i):
    if i=="+": return "-"
    if i=="-": return "+"
    if i==".": return "."
    if type(i)==type(1): return -i
def compatible_with_transcript(read,transcript,**kwargs): 
    '''
    if reads is compatible with the transcript( full length or longer than reads)
    
    !!! if read is fragment, then need add "references"= bamfile.references
    '''
    # read must in the transcript
    if isinstance(read,Fragment):
        for i0,i in enumerate(bam2x.TableIO.parse(read.reads,"bam2bed12",**kwargs)):
            if i.start < transcript.start or i.stop > transcript.stop :
                return False
            else:
                if kwargs.has_key("strand"): 
                    if kwargs["strand"]=="read1":
                        if read.reads[i0].is_read2:
                            i.strand=reverse_strand(i.strand)
                    elif kwargs["strand"]=="read2":
                        if read.reads[i0].is_read1:
                            i.strand=reverse_strand(i.strand)
                if not compatible(i,transcript,**kwargs):
                    return False
        return True
    else:
        if kwargs.has_key("strand") and kwargs["strand"]:
            if read.strand != transcript.strand: return False
        if read.start < transcript.start or read.stop > transcript.stop :
            return False
        else:
            return compatible(read,transcript,**kwargs)

def translate_fragment_to_transcript_coordinate(read,transcript,**kwargs):
    a=[]
    for i in bam2x.TableIO.parse(read.reads,"bam2bed12",**kwarts):
        a.append(_translate(transcript,i))
    return a
    

def cigar_to_coordinates(cigar,offset=0):
    '''
    demo version
    need to test
    
    deletion from genome (case 3) now consider as exon indel. 
    '''
    exon_starts=[offset]
    exon_lengths=[0]
    state=0
    for i in cigar:
        if i[0]==0 or i[0]==7 or i[0]==8:  # match
           exon_lengths[-1]+=i[1] 
           state=1 
        if i[0]==2:  # deletion from genome , need to consider this should be exon or intron? now count as exon. 
           exon_lengths[-1]+=i[1] 
           state=1 
        if i[0]==3:  # skipped region from the reference 
           if state==1:
               exon_starts.append(exon_starts[-1]+exon_lengths[-1]+i[1]);
               exon_lengths.append(0);
           else:
               exon_starts[-1]+=i[1]
           state=0
    return (exon_starts,exon_lengths)


def parse_string_to_bed(string):
    x=string.split(":")
    if len(x)!=2:
        print >>sys.stderr,"String Format should be\n chromsome:start-end"
        exit(1)
    chr=x[0]
    y=x[1].split("-")
    if len(y)!=2:
        print >>sys.stderr,"String Format should be\n chromsome:start-end"
        exit(1)
    start=int(y[0])-1
    end=int(y[1])
    return BED3(chr,start,end)

def gini_coefficient(iterator):
    '''
    Gini coefficient
    reference:
        http://en.wikipedia.org/wiki/Gini_coefficient
    '''
    l=[i for i in iterator]
    l.sort()
    s1=0.0
    s2=0.0
    y=[]
    n=len(l)
    for i,x in enumerate(l):
        s1+=(i+1)*x
        s2+=x
    if s2==0: return 0.0
    G=2*s1/(n*s2)-float(n+1)/n
    return G


def dis2entropy(iterator):
    s=0
    h={}
    for i in iterator:
        if h.has_key(i):
            h[i]+=1
        else:
            h[i]=1
        s+=1
    e=0.0
    for i in h.values():
        f=float(i)/s
        if f!=0.0:
            e-=f*math.log(f)
    return e






def test():
    logging.basicConfig(level=logging.DEBUG)
    a=Bed("chr1",100,200,"a",0.1,"+")
    b=Bed("chr1",200,300,"b",0.2,"+")
    c=Bed("chr1",400,500,"c",0.2,"+")
    d=merge_bed(a,b)
    e=merge_bed(d,c)
    print "E",e
    print "C and B",merge_bed(c,b)
    print "C and B in E",_translate_to_meta(e,merge_bed(c,b))
    print "C and B reverse translate to CHROMOSOME:",reverse_translate(e,_translate_to_meta(e,merge_bed(c,b)))
    print _translate_to_meta(e,e)
    print compatible(e,e._replace(strand="-"))
    print compatible(e,e._replace(strand="-"),unstranded=True)
if __name__=="__main__":
    test()
