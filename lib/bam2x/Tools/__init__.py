# Programmer : zhuxp
# Date:  Sep 2012
# Last-modified: 02-14-2014, 17:25:39 EST
from string import upper,lower
from bam2x.Annotation import BED6 as Bed
from bam2x.Annotation import BED12 as Bed12
from bam2x.Annotation import BED12 
from bam2x.Annotation import BED3
import bam2x
import copy
import logging
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
    m=abs(A.start-B.start)
    if( m > abs(A.start-B.stop)): m = abs(A.start-B.stop)
    if( m > abs(A.stop-B.stop)): m = abs(A.stop-B.stop)
    if( m > abs(A.stop-B.start)): m = abs(A.stop-B.start)
    return m
def translate_coordinate(coord,bed,reverse=False):
    '''
    translate bed's coordiante based on coord
    if reverse is True:
        bed is in coord's coordinates and translate it back to the chromosome coordinates
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
    Translate Bed12 Object bed based on simple annotation coord

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
def _merge_bed6(beds,id="noname"):
    '''
    a simple turing state change
    '''
    l=[]

    for i in beds:
        l.append((i.start,-1))
        l.append((i.stop,1))
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
    return BED12(chr,start,end,id,0.0,".",start,start,"0,0,0",blockCount,blockSizes,blockStarts) 

def translate(bedA,bedB):
    #TODO
    '''
    return merge bed
    and 
    the location for bedA in merge bed
    the location for bedB in merge bed
    return (merged_bed, bedA_in_merged_bed , bedB_in_merged_bed)
    '''
    meta=merge_bed(bedA,bedB,bedA.id+"_"+bedB.id+"_merged")
    new_bedA=_translate_to_meta(meta,bedA)
    new_bedB=_translate_to_meta(meta,bedB)
    return new_bedA,new_bedB,meta
def _translate_to_meta(meta,bed):
    '''
    another simple turing 
    meta is the bigger one
    meta code 1
    bed code 2
    '''
    l=[]
    for i in meta.Exons():
        l.append((i.start,-1,1))
        l.append((i.stop,1,1))
    for i in bed.Exons():
        l.append((i.start,-1,2))
        l.append((i.stop,1,2))
    l.append((bed.cds_start,0,3))
    l.append((bed.cds_stop,0,4))
    l.sort()
    logging.debug("bed id %s, l=%s",bed.id,l)
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
                logging.debug("start:%i",i[0])
                b_start=i[0]
                switch=1
                blockStarts.append(b_start)
        elif switch==1:
            if state == 0:
                logging.debug("end:%i",i[0])
                b_end=i[0]
                switch=0
                blockSizes.append(b_end-b_start)
    blockCount=len(blockSizes)
    return BED12(meta.id,blockStarts[0],blockStarts[-1]+blockSizes[-1],"id",0.0,bed.strand,cds_start,cds_stop,"0,0,0",blockCount,blockSizes,blockStarts)
    
    


    



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
        if distance(bed,result)<d:
            d=distance(bed,result)
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


def compatible(a,b,**kwargs):
    '''
    TODO: UNFINISHED IN Bam2x LIB
    VERSION: TEST
    if a and b are compatible return true
    a and b are BED12 class 
    definition of compatible
       the overlap region should be same transcript structure.
    '''
    if (not overlap(a,b)): return True; 
    if (a.strand!=b.strand): return False;
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
    a_starts_slice=[];
    #TODO REVISE IT ! 
    for i in a.exon_starts:
        if i>start and i<stop:
            a_starts_slice.append(i); 
    j=0;
    l=len(a_starts_slice);
    for i in b.exon_starts:
        if i>start and i<stop:
            if (j>l-1 or a_starts_slice[j]!=i) :
                return False
            j+=1
    if j!=l : 
        return False
    
    
    a_stops_slice=[];
    for i in a.exon_stops:
        if i>start and i<stop:
            a_stops_slice.append(i); 
    
    j=0;
    l=len(a_stops_slice);
    for i in b.exon_stops:
        if i>start and i<stop:
            if (j>l-1 or a_stops_slice[j]!=i) : 
                return False
            j+=1
    if j!=l : 
        return False
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
        if kwargs.has_key("strand") and kwargs["strand"]==True:
            if read.strand != transcript.strand: return False
        if read.start < transcript.start or read.stop > transcript.stop :
            return False
        else:
            return compatible(read,transcript,**kwargs)
    

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


def test():
    logging.basicConfig(level=logging.DEBUG)
    a=Bed("chr1",100,200,"a",0.1,"+")
    b=Bed("chr1",200,300,"b",0.2,"-")
    c=Bed("chr1",400,500,"c",0.2,"-")
    d=merge_bed(a,b)
    e=merge_bed(d,c)
    print _translate_to_meta(e,merge_bed(c,b))
    print _translate_to_meta(e,d)
if __name__=="__main__":
    test()
