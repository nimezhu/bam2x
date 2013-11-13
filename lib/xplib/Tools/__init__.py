# Programmer : zhuxp
# Date:  Sep 2012
# Last-modified: 11-13-2013, 17:17:12 EST
from string import upper,lower
from xplib.Annotation import Fragment,Bed,Bed12,GeneBed
import xplib
import copy

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
def translate_coordinate(A,B):
    '''
    translate B's coordiante based on A
    '''
    if A.chr!=B.chr: return None
    if A.strand=="+" or A.strand==".":
        return (B.start-A.start,B.stop-A.start,B.strand)
    if A.strand=="-":
        strand="."
        if B.strand=="-":strand="+"
        if B.strand=="+":strand="-"
        return (A.stop-B.stop,A.stop-B.start,strand)
def translate_coordinates(A,B): # B is Bed12 format
    '''
    Translate Bed12 Object B based on simple annotation A
    '''
    chr=A.id
    id=B.id
    (start,stop,strand)=translate_coordinate(A,B)
    score=B.score

    if isinstance(B,Bed12):
        (cds_start,cds_stop,cds_strand)=translate_coordinate(A,Bed([B.chr,B.cds_start,B.cds_stop]))
        itemRgb=B.itemRgb
        blockCount=B.blockCount
        blockSizes=copy.copy(B.blockSizes)
        if A.strand=="-": 
            blockSizes=blockSizes[::-1]
        blockStarts=copy.copy(B.blockStarts)
        if A.strand=="+" or A.strand==".":
            for i,x in enumerate(blockStarts):
                blockStarts[i]=blockStarts[i]
        elif A.strand=="-":
            for i,x in enumerate(blockStarts):
                blockStarts[i]+=blockSizes[i]
            blockStarts=blockStarts[::-1]
            for i,x in enumerate(blockStarts):
                print blockStarts[i]
                blockStarts[i]=B.stop-(blockStarts[i]+B.start)
        return Bed12([chr,start,stop,id,score,strand,cds_start,cds_stop,itemRgb,blockCount,blockSizes,blockStarts])
    elif isinstance(B,GeneBed):
        C1=Bed12(B.toBedString())
        C2=translate_coordinates(A,C1)
        return GeneBed(C2.toGenePredString())
    else:
        C=copy.copy(B)
        C.chr=chr
        C.start=start
        C.stop=stop
        C.end=stop
        C.strand=strand
        if B.__dict__.has_key('pos'):
            C.pos=C.start+1  # position is 1  index
        return C
        

        




def overlap(A,B):
    '''
    if A is overlapping with B.
    A and B are ? extends Bed class.
    '''
    if(A.chr != B.chr) : return 0
    if (A.stop < B.start) : return 0
    if (B.stop < A.start) : return 0
    return 1
from xplib.Annotation import Bed
def find_nearest(bed,dbi,extends=50000,**dict):
    start=bed.start-extends
    stop=bed.stop+extends
    chr=bed.chr
    if start<0: start=0
    new_bed=Bed([chr,start,stop])

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
    VERSION: TEST
    if a and b are compatible return true
    a and b are BED12 class or GENEBED 
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
        for i0,i in enumerate(xplib.TableIO.parse(read.reads,"bam2bed12",**kwargs)):
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
    return Bed([chr,start,end])

def suffix(string):
    x=string.split('.')
    return x[-1]


def guess_format(string):
    
    x=string.split(".")
    suffix=lower(x[-1])
    if x[-1]=="gz":
        suffix=lower(x[-2])
    if suffixToFormat.has_key(suffix):
        return suffixToFormat[suffix]
    else:
        return suffix

