#!/usr/bin/env python
# Programmer : zhuxp
# Date: 
# Last-modified: 01-08-2014, 14:40:41 EST
VERSION="0.1"
import os,sys,argparse
from xplib.Annotation import Bed
from xplib import TableIO,Tools,DBI
from xplib.Tools import IO
import signal
signal.signal(signal.SIGPIPE,signal.SIG_DFL)
import gzip
import time
import xplib.Turing.TuringCodeBook as cb
from xplib.Turing.TuringUtils import *
from xplib.Turing import TuringCode,TuringGraph
from bitarray import bitarray
from multiprocessing import  Pool
import math
'''
Query bam and splicing sites
Construct Turing Graph
and translate reads into BitStr

count reads number in each BitStr
v12:
penalty for increase the number of isoform
v14: 
improve: merging reads with splicing sites first. add mysort

v17: scort by (intron number + 1) * frag_numbers
v18: add rgb color

v20:
    DONE : sort by options
    DONE : bitarray to nice format rep

v20a plan:
    TODO : get coverage and add splicing sites into it. ( memory save not considered , or speed is not considered , just a rash pack)  
    
pv20b :
   parrallel computing
   DONE : splicing sites from query             (integrate xbam2splicingsites)  problem : memory control ? iterate once or twice ? 
   DONE : make it become a function 
pv20c :
    DONE : splicing sites from coverage and sequence        
    TODO : ( add Turing Coverage?)    
TODO.tar.gz:    
    TODO : scan model                            ( get start and end          )  --> xbam2converage splicing sites!
    TODO : start and end sites detection
'''
def sharpScore(l,n):
    '''
    minus nearby max score * fraction
    '''
    length=len(l)
    l2=list(l) #copy l
    print l,n
    for i,x in enumerate(l):
        near=nearby(i,n,length)
        max_nearby=0.0
        max_nearby_index=near[0]
        for j in near:
            if max_nearby < l2[j] and j!=i:
                max_nearby=l2[j]
                max_nearby_index=j
        if max_nearby > l2[i]:
            l[i]=l2[i]-float(n-abs(max_nearby_index-i))/n*max_nearby
        else:
            l[i]=l2[i]+float(n-abs(max_nearby_index-i))/n*max_nearby
        if l[i]<0: 
            l[i]=0.0
    '''
    for i,j in zip(l,l2):
        print i,j
    '''
    return
    
def nearby(i,n,length):
    start=i-n
    if start<0: start=0
    end=i+n
    if end>length: end=length
    return xrange(start,end)
def nearbyIter(i,n,length,include_self=True):
    if include_self:
        yield i
    for j in range(1,n+1):
        if i-j > 0: yield i-j
        if i+j < length : yield i+j 
def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency : xplib')
    p.add_argument('-v','--version',action='version',version='%(prog)s '+VERSION)
    p.add_argument('-i','--input',dest="input",default="stdin",type=str,help="input file DEFAULT: STDIN")
    p.add_argument('-I','--input_format',dest="format",default="guess",type=str,help="input file format")
    p.add_argument('-o','--output',dest="output",type=str,default="stdout",help="output file DEFAULT: STDOUT")
    # p.add_argument('-s','--splicing_sites',dest="splicing_sites",type=str,help="splicing sites [tabix file] [ output of xBamToSplicingSites.py ]")
    p.add_argument('-b','--bam',dest="bam",type=str,help="bam file")
    p.add_argument('-S','--strand',dest="strand",type=str,choices=["read1","read2"],default="read2",help="bam file")
    # p.add_argument('-t','--threshold',dest="threshold",type=float,default=0.95,help="cutoff for interpret fragments percentage, default: 0.95 ")
    p.add_argument('-m','--min_uniq_percentage',dest="min_uniq",type=float,default=0.02,help="min uniq percentage [ add an isoform only if it can interpret more than min uniq percentage fragments ], DEFAULT: %(default)f ")
    p.add_argument('-f','--min_uniq_fpk_increase',dest="min_uniq_fpk_increase",type=float,default=0.2,help="default: %(default)f ")
    p.add_argument('-p','--merge_mismatch_bp',dest="merge_bp",type=int,default=5,help="default: %(default)i ")
    p.add_argument('--sort_model',dest="sort_model",type=int,default=0,choices=[0,1,2],help="model 0: sort by abundance*intron number, model 1: sort by intron number, model 2: sort by abundance, default: %(default)i ")
    p.add_argument('-g','--genome',dest="genome",type=str,help="genome sequence in 2bit format, e.g. mm9.2bit")
    p.add_argument('-n','--num_cores',dest="num_cores",type=int, default=4 ,help="number of cpu cores , DEFAULT: %(default)i")
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()
def Main():
    global args,out,dbi_bam,g, MIN_INTRON_LENGTH, MIN_SPLICING_SITES_SCORE, MIN_FPK_RATIO,query_num
    MIN_INTRON_LENGTH=10
    MIN_SPLICING_SITES_SCORE=2
    '''
    IO TEMPLATE
    '''
    '''
    mySorts={ 0:sort_by_intron_and_abundance,
              1:sort_by_intron,
              2:sort_by_abundance
    }
    '''
    args=ParseArg()
    MIN_FPK_RATIO=args.min_uniq_fpk_increase #TO TEST
    fin=IO.fopen(args.input,"r")
    out=IO.fopen(args.output,"w")
    '''
    END OF IO TEMPLATE 
    '''
    print >>out,"# This data was generated by program ",sys.argv[0]," (version: %s)"%VERSION,
    print >>out,"in bam2x ( https://github.com/nimezhu/bam2x )"
    print >>out,"# Date: ",time.asctime()
    print >>out,"# The command line is :"
    print >>out,"#\t"," ".join(sys.argv)
    # header=["chr","start","end","id","score","strand","seq"];
    # dbi_splicing_sites=DBI.init(args.splicing_sites,"tabix",tabix="metabed",header=header);
    if args.format=="guess":
        args.format=Tools.guess_format(args.input)
    reader=TableIO.parse(fin,args.format)
    query_list=[]
    query_lists=[[] for i in range(args.num_cores)]
    query_num=0
    for i,x in enumerate(reader):
        query_lists[i%args.num_cores].append(x)
    query_num=i+1
    pool=Pool(processes=args.num_cores)
    results=pool.map(querys,query_lists)
    #print results
    output(results)


def querys(listx):
    dbi_bam=DBI.init(args.bam,"bam",method="bam2")
    genome=DBI.init(args.genome,"genome")
    results=[]
    for x in listx:
        results.append(query(x,dbi_bam,genome))
    return results

def output(results):
    for i in iter(results):
        print >>out,i
def iter(results):
    index=-1
    #print "in iter"
    #print "query_num:",query_num
    for i in range(query_num):
        j=i%args.num_cores
        if j==0:
            index+=1
        #print"in iter:",results[j][index]
        yield results[j][index]


def query(i,dbi_bam,genome): # i is query iteem
    h={}
    hc={}
    retv="" #RETURN VALUE
    retv+="BEGIN\nQR\t"+str(i)+"\n"
    '''
    PART A: scoring and adding splicing sites
    '''
    length=len(i)
    l=list()
    l.append(TuringCode(0,cb.ON))
    l.append(TuringCode(0,cb.BLOCKON))
    l.append(TuringCode(length,cb.OFF))
    l.append(TuringCode(length,cb.BLOCKOFF))
    
    array=[0.0 for j in xrange(length)]
    donorSitesScore=[0.0 for j in xrange(length+1)]
    acceptorSitesScore=[0.0 for j in xrange(length+1)]

    seq=genome.get_seq(i).upper()
    
    for j in dbi_bam.query(i,method="bam1",strand=args.strand):
        #print i
        #print j
        read=Tools.translate_coordinates(i,j)
        #print read
        #TODO : get read anc call peaks
        if read.strand=="+":
            for k in read.Exons():
                for k0 in xrange(k.start,k.stop):
                    if k0>=length:
                        break
                    if k0>=0:
                        array[k0]+=1
	    for intron in read.Introns():
	        if len(intron)< MIN_INTRON_LENGTH: continue
                if intron.start > 0  and intron.start <length:
                    donorSitesScore[intron.start]+=1
                if intron.stop > 0  and intron.stop <length:
                    acceptorSitesScore[intron.stop]+=1
    
    gt=[0 for j in range(length)]
    ag=[0 for j in range(length)]
    for j in xrange(length-1):
        if seq[j]=="G" and seq[j+1]=="T":
            gt[j]=1
        elif seq[j]=="A" and seq[j+1]=="G":
            ag[j]=1
    '''
    Score SplicingScores
    '''
    #DONE CLean Scores [ Minus the near by scores]
    sharpScore(donorSitesScore,args.merge_bp)
    sharpScore(acceptorSitesScore,args.merge_bp)
   
    hDonor={}
    hAcceptor={}
    for j in xrange(length):
        if donorSitesScore[j] > MIN_SPLICING_SITES_SCORE:
            l.append(TuringCode(j,cb.BLOCKOFF))
            hDonor[j]=1
        elif acceptorSitesScore[j] > MIN_SPLICING_SITES_SCORE:
            l.append(TuringCode(j,cb.BLOCKON))
            hAcceptor[j]=1
    '''
    closure to score sites by coverage
    '''
    def addCoverageAndSeqScore():
        '''
        two step, 
        1. phase change threshold determin
        2. find the nearest possible splicing sites 
        '''
        diff=[array[0]] # assume the pos -1  is 0.0
        for j in range(length-1):
            diff.append(math.log((array[j+1]+1)/(array[j]+1.0)))
        def mean_sd(indexes):
            n=0
            s=0.0
            sd=0.0
            a=[]
            for k in indexes:
                n+=1
                max_diff=0.0
                for k0 in nearbyIter(k,args.merge_bp,length):
                    if abs(max_diff) < abs(diff[k0]):
                        max_diff=diff[k0]
                s+=max_diff
                a.append(max_diff)
            if n>0: 
                mean=s/n
            else:
                return 0.0,0.0
            for k in a:
                sd+=(k-mean)*(k-mean)
            if n > 1:
                sd=sd/(n-1)
            if n > 1:
                sd=math.sqrt(sd/(n-1))
            else:
                sd=0.0
            return mean,sd
        (meanDonor,sdDonor)=mean_sd(hDonor.keys())
        (meanAcceptor,sdAcceptor)=mean_sd(hAcceptor.keys())
        #print meanDonor,sdDonor
        #print meanAcceptor,sdAcceptor
        '''
        adding possible splicing sites based on diff value
        '''
        def correctToNearAcceptor(j):
            signal=0
            for k in nearbyIter(j,args.merge_bp,length):
                if hAcceptor.has_key(k): 
                    return k
            for k in nearbyIter(j,args.merge_bp,length):
                if ag[k]==1:
                    return k+2
            return -1

        def correctToNearDonor(j):
            signal=0
            for k in nearbyIter(j,args.merge_bp,length):
                if hDonor.has_key(k): 
                    return k
            for k in nearbyIter(j,args.merge_bp,length):
                if gt[k]==1:
                    return k
            return -1

        for j in range(length):
            #if diff[j]  > sdAcceptor + meanAcceptor:
            if diff[j]  > meanAcceptor + sdAcceptor and diff[j] > 0.0:
                k=correctToNearAcceptor(j)
                if k!=-1 and not hAcceptor.has_key(k):
                    l.append(TuringCode(k,cb.BLOCKON))


            #if  diff[j] < meanDonor - sdDonor:
            if  diff[j] < meanDonor -  sdDonor and diff[j] < 0.0:
                k=correctToNearDonor(j)
                if k!=-1 and not hDonor.has_key(k):
                    l.append(TuringCode(k,cb.BLOCKOFF))
    addCoverageAndSeqScore()
    l.sort()
    '''
    PART B : Report wig
    '''
    retv+="WIG"+"\n"
    retv+= "index\tnt\tcoverage\tdonorSitesScore\tacceptorSitesScore\tGT\tAG\n"
    for j in xrange(len(i)):
        retv+= str(j)+"\t"+str(seq[j])+"\t"+str(array[j])+"\t"+str(donorSitesScore[j])+"\t"+str(acceptorSitesScore[j])+"\t"+str(gt[j])+"\t"+str(ag[j])+"\n"
    '''
    end of adding wig
    '''
    
    
    
    '''
    PART C: GREEDY ADDING ISOFORMS
    
    '''
    
    g=TuringGraph(l)
    bitarray_path=bitarray(2*len(g))
    bitarray_path.setall(True)
    #print >>sys.stderr,"processing",i0," entry:",i
    retv+="FIGURE\t"+g.graph_str(600)+"\n"
    paths_number=g.count_paths_number()
    retv+="PATH_NUMBER\t"+str(paths_number)
    print >>sys.stderr,"path number:",paths_number
    h={}
    hc={}
    j0=0;
    total_frag=0
    for j in dbi_bam.query(i,method="bam2",strand=args.strand):
        p=[]
        for k in j:
            p.append(TuringFactory(Tools.translate_coordinates(i,k)))
        a=g.translate_paths_into_bits(p,args.merge_bp)
        if isSelfIncompatible(a): continue
        if h.has_key(a.tobytes()):
            h[a.tobytes()]+=1
            if hc[a.tobytes()]!=a:
                print >>sys.stderr,"WARNING"
        else:
            h[a.tobytes()]=1
            hc[a.tobytes()]=a
        j0+=1
    retv+="PATTERN_NUMBER\t"+str(len(h.keys()))+"\n"
    retv+="FRG_NUMBER\t"+str(j0)+"\n"
    total_frag=j0
    print >>sys.stderr,"fragments number:",j0
    print >>sys.stderr,"path pattern number:",len(h.keys())
    if j0==0: return retv
    #sorted_keys=sorted(h,key=h.get,reverse=True)    
    #TODO 
    #CLOUSURE
    # sorted_keys=sorted(h.keys(),cmp=mySorts[args.sort_model],reverse=True)    
    sorted_keys=h.keys()
    if args.sort_model==0:
        sorted_keys.sort(lambda x,y:bitarray_to_intron_number(hc[x])*h[x]-bitarray_to_intron_number(hc[y])*h[y]  or h[x]-h[y], reverse=True)
    elif args.sort_model==1:
        sorted_keys.sort(lambda x,y:bitarray_to_intron_number(hc[x])-bitarray_to_intron_number(hc[y]) or h[x]-h[y], reverse=True)
    else:
        sorted_keys.sort(lambda x,y:h[x]-h[y] or bitarray_to_intron_number(hc[x])-bitarray_to_intron_number(hc[y]), reverse=True) 
    clique=[]
    cliques=[clique]
    bits=bitarray(len(g)*2)
    bits.setall(True)
    cliques_pattern=[bits]
    for j,key in enumerate(sorted_keys):
        retv+="No."+str(j)+"\t"+str(bitarray_to_rep(hc[key]))+" "+str(h[key])+" "+str(bitarray_to_intron_number(hc[key]))+"\n"
        joined_clique=False
        for m,clique in enumerate(cliques):
            if isCompatible(cliques_pattern[m],hc[sorted_keys[j]]):
                joined_clique=True
                cliques[m].append(j)
                cliques_pattern[m]=bitarray_and(cliques_pattern[m],hc[sorted_keys[j]])
                break
        if not joined_clique:
            clique=[]
            bits=bitarray(len(g)*2)
            bits.setall(True)
            max_index=0
            clique.append(j)
            bits=bitarray_and(bits,hc[sorted_keys[j]])
            cliques.append(clique)
            cliques_pattern.append(bits)
    retv+="OUTPUT\n"

    cumulative_score=0
    #update all cliques pattern.
    j0=0
    #buffer_keys=[] # those key are ignored in previous uniqs
    max_uniq_fpk=0.0
    max_uniq=0
    for j,x in enumerate(cliques):
        score=0
        c=[] 
        for k,y in enumerate(sorted_keys):
            if isCompatible(cliques_pattern[j],hc[y]):
                cliques_pattern[j]=bitarray_and(cliques_pattern[j],hc[y])
                score+=h[y]
                c.append(k)
        cliques_pattern[j][-1]=True
        cliques_pattern[j][-2]=True
        uniq_score=0
        for k,y in enumerate(x):
            uniq_score+=h[sorted_keys[y]]
        '''
        for k,y in enumerate(buffer_keys):
            if isCompatible(cliques_pattern[buffer_keys[k]],
                cliques[j].append(k)
        '''
        if float(uniq_score)/total_frag < args.min_uniq :
            print >>sys.stderr,"ignore due to small uniq frags:",float(uniq_score)/total_frag
            continue
        else:
            bed=g.translate_bits_into_bed(cliques_pattern[j])
            cdna_length=bed.cdna_length()
            if cdna_length==0:
                print >>sys.stderr,"ignore due to cdna_length",bed
                continue
            uniq_fpk=float(uniq_score)/cdna_length*1000.0
            if  uniq_fpk > max_uniq_fpk:
                max_uniq_fpk=uniq_fpk
            if uniq_fpk  < max_uniq_fpk * MIN_FPK_RATIO:
                print >>sys.stderr,"ignore due to fpk",bed,"\tuniq_fpk:",uniq_fpk,"\t current max:",max_uniq_fpk
                continue
            rgb=255-uniq_score*240/total_frag  
            cumulative_score+=uniq_score
            j0+=1
            retv+="NO."+str(j0)+" CLIQUE"+"\n"    
            retv+="CLIQUE\t"+str(c)+"\nUNIQ\t"+str(cliques[j])+"\nPATTERN\t"+bitarray_to_rep(cliques_pattern[j])+"\n"
            #retv+="CLIQUE\t",c,"\nUNIQ\t",cliques[j],"\nPATTERN\t",cliques_pattern[j]
            
            bed=g.translate_bits_into_bed(cliques_pattern[j])
            bed.score=score*1000.0/bed.cdna_length()
            bed.chr=i.id
            bed.id=i.id+"_"+"NO."+str(j0)
            bed.itemRgb=str(rgb)+","+str(rgb)+","+str(rgb)
            retv+="UNIQ_FRG\t"+str(uniq_score)+"\n"
            retv+="TOTAL_FRG\t"+str(score)+"\n"
            retv+="UNIQ_FPKi\t"+str(uniq_fpk)+"\n"
            retv+="FPK(SCORE)\t"+str(bed.score)+"\n"
            retv+="TR\t"+str(bed)+"\n"

            retv+="UNIQ\t"+str(uniq_score)+"\tBED\t"+str(Tools.translate_coordinates(i,bed,True))+"\n"
            #print >>sys.stderr,"UNIQ_FPK/CURRENT_MAX_UNIQ_FPK",uniq_fpk/max_uniq_fpk

            retv+="\n"
    retv+="INTEPRET FRG\t"+str(float(cumulative_score)/total_frag)+"\n"
    print >>sys.stderr,"INTEPRET FRG\t",float(cumulative_score)/total_frag
    retv+="END\n"
    retv+="\n"
    retv+="\n"
    return retv

    
'''        
def sort_by_intron_and_abundance(x,y,h,hc):
    return bitarray_to_intron_number(hc[x])*h[x]-bitarray_to_intron_number(hc[y])*h[y]  or h[x]-h[y]
def sort_by_intron(x,y):
    return bitarray_to_intron_number(hc[x])-bitarray_to_intron_number(hc[y]) or h[x]-h[y]
def sort_by_abundance(x,y):
    return  h[x]-h[y] or bitarray_to_intron_number(hc[x])-bitarray_to_intron_number(hc[y]) 
'''



    

def bedToID(bed):
    s=bed.chr+"\t"+str(bed.start)+"\t"+str(bed.stop)+"\t"+bed.strand
    return s;

if __name__=="__main__":
    Main()

