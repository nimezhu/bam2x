#!/usr/bin/env python
# Programmer : zhuxp
# Date: 
# Last-modified: 02-18-2014, 14:33:50 EST
VERSION="bam2x 0.1.1"
import os,sys,argparse
from bam2x.Annotation import BED6 as Bed
from bam2x.Annotation import BED12
from bam2x import TableIO,Tools,DBI,IO
import signal
signal.signal(signal.SIGPIPE,signal.SIG_DFL)
import gzip
import time
import bam2x.Turing.TuringCodeBook as cb
from bam2x.Turing.TuringUtils import *
from bam2x.Turing import TuringCode,TuringGraph
from bitarray import bitarray
from multiprocessing import  Pool
import math
import itertools
import logging
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
pv20d :
    DONE : using Turing Tuple Code instead of Turing Code
    DONE : using Bed Tuple IO ,instead of Bed IO, write Bed Tuple IO.
    DONE : simple query bam paired end (ignore the reads are not in this region)
pv20e :
    DONE: change result to json or xml
pv20f: 
    DONE : fix bugs of same location have two or more codes , revise codes length into correct one using hash [ only revise the tuple version now ] 
    DONE : strand correction
pv20g:
    DONE : adding exon structure in bam2peak?
pv20h:
    DONE : start with the exon structure.
    DONE : fix the some result entry don't have pattern error (ignore these entries)
    DONE : increase the dark of output
    DONE : ignore the pattern that is not compatible with query bed12
    DONE : add options to report seq or not
pv22a:
    DONE : report singleton as singleton , score is still cliques score
    TODO : link undetected exon.
v0.0.4:
    DONE: rm xplib and convert to bam2x lib, 
    DONE: using turing namedtuple instead of tuple
    DONE: no bam2 query but bam1 query ( fragment query in lib )
v0.1.1:
    DONE: merge only overlapping reads.  link stragdy is changed
'''
def help():
    return "detect isoforms."
def report_format(x,**kwargs):
    #TOD
    #logging.debug(x)
    try:
        S=""
        S+="QR\t"+str(x["QR"])+"\n"
        S+="FIGURE %s\n"%x["FIGURE"]
        if args.report_seq:
            if x["QR"][2]-x["QR"][1] < 100:
                S+="SQ\t"+x["WIG_TABLE"][0]
            else:
                S+="PARTIAL_SQ\t"+x["WIG_TABLE"][0][0:100]+"......."
            S+="\n"
        S+="FRG_NUMBER %s\n"%x["FRG_NUMBER"]
        S+="PATTERN_NUMBER %s\n"%x["PATTERN_NUMBER"]
        S+="  INDEX\tPATTERN\tFRGs\tINTRON_NUMBER\n"
        for i,y in enumerate(x["PATTERNS"]):
            S+="  "
            S+="No.%s\t%s\t%s\t%s\n"%(i,y[0],y[1],y[2])
        S+="CLIQUES\n"
        for i,y in enumerate(x["CLIQUES"]):
            S+="  NO."+str(i+1)+"\n"
            S+="    REP"+"\t"+y["REP"]+"\n"
            S+="    UNIQ_SCORE\t"+str(y["UNIQ_SCORE"])+"\n"
            if y.has_key("BED"):
                S+="    BED"+"\t"+str(y["BED"])+"\n"
            elif y.has_key("BEDS"):
                for bed in y["BEDS"]:
                    S+="    BED"+"\t"+str(bed)+"\n"
            if args.report_seq:
                S+="    ALL_COMPATIBLE_PATTERN\t"+str(y["ALL_GROUP"])+"\n"
                S+="    UNIQ_COMPATIBLE_PATTERN\t"+str(y["UNIQ_GROUP"])+"\n"
        S+="INTERPRET\t"+str(x["INTERPRET_FRG"])+"\n"
        S+="//\n"
        return S
    except:
        raise
        if x.has_key("QR"):
            return "# ERROR RESULT FORMAT FOR "+str(x["QR"])
        else:
            return "# UNKNOWN QR"
    


def sharpScore(l,n):
    '''
    minus nearby max score * fraction
    '''
    length=len(l)
    l2=list(l) #copy l
    #print l,n
    for i,x in enumerate(l):
        #print i,x
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
    #print "end"
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
def set_parser(p):
    ''' This Function Parse the Argument '''
    p.add_argument('-I','--input_format',dest="format",default="bed12",choices=["bed6","bed12"],type=str,help="input file format DEFAULT:%(default)s , all possible exons")
    p.add_argument('-b','--bam',dest="bam",type=str,help="bam file")
    p.add_argument('-S','--strand',dest="strand",type=str,choices=["read1","read2"],default="read2",help="strand information")
    p.add_argument('-m','--min_uniq_percentage',dest="min_uniq",type=float,default=0.02,help="min uniq percentage [ add an isoform only if it can interpret more than min uniq percentage fragments ], DEFAULT: %(default)f ")
    p.add_argument('-f','--min_uniq_fpk_increase',dest="min_uniq_fpk_increase",type=float,default=0.05,help="default: %(default)f ")
    p.add_argument('-p','--merge_mismatch_bp',dest="merge_bp",type=int,default=5,help="default: %(default)i ")
    p.add_argument('--sort_model',dest="sort_model",type=int,default=0,choices=[0,1,2],help="model 0: sort by abundance*intron number, model 1: sort by intron number, model 2: sort by abundance, default: %(default)i ")
    p.add_argument('-g','--genome',dest="genome",type=str,help="genome sequence in 2bit format, e.g. mm9.2bit")
    p.add_argument('--report_seq',default=False,dest="report_seq",action="store_true",help="report sequence [ warning: cost mem ]")
def run(local_args):
    logging.basicConfig(level=logging.WARNING)
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
    args=local_args
    #print "debug:",args.report_seq
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
    #if args.format=="guess":
    #    args.format=IO.guess_format(args.input)
    reader=TableIO.parse(fin,args.format)
    query_list=[]
    query_lists=[[] for i in range(args.num_cpus)]
    query_num=0
    for i,x in enumerate(reader):
        query_lists[i%args.num_cpus].append(x)
    query_num=i+1
    #querys(query_lists[0])  #DEBUG
    pool=Pool(processes=args.num_cpus)
    results=pool.map(querys,query_lists)
    #print results
    output(results)


def querys(listx):
    dbi_bam=DBI.init(args.bam,"bam")
    genome=DBI.init(args.genome,"genome")
    results=[]
    for x in listx:
        results.append(query(x,dbi_bam,genome))
    return results

def output(results):
    for i in iter(results):
        print >>out,report_format(i,indent="  ")
def iter(results):
    index=-1
    #print "in iter"
    #print "query_num:",query_num
    for i in range(query_num):
        j=i%args.num_cpus
        if j==0:
            index+=1
        #print"in iter:",results[j][index]
        yield results[j][index]


def query(i,dbi_bam,genome): # i is query iteem
    h={}
    hc={}
    #retv="" #RETURN VALUE
    #retv+="BEGIN\nQR\t"+str(i)+"\n"
    x=i
    ret_dict={}
    ret_dict["QR"]=i
    '''
    PART A: scoring and adding splicing sites
    '''
    print >>sys.stderr,"QR",str(i)
    length=i.stop-i.start
    g=list()
    g.append(TuringCode(0,cb.ON))
    g.append(TuringCode(0,cb.BLOCKON))
    g.append(TuringCode(length,cb.OFF))
    g.append(TuringCode(length,cb.BLOCKOFF))
    #logging.debug(g)
    array=[0.0 for j in xrange(length)]
    donorSitesScore=[0.0 for j in xrange(length+1)]
    acceptorSitesScore=[0.0 for j in xrange(length+1)]
    seq=genome.get_seq(i).upper()
    for j in dbi_bam.query(i,method="bam1",strand=args.strand):
        read=Tools.translate_coordinates(i,j)
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
    hAcceptor[0]=1
    hDonor[length]=1
    for j in xrange(length):
        if donorSitesScore[j] > MIN_SPLICING_SITES_SCORE:
            g.append(TuringCode(j,cb.BLOCKOFF))
            hDonor[j]=1
        elif acceptorSitesScore[j] > MIN_SPLICING_SITES_SCORE:
            g.append(TuringCode(j,cb.BLOCKON))
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
        #logging.debug(length)
        #logging.debug(g)
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
                    if k0>=length: continue
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
                try:
                    if ag[k]==1:
                        return k+2
                except:
                    pass
            return -1
        def correctToNearDonor(j):
            signal=0
            for k in nearbyIter(j,args.merge_bp,length):
                if hDonor.has_key(k): 
                    return k
            for k in nearbyIter(j,args.merge_bp,length):
                try:
                    if gt[k]==1:
                        return k
                except:
                    pass
            return -1


        for j in range(length):
            #if diff[j]  > sdAcceptor + meanAcceptor:
            if diff[j]  > meanAcceptor + sdAcceptor and diff[j] > 0.0:
                k=correctToNearAcceptor(j)
                if k!=-1 and not hAcceptor.has_key(k):
                    hAcceptor[k]=1
                    g.append(TuringCode(k,cb.BLOCKON))
                # add to revise small fragrments ?
                '''
                elif k==-1 and diff[j] > meanAcceptor + 3 * sdAcceptor and not hAcceptor.has_key(j):
                    hAcceptor[j]=1
                    l.append((j,cb.BLOCKON))
                '''


            #if  diff[j] < meanDonor - sdDonor:
            if  diff[j] < meanDonor -  sdDonor and diff[j] < 0.0:
                k=correctToNearDonor(j)
                if k!=-1 and not hDonor.has_key(k):
                    hDonor[k]=1
                    g.append(TuringCode(k,cb.BLOCKOFF))
                '''
                elif k==-1 and diff[j] <  meanDonor - 3 * sdDonor and not hDonor.has_key(j):
                    hDonor[j]=1
                    l.append((j,cb.BLOCKOFF))
                '''
        if isinstance(x,BED12): #strange could not use i?
            ti=Tools.translate_coordinates(x,x)
            #print "debug ti:",tia
            ti_codes=[]
            for j in ti.blockStarts:
                k=correctToNearAcceptor(j)
                if k==-1:
                    ti_codes.append(TuringCode(j,cb.BLOCKON))
                else:
                    ti_codes.append(TuringCode(k,cb.BLOCKON))
                if k!=-1 and not hAcceptor.has_key(k):
                    hAcceptor[k]=1
                    g.append(TuringCode(k,cb.BLOCKON))
                if k==-1 and not hAcceptor.has_key(j):
                    hAcceptor[j]=1
                    g.append(TuringCode(j,cb.BLOCKON))
            last_stop=0
            for j,m in itertools.izip(ti.blockStarts,ti.blockSizes):
                blockStop=j+m
                #print "debug blockStop",blockStop
                k=correctToNearDonor(blockStop)
                if k==-1:
                    ti_codes.append(TuringCode(blockStop,cb.BLOCKOFF))
                    last_stop=j
                else:
                    ti_codes.append(TuringCode(k,cb.BLOCKOFF))
                    last_stop=k
                if k!=-1 and not hDonor.has_key(k):
                    hDonor[k]=1
                    g.append(TuringCode(k,cb.BLOCKOFF))
                elif k==-1 and not hDonor.has_key(blockStop):
                    hDonor[blockStop]=1
                    g.append(TuringCode(blockStop,cb.BLOCKOFF))

            ti_codes.append(TuringCode(correctToNearAcceptor(ti.blockStarts[0]),cb.ON))
            ti_codes.append(TuringCode(last_stop,cb.OFF))
            ti_codes.sort()
            #print "debug ti_codes",ti_codes
            #print ti_codes
            #TODO  GENERate start path.
            # l.sort()
            graph=TuringGraph(g)
            g_len=len(graph)
            initial_bits=graph.translate_path_into_bits(TuringGraph(ti_codes),args.merge_bp)
            #print "debug ", initial_bits
            for i in range(0,2*g_len,2):
                if initial_bits[i]==True and initial_bits[i+1]==False:
                    initial_bits[i+1]=True
            initial_bits[-1]=True
            initial_bits[-2]=True
            return initial_bits,graph
        # l.sort()
        return None,TuringGraph(g)




            


    
    initial_bits,g=addCoverageAndSeqScore()
    '''
    PART B : Report wig
    retv+="WIG"+"\n"
    retv+= "index\tnt\tcoverage\tdonorSitesScore\tacceptorSitesScore\tGT\tAG\n"
    for j in xrange(len(i)):
        retv+= str(j)+"\t"+str(seq[j])+"\t"+str(array[j])+"\t"+str(donorSitesScore[j])+"\t"+str(acceptorSitesScore[j])+"\t"+str(gt[j])+"\t"+str(ag[j])+"\n"
    '''
    #logging.debug(initial_bits)
    if args.report_seq:
        ret_dict["WIG_TABLE"]=(seq,array,donorSitesScore,acceptorSitesScore,gt,ag)
    
    '''
    end of adding wig
    '''
    
    
    
    '''
    PART C: GREEDY ADDING ISOFORMS
    
    '''
    
    
    # bitarray_path=bitarray(2*len(g))
    # bitarray_path.setall(True)
    #TODO change bitarray?
    #logging.debug(g)
    ret_dict["FIGURE"]=g.graph_str(600)
    h={}
    hc={}
    j0=0;
    total_frag=0
    #for j in g: print "debug g:",j
    #print "INIT",bitarray_to_rep(initial_bits)
    g_len=len(g)
    if initial_bits is None:
        initial_bits=bitarray(2*g_len)
        initial_bits.setall(True)

    for j in dbi_bam.query(i,method="bam2fast",strand=args.strand):
        p=[]
        #logging.debug(j)
        #print "debug",j
        
        if j[0].strand!=i.strand: continue
        for k in j:
            p.append(TuringFactory(Tools.translate_coordinates(i,k))) #TODO check this
        # p.append(TuringFactory(Tools.translate_coordinates(i,j))) #TODO check this

        #print "debug glen:",g_len

        a=g.translate_paths_into_bits(p,args.merge_bp)
        #print "debug bits: ",a
        #print "debug bed :",translate_bits_into_bed(g,a)
        if isSelfIncompatible(a): continue
        if h.has_key(a.tobytes()):
            h[a.tobytes()]+=1
            if hc[a.tobytes()]!=a:
                print >>sys.stderr,"WARNING"
        else:
            h[a.tobytes()]=1
            hc[a.tobytes()]=a
        j0+=1
     #retv+="PATTERN_NUMBER\t"+str(len(h.keys()))+"\n"
    ret_dict["PATTERN_NUMBER"]=len(h.keys())
    # print >>sys.stderr,"debug pattern number",ret_dict["PATTERN_NUMBER"]
    #retv+="FRG_NUMBER\t"+str(j0)+"\n"
    ret_dict["FRG_NUMBER"]=j0
    total_frag=j0
    if j0==0: return ret_dict
    sorted_keys=h.keys()
    if args.sort_model==0:
        sorted_keys.sort(lambda x,y:bitarray_to_intron_number(hc[x])*h[x]-bitarray_to_intron_number(hc[y])*h[y]  or h[x]-h[y], reverse=True)
    elif args.sort_model==1:
        sorted_keys.sort(lambda x,y:bitarray_to_intron_number(hc[x])-bitarray_to_intron_number(hc[y]) or h[x]-h[y], reverse=True)
    else:
        sorted_keys.sort(lambda x,y:h[x]-h[y] or bitarray_to_intron_number(hc[x])-bitarray_to_intron_number(hc[y]), reverse=True) 
    clique=[]
    cliques=[]
    bits_no_prior=bitarray(g_len*2)
    bits_no_prior.setall(True)
    cliques_pattern=[]
    ret_dict["PATTERNS"]=list()

    keys_len=len(sorted_keys)
    
    
    joined_cliques=[False for i0 in range(keys_len)]
    unjoined_pattern=range(0,keys_len)

    for j,key in enumerate(sorted_keys):
        ret_dict["PATTERNS"].append((bitarray_to_rep(hc[key]),h[key],bitarray_to_intron_number(hc[key])))
        if joined_cliques[j] : continue
        joined_clique=False
        joined_cliques[j]=True
        bits_no_prior=bitarray(g_len*2)
        bits_no_prior.setall(True)  
        h_clique={}
        if isCompatible(initial_bits,hc[sorted_keys[j]]):
            max_index=0
            clique=[]
            clique.append(j)
            cliques.append(clique)
            cliques_pattern.append(hc[sorted_keys[j]])
            h_clique[j]=1
            update=True
            while update:
                update=False
                logging.debug("update:%s",update)
                for k0,k in enumerate(unjoined_pattern):
                    logging.debug("del %i",k)
                    if k==j:
                        del unjoined_pattern[k0]
                    elif isOverlapCompatible(cliques_pattern[-1],hc[sorted_keys[k]]):
                        update=True
                        cliques[-1].append(k)
                        cliques_pattern[-1]=bitarray_and(cliques_pattern[-1],hc[sorted_keys[k]])
                        joined_cliques[k]=True
                        del unjoined_pattern[k0]
                        h_clique[k]=1
                        break
                if update==False:
                    if k==j:
                        continue
                    for k in range(keys_len):
                        if h_clique.has_key(k): continue
                        if isOverlapCompatible(cliques_pattern[-1],hc[sorted_keys[k]]):
                            update=True
                            cliques[-1].append(k)
                            cliques_pattern[-1]=bitarray_and(cliques_pattern[-1],hc[sorted_keys[k]])
                            h_clique[k]=1
                            break



    ret_dict["CLIQUES"]=list()

    cumulative_score=0
    #update all cliques pattern.
    j0=0
    #buffer_keys=[] # those key are ignored in previous uniqs
    max_uniq_fpk=0.0
    max_uniq=0
    for j,x in enumerate(cliques):
        score=0
        
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
            #print >>sys.stderr,"ignore due to small uniq frags:",float(uniq_score)/total_frag
            continue
        else:
            logging.debug(cliques_pattern)
            bed=g.translate_bits_into_bed(cliques_pattern[j])
            cdna_length=bed.cdna_length()
            if cdna_length==0:
                #print >>sys.stderr,"ignore due to cdna_length",bed
                continue
            uniq_fpk=float(uniq_score)/cdna_length*1000.0
            if  uniq_fpk > max_uniq_fpk:
                max_uniq_fpk=uniq_fpk
            if uniq_fpk  < max_uniq_fpk * MIN_FPK_RATIO:
                #print >>sys.stderr,"ignore due to fpk",bed,"\tuniq_fpk:",uniq_fpk,"\t current max:",max_uniq_fpk
                continue
            rgb=200-uniq_score*200/total_frag  
            cumulative_score+=uniq_score
            j0+=1
            #retv+="NO."+str(j0)+" CLIQUE"+"\n"    
            #retv+="CLIQUE\t"+str(c)+"\nUNIQ\t"+str(cliques[j])+"\nPATTERN\t"+bitarray_to_rep(cliques_pattern[j])+"\n"
            ret_dict["CLIQUES"].append(dict())
            #ret_dict["CLIQUES"][-1]["CLIQUE"]=(c,cliques[j],bitarray_to_rep(cliques_pattern[j]))
            ret_dict["CLIQUES"][-1]["ALL_GROUP"]=cliques[j]
            ret_dict["CLIQUES"][-1]["UNIQ_GROUP"]=cliques[j]
            ret_dict["CLIQUES"][-1]["REP"]=bitarray_to_rep(cliques_pattern[j])
            #retv+="CLIQUE\t",c,"\nUNIQ\t",cliques[j],"\nPATTERN\t",cliques_pattern[j]
            '''
            pattern
            need to revise
            '''
            pattern=cliques_pattern[j]
            pattern[-1]=True
            pattern[-2]=True
            beds=[bed for bed in g.translate_bits_into_beds(cliques_pattern[j])]
            ret_dict["CLIQUES"][-1]["BEDS"]=[]
            for j1,bed in enumerate(beds):
                #logging.debug(i)
                score=score*1000.0/cdna_length
                chr=i.id
                id=i.id+"_"+"NO."+str(j0)+"_"+str(j1)
                itemRgb=str(rgb)+","+str(rgb)+","+str(rgb)
                bed=bed._replace(chr=chr,id=id,itemRgb=itemRgb)
                ret_dict["CLIQUES"][-1]["BEDS"].append(str(Tools.translate_coordinates(i,bed,True)))

           


            ret_dict["CLIQUES"][-1]["UNIQ_SCORE"]=uniq_score
            ret_dict["CLIQUES"][-1]["UNIQ_FPK"]=uniq_fpk
            ret_dict["CLIQUES"][-1]["SCORE"]=score

            #retv+="UNIQ\t"+str(uniq_score)+"\tBED\t"+str(Tools.translate_coordinates(Bed(i),bed,True))+"\n"
            #print >>sys.stderr,"UNIQ_FPK/CURRENT_MAX_UNIQ_FPK",uniq_fpk/max_uniq_fpk

            #retv+="\n"
    ret_dict["INTERPRET_FRG"]=float(cumulative_score)/total_frag
    return ret_dict

    



    


if __name__=="__main__":
    run()

