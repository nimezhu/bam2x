from __future__ import print_function
import os
import sys
import logging
import argparse
from  bam2x import IO,TableIO,DBI
from bam2x.Tools import compatible,merge_beds
import re
def help():
    return "group bed12 into clusters and groups."
def set_parser(parser):
    pass
    #parser.add_argument("-m",type=str,choices=("seq","cDNA","cdna","cds","utr5","utr3"),dest="method")
    
    
def run(args):
    fin=IO.fopen(args.input,"r")
    out=IO.fopen(args.output,"w")
    beds=[i for i in TableIO.parse(fin,"bed12")]
    beds.sort()
    for i,x in enumerate(iter_cluster(beds)):
        id=find_prefix_consensus([i0.id for i0 in x[1]])
        strand=find_consensus_strand([i0.strand for i0 in x[1]])
        print("REGION\tCL_{index}\t{chr}\t{start}\t{end}\t{id}\t{score}\t{strand}".format(strand=strand,score=len(x[1]),chr=x[1][0].chr,start=x[1][0].start,end=x[0],index=str(i+1),id=id),file=out)
        
        for j,y in enumerate(greedy_iter_compatible_group(x[1])):
            print("\tGROUP{j}\t{bed}".format(j=j+1,bed=merge_beds(y,id="CL.{i}_GP.{j}".format(i=i+1,j=j+1))),file=out)
            for k,z in enumerate(sorted(y,key= lambda x0:x0.cdna_length(), reverse=True)):
                print("\t\tCL.{i}_GP.{j}_TR.{k}\t{l}\t{z}".format(i=i+1,j=j+1,k=k+1,l=z.cdna_length(),z=z),file=out)
def find_consensus_strand(l):
    strand=l[0]
    for i in l[1:]:
        if i!=strand: return "."
    return strand


def find_prefix_consensus(l):
    h={}
    sep=re.compile("[_\.]")
    for i in l:
        x=sep.split(i)
        if h.has_key(x[0]): 
            h[x[0]]+=1
        else:
            h[x[0]]=1
    a=sorted(h.keys())
    return "_".join(a)
def greedy_iter_compatible_group(beds):
    local_beds=[i for i in beds]
    compatible_groups=[[]]
    for i in local_beds:
        sign=False
        for j,c in enumerate(compatible_groups):
            if compatible_with_group(i,c):
                c.append(i)
                sign=True
                break
        if not sign:
            compatible_groups.append([i])

    for i in sorted(compatible_groups):
        yield i








def compatible_with_group(i,group):
    for j in group:
        if not compatible(i,j):
            return False
    return True




def iter_cluster(beds):
    cluster=[]
    last=beds[0]
    for i in beds:
        if i.chr==last.chr and i.start < last.stop:
            cluster.append(i)
            if i.stop > last.stop:
                last = i
        else:
            if len(cluster) > 0 : yield last.stop,cluster
            cluster=[i]
            last=i
    if len(cluster) > 0 : yield last.stop,cluster


if __name__=="__main__":
    from bam2x.IO import parser_factory
    p=parser_factory(description=help())
    set_parser(p)
    if len(sys.argv)==1:
        print(p.print_help())
        exit(0)
    run(p.parse_args())







