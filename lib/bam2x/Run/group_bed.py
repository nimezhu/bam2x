from __future__ import print_function
import os
import sys
import logging
import argparse
from  bam2x import IO,TableIO,DBI
from bam2x.Tools import compatible
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
        print("CL_"+str(i+1),file=out)
        for j,y in enumerate(greedy_iter_compatible_group(x)):
            print("\tGROUP{j}".format(j=j+1),file=out)
            for k,z in enumerate(sorted(y,key= lambda x0:x0.cdna_length(), reverse=True)):
                print("\tCL.{i}_GP.{j}_TR.{k}\t{l}\t{z}".format(i=i+1,j=j+1,k=k+1,l=z.cdna_length(),z=z),file=out)


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
            if len(cluster) > 0 : yield cluster
            cluster=[i]
            last=i
    if len(cluster) > 0 : yield cluster


if __name__=="__main__":
    from bam2x.IO import parser_factory
    p=parser_factory(description=help())
    set_parser(p)
    if len(sys.argv)==1:
        print(p.print_help())
        exit(0)
    run(p.parse_args())







