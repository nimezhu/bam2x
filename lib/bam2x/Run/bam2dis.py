from __future__ import print_function
import os
import sys
import logging
import argparse
from  bam2x import IO,TableIO,DBI
import itertools
from bam2x.Tools import compatible,extend_slice,compatible_with_transcript,_translate,overlap
import logging
'''
Step 1. Group Bed12 into Gene Models(Cluster Isoforms)
Step 2. Compare reads with Gene Models.
'''
def help():
    return "compare bam with bed12. input is bed12 file. aggregation plot nearby tss, nearby tts , and gene body"
def set_parser(parser):
    pass
    #parser.add_argument("-m",type=str,choices=("seq","cDNA","cdna","cds","utr5","utr3"),dest="method")
    parser.add_argument("-b",type=str,dest="bam",help="bam file");
    parser.add_argument("-s","--strand",type=str,dest="strand",choices=("read1","read2"),default="read1",help="strand default:%(default)s");
    parser.add_argument("--bp",type=int,dest="bp",default=2000,help="upstream/downstream bp default:%(default)s");
    #parser.add_argument("-R","--RNASeq",action='store_true',dest="RNASeq",help="if it is RNASeq");
def getNM(bed):
    return int(bed.itemRgb.split(",")[0])
def run(args):
    logging.basicConfig(level=logging.DEBUG) 
    fin=IO.fopen(args.input,"r")
    out=IO.fopen(args.output,"w")
    bam=DBI.init(args.bam,"bam");
    beds=[i for i in TableIO.parse(fin,"bed12")]
    beds.sort()
    bp=args.bp
    print("mapped:{}".format(bam.mapped))
    print("unmapped:{}".format(bam.unmapped))
    data={}
    for i,x in enumerate(iter_cluster(beds)):
        print("{}\t{}:{}-{}".format(i+1,x["chr"],x["start"]+1,x["stop"]))
        '''
        cds=[z.cds() for z in x["beds"] if z.cds()]
        utr3=[z.utr3() for z in x["beds"] if z.utr3()]
        utr5=[z.utr5() for z in x["beds"] if z.utr5()]
        '''
        
        coords = [ up_down_coordinate(gene,args.bp,args.bp) for gene in x["beds"] ]
        for j,y in enumerate(coords):
            data[y.id]={}
            data[y.id]["coord"]=y
            data[y.id]["values"]=[0.0 for l in range(y.cdna_length())];
        coord_beds = [ _translate(coord,bed) for coord,bed in itertools.izip(coords,x["beds"])]
        for j,read in enumerate(bam.query(method="bam1",chr=x["chr"],start=x["start"]-args.bp,stop=x["stop"]+args.bp,strand=args.strand)):
            NM=getNM(read)  # number of hits
            NC=0            # number of compatible 
            c_coords=[]
            for k,coord in enumerate(coords):
                if overlap(read,coord) and compatible(read,coord): # don't consider the reads extend out of coords.
                    NC+=1
                    c_coords.append(k)
            for k,c in enumerate(c_coords):
                coord=coords[c]
                if read.start < coord.start or read.stop > coord.stop:
                    start=max(read.start,coord.start)
                    stop=min(read.stop,coord.stop)
                    read=read._slice(start,stop)
                read_in_coord = _translate(coord,read)
                for l in xrange(read_in_coord.start,read_in_coord.stop):
                    data[coord.id]["values"][l]+=1.0/NC/NM
        for j,y in enumerate(coords):
            print(data[y.id]["coord"])
            print(data[y.id]["values"])

    '''
    for i in data.keys():
        print(data[i]["coord"])
        print(data[i]["values"])
    ''' 


                    
            




def up_down_coordinate(gene,up=1000,down=1000,chromSizes={}):
    new_start=gene.start-up
    if new_start<0: 
        new_start=0
        if gene.strand=="+" or gene.strand==".":
            up=start
        elif gene.strand=="-":
            down=start
    new_stop=gene.stop+down
    if chromSizes.has_key(gene.chr):
        chromSize=chromSizes[gene.chr]
        if new_stop>chromSize:
            new_stop=chromSize
            if strand=="+" or strand==".":
                down=new_stop-gene.stop
            elif strand=="-":
                up=new_stop-gene.stop
    coord=extend_slice(gene,new_start,new_stop)
    coord=coord._replace(id=coord.id+"_up"+str(up)+"_down"+str(down))
    return coord


def iter_cluster(beds):
    retv={}
    cluster=[]
    last=beds[0]
    retv["start"]=last.start
    retv["chr"]=last.chr
    for i in beds:
        if i.chr==last.chr and i.start < last.stop:
            cluster.append(i)
            if i.stop > last.stop:
                last = i

        else:
            if len(cluster) > 0 :
                retv["stop"]=last.stop
                retv["beds"]=cluster
                yield retv
            cluster=[i]
            last=i
            retv={}
            retv["start"]=last.start
            retv["chr"]=last.chr
    if len(cluster) > 0 : 
        retv["stop"]=last.stop
        retv["beds"]=cluster
        yield retv 


if __name__=="__main__":
    from bam2x.IO import parser_factory
    p=parser_factory(description=help())
    set_parser(p)
    if len(sys.argv)==1:
        print(p.print_help())
        exit(0)
    run(p.parse_args())


