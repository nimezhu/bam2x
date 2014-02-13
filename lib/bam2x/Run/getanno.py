import os
import sys
import logging
import argparse
from bam2x import TableIO,IO,DBI
from bam2x.Tools import seq_wrapper
def help():
    return "get [utr,exon,intron] annotations from bed12 file"
def set_parser(parser):
    parser.add_argument("-a","--anno",type=str,choices=("cds","utr5","utr3","exon","intron","utr"),dest="annotation",default="cdna")
    
def run(args):
    out=IO.fopen(args.output,"w")
    if args.annotation=="exon":
        for i in TableIO.parse(IO.fopen(args.input,"r"),"bed12"):
            for j in i.Exons():
                print >>out,j
    elif args.annotation=="intron":
        for i in TableIO.parse(IO.fopen(args.input,"r"),"bed12"):
            for j in i.Introns():
                print >>out,j
    elif args.annotation=="cds":
        for i in TableIO.parse(IO.fopen(args.input,"r"),"bed12"):
            j=i.cds()
            if j is not None and j.cdna_length() > 0:
                print >>out,j

    elif args.annotation=="cds":
        for i in TableIO.parse(IO.fopen(args.input,"r"),"bed12"):
            j=i.cds()
            if j is not None and j.cdna_length() > 0:
                print >>out,j

    elif args.annotation=="utr5":
        for i in TableIO.parse(IO.fopen(args.input,"r"),"bed12"):
            j=i.utr5()
            if j is not None and j.cdna_length() > 0:
                print >>out,j

    elif args.annotation=="utr3":
        for i in TableIO.parse(IO.fopen(args.input,"r"),"bed12"):
            j=i.utr3()
            if j is not None and j.cdna_length() > 0:
                print >>out,j
    elif args.annotation=="utr":
        for i in TableIO.parse(IO.fopen(args.input,"r"),"bed12"):
            j=i.utr5()
            if j is not None and j.cdna_length() > 0:
                print >>out,j
            j=i.utr3()
            if j is not None and j.cdna_length() > 0:
                print >>out,j




    








