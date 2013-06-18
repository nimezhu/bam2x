#!/usr/bin/python
# Wei Guifeng, <guifengwei@gmail.com>
# -- coding: utf-8 --

import sys, os, argparse, string
from zbed import *
from xplib.Annotation.Utils import *
from xplib import TableIO

def parse_argument():
    ''' argument parser'''
    p=argparse.ArgumentParser(description='example: %(prog)s -i *.bed --gene *.tab -o output',
                    epilog="dependency: python2.7, zbed")
    # *.bed
    p.add_argument('-i', '--input', dest='bed', metavar='interval', type=str, required=True, 
                    help = "the input of bed of interst, which will be annotated as intergenic,\
                            intronic, antisense, SenseOverlap with protein-coding gene ")
    # genefile
    p.add_argument('-G', '-g', '--gene', dest='genefile', metavar='Gene', type=str,
                    help="The reference annotation protein-coding gene file. Format: GenePred(genetab) ")
    # output
    p.add_argument('-o', '--output', dest='output', type=str, default='output', help=" the Output file")
    if len(sys.argv) == 1 :
        sys.exit(p.print_help())
    args=p.parse_args()
    return args

def main():
    ''' main scripts '''
    args = parse_argument()
    bed  = args.bed
    gene = readIntoBinIndex(TableIO.parse( args.genefile, "genebed") )
    for i in TableIO.parse(args.bed, 'bed'):
        if i.strand not in ['+','-']: continue
        else:
            OverlapGene    = getOverlapFeatures(i, gene)
            Overlap_dict   = Classify_Overlap(i, OverlapGene)
            overlap_string = ''
            for k,v in Overlap_dict.iteritems():
                if v:
                    overlap_string += "".join([ str(k+'_'+each)+';' for each in v])
            if not overlap_string:
                overlap_string = 'intergenic'
            print i, overlap_string

if __name__== "__main__":
    main()

