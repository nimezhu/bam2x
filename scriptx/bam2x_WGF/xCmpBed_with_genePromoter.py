#!/usr/bin/python
# Wei Guifeng, <guifengwei@gmail.com>
# -- coding: utf-8 --
'''
Usage: 
    To annotate the ChIPseq PEAKs to overlap which gene promoter and/or lncRNA promoter
    To annotate the PEAKs to overlap which intersted bed regions of yours
    
    python xCmpBed.py --gene refSeqGene.tab --RNA lincRNA.tab -f other.bed
    
    Copyright (c) 2012 - 2013, Wei Guifeng, <guifengwei@gmail.com>

'''
import sys, os, argparse, string
from zbed import *
from xplib.Annotation.Utils import *
from xplib import TableIO


def parse_argument():
    ''' argument parser'''
    p=argparse.ArgumentParser(description='example: %(prog)s -i *.bed --gene *.tab --RNA *.tab -o output',
                    epilog="dependency: python2.7, zbed")
    # *.bed
    p.add_argument('-i', '--input', dest='bed', metavar='interval', type=str, required=True, 
                    help = "the input of bed of interst, which will be annotated as intergenic,\
                            intronic, antisense, SenseOverlap with protein-coding gene ")
    # genefile
    p.add_argument('-G', '-g', '--gene', dest='genefile', metavar='Gene', type=str,
                    help="The reference annotation protein-coding gene file. Format: GenePred(genetab) ")
    # lncRNA file
    p.add_argument("-R", '-r', '--RNA', dest='rna', metavar='LncRNA', type=str, 
                    help='The reference annotation long noncoding RNAs file. Format: GenePred(genetab) ')
    # other featurefile
    p.add_argument('-F', '-f', '--feature', dest='feature', metavar='Feature', type=str,
                    help='The other feature of interest. Format: Bed ')
    # promoter size
    p.add_argument('-P', '--promoterSize', dest='bp', metavar='bp', type=int, default=2500,
                    help='The promoter region around the TSS. e.g. promoterSize=2500bp,\
                          in fact, the real promoter region is [-2500, 2500]')
    # output
    p.add_argument('-o', '--output', dest='output', type=str, default='output', help=" the Output file")
    if len(sys.argv) == 1 :
        sys.exit(p.print_help())
    args=p.parse_args()
    return args

def main():
    ''' main scripts '''
    args = parse_argument()
    
    PromoterList =[]
    feature = []

    print "# loading and reading the promoter of Gene ... "
    if args.genefile:
        for g in TableIO.parse(args.genefile, "genebed"):
            a = g.promoter(args.bp)
            PromoterList.append(a)
    if args.rna:
        for g in TableIO.parse(args.rna, "genebed"):
            a = g.promoter(args.bp)
            PromoterList.append(a)
    if args.feature:
        feature = TableIO.parse(args.feature, 'bed')
    print "# loading and reading Done !" 
    
    PromoterData = readIntoBinIndex(PromoterList)
    FeatureData  = readIntoBinIndex(feature)
    
    for i in TableIO.parse(args.bed, 'bed'):
        overlapGene, overlapFeature, overlap_string = [], [], ''
        if i.strand not in ['+','-']: continue
        else:
            OverlapGene = getOverlapFeatures(i, PromoterData)
            if FeatureData:
                overlapFeature = getOverlapFeatures(i, FeatureData)
            for g in overlapGene +  overlapFeature:
                overlap_string += g.id+';'
            print i, "\t", overlap_string


if __name__== "__main__":
    main()

