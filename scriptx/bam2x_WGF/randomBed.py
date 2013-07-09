#!/usr/bin/python
# Wei Guifeng, <guifengwei@gmail.com>
# -- coding: utf-8 --

import sys,random,argparse

_hg19_references=('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM')

_hg19_lengths=(249250621L, 243199373L, 198022430L, 191154276L, 180915260L, 171115067L, 159138663L, 146364022L, 141213431L, 135534747L, 135006516L, 133851895L, 115169878L, 107349540L, 102531392L, 90354753L, 81195210L, 78077248L, 59128983L, 63025520L, 48129895L, 51304566L, 155270560L, 59373566L, 16571L)


def parse_argument():
    """ argument parser"""
    p=argparse.ArgumentParser(description="example: randomBed.py --length 100 -n 500")
    p.add_argument("-l","--length",dest="length",metavar="length",type=int,required=True,
            help="length of the random bed region")
    p.add_argument("-n",dest="n",metavar="bed_number",default=100,type=int,
            help="how many beds you want? default=%(default)d")
    if len(sys.argv) == 1:
        sys.exit(p.print_help())
    else:
        args= p.parse_args()
        return args

def randomBed():
    ''' generate a random bed '''
    chrs=''
    start=0
    i=1
    orgnism = 'hg19'
    if orgnism == 'hg19':
        while i< args.n+1:
            chrs = random.choice(_hg19_references)
            start= random.randint(0, _hg19_lengths[ _hg19_references.index(chrs) ] - 1 - args.length )
            print "%s\t%d\t%d\t%s\t%d\t%s" %(chrs, start, start+args.length, "randomBed_"+str(i),0,"." )
            i+=1
    else:
        print >>sys.stderr,"genome size worng"

def main():
    """ main scripts """
    global args
    args = parse_argument()
    randomBed()


if __name__=="__main__":
    main()

