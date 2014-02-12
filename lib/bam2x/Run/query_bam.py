import os
import sys
import logging
import argparse
def help():
    return "query_bam help"
def set_parser(parser):
    parser.add_argument("-m",type=str,choices=("pileup","fetch"),dest="method")
    
    
def run(args):
    print "in query_bam module now"
    print args







