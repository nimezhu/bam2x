#!/usr/bin/env python
from __future__ import print_function
# Programmer : zhuxp
# Date: 
# Last-modified: 02-25-2014, 16:29:20 EST
import os,sys,argparse
from bam2x import TableIO,Tools
from bam2x import IO
from bam2x.Annotation import BED12 as Bed12
import sqlite3
import string

template=string.Template("""select * from $table_name where name=\"$name\"""")

def set_parser(p):
    ''' This Function Parse the Argument '''
    p.add_argument("-D","--database",dest="db",type=str,help="database file name",default="guess")
    p.add_argument("-t","--table_name",dest="table_name",type=str,help="table name",default="test")
def help():
    return "query sqlite3 database to get a transcipt, and then do further analysis. this is a demo program."
def factory(cursor,r):
    return Bed12._make(Bed12._types(r[1:]))
def run(args):
    db_filename=args.db
    out=IO.fopen(args.output,"w")
    if os.path.exists(args.input):
        fin=IO.fopen(args.input,"r")
    else:
        fin=(args.input,)
    if not os.path.exists(db_filename):
        print("can't find database %s"%db_filename)
        exit(1)
    print("Database file : %s"%db_filename,file=out)
    with sqlite3.connect(db_filename) as conn:
        conn.row_factory=factory
        cursor=conn.cursor()
        for i in fin:
            i=i.strip()
            i=i.strip(" ")
            print("query %s"%i,file=out)
            s=template.substitute({"table_name":args.table_name,"name":i})
            print(s,file=out)
            cursor.execute(s)
            r=None
            try:
                r=cursor.fetchone()
                print(r,file=out)
            except:
                raise
        
if __name__=="__main__":
    from bam2x.IO import parser_factory
    p=parser_factory(description=help())
    set_parser(p)
    if len(sys.argv)==1:
        print(p.print_help())
        exit(0)
    run(p.parse_args())
