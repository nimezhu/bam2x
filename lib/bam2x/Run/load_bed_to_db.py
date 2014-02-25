#!/usr/bin/env python
from __future__ import print_function
# Programmer : zhuxp
# Date: 
# Last-modified: 02-25-2014, 15:00:54 EST
import os,sys,argparse
from bam2x import TableIO,Tools
from bam2x import IO
import sqlite3
import string
'''
table schema template
'''
schema_template= string.Template("""create table if not exists $table_name (
    id      integer primary key autoincrement not null,
    chr     varchar(20),
    start   integer,
    stop    integer,
    name    varchar(200),
    score   float,
    strand  character(1),
    cds_start   integer,
    cds_stop    integer,
    itemRgb     varchar(20),
    blockCount  integer,
    blockSizes  varchar,
    blockStarts varchar
)
""")

SQL_template=string.Template("""insert into $table (chr,start,stop,name,score,strand,cds_start,cds_stop,itemRgb,blockCount,blockSizes,blockStarts)
                            values (?,?,?,?,?,?,?,?,?,?,?,?)
                            """)
   

def set_parser(p):
    ''' This Function Parse the Argument '''
    p.add_argument("-D","--database",dest="db",type=str,help="database file name",default="test.db")
    p.add_argument("-t","--table_name",dest="table_name",type=str,help="table name",default="test")
def help():
    return "read bed12 file and generate a sqlite3 db."
    
def run(args):
    db_filename=args.db
    db_is_new = not os.path.exists(db_filename)
    out=IO.fopen(args.output,"w")
    with sqlite3.connect(db_filename) as conn:
        cursor=conn.cursor()
        if db_is_new:
            print ('Creating schema',file=out)
            S=schema_template.substitute({"table_name":args.table_name})
            print (S,file=out)
            cursor.execute(S)
        fin=IO.fopen(args.input,"r")
        S1=SQL_template.substitute({"table":args.table_name})
        print(S1,file=out)
        s=TableIO.parse(args.input,"simple")
        for i in s:
            cursor.execute(S1,i)
        conn.commit()

        
if __name__=="__main__":
    from bam2x.IO import parser_factory
    p=parser_factory(description=help())
    set_parser(p)
    if len(sys.argv)==1:
        print(p.print_help())
        exit(0)
    run(p.parse_args())
