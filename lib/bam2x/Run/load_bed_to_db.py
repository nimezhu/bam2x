#!/usr/bin/env python
from __future__ import print_function
# Programmer : zhuxp
# Date: 
# Last-modified: 09-09-2014, 13:48:35 EDT
import os,sys,argparse
from bam2x import TableIO,Tools
from bam2x import IO
import sqlite3
import string
from bam2x.DBI.Templates import schema_templates
from bam2x.DBI.Templates import insert_templates

def set_parser(p):
    ''' This Function Parse the Argument '''
    p.add_argument("-D","--database",dest="db",type=str,help="database file name",default="guess")
    p.add_argument("-I","--input_format",dest="input_format",choices=schema_templates.keys(),type=str,help="input format , default: %(default)s",default="bed12")
    p.add_argument("-t","--table_name",dest="table_name",type=str,help="table name, default: %(default)s",default="test")
def help():
    return "read bed12 file into a sqlite3 db. if sqlite3 db is not exists , it will generate one"
    
def run(args):
    schema_template=schema_templates[args.input_format]
    SQL_template=insert_templates[args.input_format]
    db_filename=args.db
    out=IO.fopen(args.output,"w")
    if db_filename=="guess":
        db_filename=args.input.strip(".gz")+".db"
    db_is_new = not os.path.exists(db_filename)
    print("Database file : %s"%db_filename,file=out)
    with sqlite3.connect(db_filename) as conn:
        cursor=conn.cursor()
        if db_is_new:
            print ('Creating table %s if not exists\n________________________________'%args.table_name,file=out)
            S=schema_template.substitute({"table_name":args.table_name})
            print (S,file=out)
            print ("_______________________________",file=out)
            cursor.execute(S)
        else:
            S=schema_template.substitute({"table_name":args.table_name})
            cursor.execute(S)

        fin=IO.fopen(args.input,"r")
        S1=SQL_template.substitute({"table_name":args.table_name})
        print(S1,file=out)
        s=TableIO.parse(args.input,"simple")
        cursor.executemany(S1,s)
        conn.commit()
        print("loaded",file=out)

        
if __name__=="__main__":
    from bam2x.IO import parser_factory
    p=parser_factory(description=help())
    set_parser(p)
    if len(sys.argv)==1:
        print(p.print_help())
        exit(0)
    run(p.parse_args())
