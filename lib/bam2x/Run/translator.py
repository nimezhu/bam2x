from __future__ import print_function
import os
import sys
import logging
import argparse
import sqlite3
from os.path import splitext
import string
from bam2x.Tools import reverse_translate
from bam2x.DBI.Templates import bed12_schema_template as schema_t
from bam2x.DBI.Templates import bed12_insert_template as insert_t
import logging
from bam2x import TableIO,IO
from bam2x.Annotation import BED12 as Bed12
def help():
    return "translate bed from gene coordinates to chromosome coordinates"
def set_parser(parser):
    parser.add_argument("-t","--translator",type=str,dest="translator",help="translator bed file or db[sqlite3] (gene coordinates in chromosome)")
    parser.add_argument("-T","--table_name",type=str,dest="table_name",default="test")
    
def _generate_db(filename,db_filename,table_name):
    with sqlite3.connect(db_filename) as conn:
        cursor=conn.cursor()
        S=schema_t.substitute({"table_name":table_name})
        cursor.execute(S)
        LOAD_S=insert_t.substitute({"table_name":table_name})
        s=TableIO.parse(IO.fopen(filename,"r"),"simple")
        cursor.executemany(LOAD_S,s)
        conn.commit()

template=string.Template("""select * from $table_name where name='$name'""")
def run(args):
    logging.basicConfig(level=logging.INFO)
    db_filename=args.translator
    t_name,t_ext=splitext(args.translator)
    '''
    test if it is db file
    generate db file if it doesn't exists.
    '''
    if t_ext!=".db":
        #possible_db=args.translator.strip("\\.gz")+".db"
        possible_db=args.translator+".db"
        if os.path.exists(possible_db):
            db_filename=possible_db
        else:
            _generate_db(args.translator,possible_db,args.table_name)
            db_filename=possible_db
    
    '''
    query db file
    '''
    out=IO.fopen(args.output,"w")
    with sqlite3.connect(db_filename) as conn:
        conn.row_factory=lambda conn,x: Bed12._make(Bed12._types(x[1:]))
        cursor=conn.cursor()
        for i in TableIO.parse(IO.fopen(args.input,"r"),"bed"):
            s=template.substitute({"table_name":args.table_name,"name":i.chr.strip()})
            #print(s)
            cursor.execute(s)
            gene=None
            try:
                gene=cursor.fetchone()
                logging.debug(i)
                logging.debug(i.cdna_length())
                logging.debug(gene)
                logging.debug(gene.cdna_length())
            except:
                raise
                logging.warning("can't find gene %s"%i.chr)
                continue
            assert gene.cdna_length() > i.cdna_length()
            print(reverse_translate(gene,i),file=out) 
    

if __name__=="__main__":
    from bam2x.IO import parser_factory
    p=parser_factory(description=help())
    set_parser(p)
    if len(sys.argv)==1:
        print(p.print_help())
        exit(0)
    run(p.parse_args())







