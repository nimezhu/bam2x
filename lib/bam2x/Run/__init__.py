__all__=["getseq","query_bam","read","sort","cmpgene","getanno","query_RNASeq","load_bed_to_db","query_db","query_RNA_RPKM","tabix"]
from . import *

'''
command lines
'''
commands= {
    "getseq":getseq,
    "getanno":getanno,
    "query_bam":query_bam,
    "query_RNASeq":query_RNASeq,
    "query_RNA_RPKM":query_RNA_RPKM,
    "read":read,
    "sort":sort,
    "cmpgene":cmpgene,
    "load_bed_to_db":load_bed_to_db,
    "query_db":query_db,
    "tabix":tabix
}


