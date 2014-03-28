import string
from bam2x.Annotation import BED12,BED6,BED3
import itertools

hTypesToSqlite3 = {
    "int":"INTEGER",
    "long":"INTEGER",
    "float":"REAL",
    "double":"REAL",
    "text":"TEXT",
    "blob":"NONE",
    "boolean":"NUMERIC",
    "date":"NUMERIC"
}
def translate_to_sqlite3_format(f):
    if hTypesToSqlite3.has_key(f):
        return hTypesToSqlite3[f]
    else:
        return "TEXT"


def sqlite3_template_factory(ids,formats):
    schema="create table if not exists $table_name (\n"
    id_str=""
    v=""
    for i,f in itertools.izip(ids,formats):
        schema+="{id} {format},\n".format(id=i,format=translate_to_sqlite3_format(f.__name__))
        if i=="id" or i=="name":
            schema=schema.strip("\n")
            schema=schema.strip(",")
            schema+=" primary key,\n"
        id_str+=i+","
        v+="?,"
    id_str=id_str.strip(",")
    v=v.strip(",")
    schema=schema.strip("\n")
    schema=schema.strip(",")
    schema+=")"
    insert="insert into $table_name ({id_str}) values ({v})".format(id_str=id_str,v=v)
    return string.Template(schema),string.Template(insert)


bed12_schema_template= string.Template("""create table if not exists $table_name (
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

bed12_insert_template=string.Template("""insert into $table_name (chr,start,stop,name,score,strand,cds_start,cds_stop,itemRgb,blockCount,blockSizes,blockStarts)
                            values (?,?,?,?,?,?,?,?,?,?,?,?)
                            """)
 
bed12_select_template=string.Template("""select * from $table_name where name==$name""")


bed6_schema_template= string.Template("""create table if not exists $table_name (
    id      integer primary key autoincrement not null,
    chr     varchar(20),
    start   integer,
    stop    integer,
    name    varchar(200),
    score   float,
    strand  character(1)
)
""")

bed3_schema_template= string.Template("""create table if not exists $table_name (
    id      integer primary key autoincrement not null,
    chr     varchar(20),
    start   integer,
    stop    integer
)
""")

bed6_insert_template=string.Template("""insert into $table_name (chr,start,stop,name,score,strand)
                            values (?,?,?,?,?,?)
                            """)
 

bed3_insert_template=string.Template("""insert into $table_name (chr,start,stop)
                            values (?,?,?)
                            """)
 

schema_templates={
    "bed12":bed12_schema_template,
    "bed6":bed6_schema_template,
    "bed3":bed3_schema_template,
}
insert_templates={   
    "bed12":bed12_insert_template,
    "bed6":bed6_insert_template,
    "bed3":bed3_insert_template,
}

factories = {
"bed12":lambda c,r:BED12._make(BED12._types(r[1:])),
"bed6":lambda c,r:BED6._make(BED6._types(r[1:])),
"bed3":lambda c,r:BED3._make(BED3._types(r[1:])),
}


select_template=string.Template("""select * from $table_name where name=\"$name\"""")



