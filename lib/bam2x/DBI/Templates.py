import string
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



