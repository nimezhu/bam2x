import struct
import array
import collections
from bam2x.Turing import  TuringCode
import TuringCodeBook as cb
def write_int_array(ints,f):
    sign=0
    if type(f)==type("str"):
        try:
            f=open(f,"wb")
            sign=1
        except:
            print >>sys.stderr,"can't open file to write"
            exit(0)
    for i in ints:
        f.write(struct.pack('i',i))
    if sign==1:
        f.close()


def read_int_array(f):
    '''
    read any iterable int array or read from file
    '''
    sign=0
    if type(f)==type("str"):
        try:
            f=open(f,"rb")
            sign=1
        except:
            print >>sys.stderr,"can't open file to read"
            exit(0)
    if isinstance(f,file):
        while True:
            a = array.array('i')
            a.fromstring(f.read(4096))
            if not a:
                break
            for i in a:
                yield i
        if sign==1:
            f.close()
        raise StopIteration
    elif isinstance(f,collections.Iterable):
        for i in f:
            yield i
        raise StopIteration

def write_turing_array(turings,f,mode="n"):
    '''
    HEADERCODE: normal 3721
    HEADERCODE: compress 9981
    '''
    if mode=="n":
        header=cb.NORMAL_CODES_ARRAY_HEADER
    elif mode=="c":
        header=cb.COMPRESS_CODES_ARRAY_HEADER

    sign=0
    if type(f)==type("str"):
        try:
            f=open(f,"wb")
            sign=1
        except:
            print >>sys.stderr,"can't open file to write"
            exit(0)
    f.write(struct.pack("i",header))
    if mode=="n":
        for i in turings:
            f.write(struct.pack('i',i.pos))
            f.write(struct.pack('i',i.code))
    elif mode=="c":
        last=turings.next()
        number=1
        for i in turings:
            if i.pos==last.pos and i.code==last.code:
                number+=1
            else:
                f.write(struct.pack('i',last.pos))
                f.write(struct.pack('i',last.code))
                f.write(struct.pack('i',number))
                last=i
                number=1
        f.write(struct.pack('i',last.pos))
        f.write(struct.pack('i',last.code))
        f.write(struct.pack('i',number))
    if sign==1:
        f.close()



def read_turing_array(f,mode="n"):
    '''
    read any turing encode array 
    turing encode array:
      sorting array and has a header
    parse it into NORMAL MODE or COMPRESS MODE
    '''
    r=read_int_array(f)
    header=r.next()
    if mode=="n":
        if header==cb.NORMAL_CODES_ARRAY_HEADER:
            while True:
                try:
                    pos=r.next()
                    code=r.next()
                    yield TuringCode(pos,code)
                except StopIteration:
                    raise StopIteration
        elif header==cb.COMPRESS_CODES_ARRAY_HEADER:
            while True:
                try:
                    pos=r.next()
                    code=r.next()
                    number=r.next()
                    for i in range(number):
                        yield TuringCode(pos,code)
                except StopIteration:
                    raise StopIteration
        raise StopIteration
    elif mode=="c":
        last_pos=None
        last_code=None
        number=None
        if header==cb.COMPRESS_CODES_ARRAY_HEADER:
            while True:
                try:
                    pos=r.next()
                    code=r.next()
                    number=r.next()
                    yield TuringCode(pos,code),number
                except:
                    raise StopIteration
        elif header==cb.NORMAL_CODES_ARRAY_HEADER:
            try:
                last_pos=r.next()
                last_code=r.next()
                number=1
            except StopIteration:
                raise StopIteration
            while True:
                try:
                    pos=r.next()
                    code=r.next()
                    if pos==last_pos and code==last_code:
                        number+=1
                    else:
                        yield TuringCode(last_pos,last_code),number
                        number=1
                        last_pos=pos
                        last_code=code
                except StopIteration:
                    break
            yield TuringCode(last_pos,last_code),number
            raise StopIteration
        raise StopIteration



