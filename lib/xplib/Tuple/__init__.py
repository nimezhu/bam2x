

class Parser:
    def __init__(self,header):
        self.header=header
    def get(self,x,attr):
        return x[self.header[attr]]
    def set(self,x,attr,value):
        '''
        tuple can't change in function
        usage:
        x=set(x,"chr","chr01")
        '''
        a=list(x)
        a[self.header[attr]]=value
        return tuple(a)




import BedTuple
def test():
    a=("chr1",1,2)
    b=("chr1",2,3)
    c=("chr2",2,3)
    l=[a,c,b]
    print l
    BedTuple.tuple_sort(l)
    print l
    parser=Parser(BedTuple.tuple_header)
    for i in l:
        print parser.get(i,"chr")
        print BedTuple.tuple_factory(i)
if __name__=="__main__":
    test()



