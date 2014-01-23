from xplib.Annotation import Bed
tuple_header={"chr":0,"start":1,"stop":2,"id":3,"score":4,"strand":5}

def tuple_cmp(a,b):
    '''
    cmp two bed tuples
    '''
    return cmp(a[0],b[0]) or cmp(a[1],b[1]) or cmp(a[2],b[2])
def tuple_factory(a):
    return Bed(a)
def tuple_sort(l):
    l.sort(cmp=tuple_cmp)
def tuple_len(a):
    return a[2]-a[1]
'''
below is function for bed
'''
def get_exons(a):
    '''
    return bed tuple or a list of bed tuple?
    or an iterator?
    '''
    #TODO
    pass
def get_introns(a):
    '''
    same as get_exons
    '''
    #TODO
    pass

def is_overlap(a,b):
    '''
    '''
    #TODO
    pass
def get_upstream(a):
    '''
    '''
    #TODO
    pass
def get_downstream(a):
    '''
    '''
    #TODO
    pass
