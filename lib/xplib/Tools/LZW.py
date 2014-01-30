def compress(seq):
    w=""
    d={"A":"A","C":"C","G":"G","T":"T"}
    i=0
    S=""
    for k in seq:
        if d.has_key(w+k):
            w=w+k
        else:
            d[w+k]=i
            i+=1
            S+=str(d[w])+","
            w=k
    for i in d.keys():
        print i,d[i]
    return S


def complexity(seq):
    w=""
    d={"A":"A","C":"C","G":"G","T":"T"}
    i=0
    S=0
    for k in seq.upper():
        if d.has_key(w+k):
            w=w+k
        else:
            d[w+k]=i
            i+=1
            S+=1
            w=k
    return S


if __name__=="__main__":
    b="ACGATACGGAATT"
    a="ACACACACACACACACACACACACACACACACACAC"
    c="ACGACGACGACGA"
    print a,len(a),complexity(a)
    print b,len(b),complexity(b)
    print c,len(c),complexity(c)



