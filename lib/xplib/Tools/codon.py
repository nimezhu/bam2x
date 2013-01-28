codon2acid={
"GCT":"A",
"GCC":"A",
"GCA":"A",
"GCG":"A",
"TTA":"L",
"TTG":"L",
"CTT":"L",
"CTC":"L",
"CTA":"L",
"CTG":"L",
"CGT":"R",
"CGC":"R",
"CGA":"R",
"CGG":"R",
"AGA":"R",
"AGG":"R",
"AAA":"K",
"AAG":"K",
"AAT":"N",
"AAC":"N",
"ATG":"M",
"GAT":"D",
"GAC":"D",
"TTT":"F",
"TTC":"F",
"TGT":"C",
"TGC":"C",
"CCT":"P",
"CCC":"P",
"CCA":"P",
"CCG":"P",
"CAA":"Q",
"CAG":"Q",
"TCT":"S",
"TCC":"S",
"TCA":"S",
"TCG":"S",
"AGT":"S",
"AGC":"S",
"GAA":"E",
"GAG":"E",
"ACT":"T",
"ACC":"T",
"ACA":"T",
"ACG":"T",
"GGT":"G",
"GGC":"G",
"GGA":"G",
"GGG":"G",
"TGG":"W",
"CAT":"H",
"CAC":"H",
"TAT":"Y",
"TAC":"Y",
"ATT":"I",
"ATC":"I",
"ATA":"I",
"GTT":"V",
"GTC":"V",
"GTA":"V",
"GTG":"V",
"TAG":")",
"TGA":")",
"TAA":")"
}

acid2codon={
"A":("GCT","GCC","GCA","GCG"),
"E":("GAA","GAG"),
"S":("TCT","TCC","TCA","TCG","AGT","AGC"),
"L":("TTA","TTG","CTT","CTC","CTA","CTG"),
"C":("TGT","TGC"),
"P":("CCT","CCC","CCA","CCG"),
"F":("TTT","TTC"),
"G":("GGT","GGC","GGA","GGG"),
")":("TAG","TGA","TAA"),
"H":("CAT","CAC"),
"V":("GTT","GTC","GTA","GTG"),
"K":("AAA","AAG"),
"N":("AAT","AAC"),
"T":("ACT","ACC","ACA","ACG"),
"Q":("CAA","CAG"),
"I":("ATT","ATC","ATA"),
"W":("TGG"),
"D":("GAT","GAC"),
"M":("ATG"),
"R":("CGT","CGC","CGA","CGG","AGA","AGG"),
"Y":("TAT","TAC"),
"(":("ATG")
}
acid2codon_number={
"A":4,
"E":2,
"S":6,
"L":6,
"C":2,
"P":4,
"F":2,
"G":4,
")":3,
"H":2,
"V":4,
"K":2,
"N":2,
"T":4,
"Q":2,
"I":3,
"W":1,
"D":2,
"M":1,
"R":6,
"Y":2,
"(":1,
}


def get_syn_codon_seq(seq):
    from random import random
    seq=seq.upper()
    new_codons=[]
    for i in range(len(seq)/3):
        codon=seq[i*3:i*3+3]
        if codon2acid.has_key(codon):
            acid=codon2acid[codon]
            c=int(random()*120000)
            c=c%acid2codon_number[acid]
            new_codon=acid2codon[acid][c]
        else:
            new_codon=codon
        new_codons.append(new_codon)
    return "".join(new_codons)

def get_seq_codon_bit(seq):
    from math import log
    seq=seq.upper()
    bit=0.0
    for i in range(len(seq)/3):
        codon=seq[i*3:i*3+3]
        if codon2acid.has_key(codon):
            c=acid2codon_number[codon2acid[codon]]
            bit+=log(1/float(c))/log(2)
    return -bit
def get_bitseq(seq):
    seq=seq.upper()
    bitseq=""
    for i in range(len(seq)/3):
        codon=seq[i*3:i*3+3]
        if codon2acid.has_key(codon):
            c=acid2codon_number[codon2acid[codon]]
            bitseq+=str(c)
    return bitseq




def pro_seq_codon_bit(seq):
    from math import log
    seq=seq.upper()
    bit=0.0
    for i in range(len(seq)):
        acid=seq[i]
        if acid2codon_number.has_key(acid):
            c=acid2codon_number[acid]
            bit+=log(1/float(c))/log(2)
    return -bit
