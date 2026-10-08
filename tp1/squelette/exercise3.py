import sys
import utils

nuc2aa = {
    "TTT": "F", "TCT": "S", "TAT": "Y", "TGT": "C",
    "TTC": "F", "TCC": "S", "TAC": "Y", "TGC": "C",
    "TTA": "L", "TCA": "S", "TAA": "*", "TGA": "*",
    "TTG": "L", "TCG": "S", "TAG": "*", "TGG": "W",

    "CTT": "L", "CCT": "P", "CAT": "H", "CGT": "R",
    "CTC": "L", "CCC": "P", "CAC": "H", "CGC": "R",
    "CTA": "L", "CCA": "P", "CAA": "Q", "CGA": "R",
    "CTG": "L", "CCG": "P", "CAG": "Q", "CGG": "R",

    "ATT": "I", "ACT": "T", "AAT": "N", "AGT": "S",
    "ATC": "I", "ACC": "T", "AAC": "N", "AGC": "S",
    "ATA": "I", "ACA": "T", "AAA": "K", "AGA": "R",
    "ATG": "M", "ACG": "T", "AAG": "K", "AGG": "R",

    "GTT": "V", "GCT": "A", "GAT": "D", "GGT": "G",
    "GTC": "V", "GCC": "A", "GAC": "D", "GGC": "G",
    "GTA": "V", "GCA": "A", "GAA": "E", "GGA": "G",
    "GTG": "V", "GCG": "A", "GAG": "E", "GGG": "G"
}

def traduire_seq2aa_ORF(seq):

    # d'abord on trouve le plus grand cadre de lecture
    startindex = 0
    maxindex = 0
    maxlength = 0
    start = False
    for i in range(0, len(seq), 3):
        nucs = seq[i:i+3]
        if len(nucs) < 3: break
        aa = nuc2aa[nucs]
        if not start and aa == "M": 
            start = True
            startindex = i

        if start and aa == "*": 
            start = False
            orf_len = i - startindex
            if orf_len > maxlength: 
                maxlength = orf_len
                maxindex = startindex    

    aaseq = []
    start = False
    for i in range(maxindex, len(seq), 3):
        nucs = seq[i:i+3]
        if len(nucs) < 3: break

        aa = nuc2aa[nucs]
        if not start and aa == "M": start = True
        if start and aa == "*": break # le ribosome arrete sa lecture 
        if start : aaseq.append(aa)

        
    return "".join(aaseq)

def traduire_seq2aa_3ORF(seq):
    ORF1 = traduire_seq2aa_ORF(seq)
    ORF2 = traduire_seq2aa_ORF(seq[1:])
    ORF3 = traduire_seq2aa_ORF(seq[2:])
    return(ORF1,ORF2,ORF3)
    
protXnucseq = utils.read_single_fasta_sequence("../donnees/sequence.fasta")
# a) le cadre de lecture de la prot
print(traduire_seq2aa_3ORF(protXnucseq)[1])
print(utils.read_fasta_sequences("../donnees/geneX.fasta")["geneX"])