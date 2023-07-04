codon_freq_table = {
    "UUU": 5.0, "UUC": 27.1, "UUA": 0.6, "UUG": 4.0, 
    "UCU": 4.7, "UCC": 16.1, "UCA": 3.2, "UCG": 16.1, 
    "UAU": 2.6, "UAC": 22.8, "UAA": 1.0, "UAG": 0.4, 
    "UGU": 1.4, "UGC": 13.1, "UGA": 0.5, "UGG": 13.2, 
    "CUU": 4.4, "CUC": 13.0, "CUA": 2.6, "CUG": 65.2, 
    "CCU": 8.1, "CCC": 29.5, "CCA": 5.1, "CCG": 20.7, 
    "CAU": 2.2, "CAC": 17.2, "CAA": 4.2, "CAG": 36.3, 
    "CGU": 4.9, "CGC": 34.9, "CGA": 2.0, "CGG": 11.2, 
    "AUU": 8.0, "AUC": 26.6, "AUA": 1.1, "AUG": 25.7, 
    "ACU": 5.2, "ACC": 27.7, "ACA": 4.1, "ACG": 15.9, 
    "AAU": 2.8, "AAC": 28.5, "AAA": 2.4, "AAG": 43.3, 
    "AGU": 2.6, "AGC": 22.8, "AGA": 0.7, "AGG": 2.7, 
    "GUU": 5.1, "GUC": 15.4, "GUA": 2.0, "GUG": 46.5, 
    "GCU": 16.7, "GCC": 54.6, "GCA": 10.6, "GCG": 44.4, 
    "GAU": 6.7, "GAC": 41.7, "GAA": 2.8, "GAG": 53.5, 
    "GGU": 9.5, "GGC": 62.0, "GGA": 5.0, "GGG": 9.7
}

codon_freq_table = {
    "TTT": 5.0, "TTC": 27.1, "TTA": 0.6, "TTG": 4.0, 
    "TCT": 4.7, "TCC": 16.1, "TCA": 3.2, "TCG": 16.1, 
    "TAT": 2.6, "TAC": 22.8, "TAA": 1.0, "TAG": 0.4, 
    "TGT": 1.4, "TGC": 13.1, "TGA": 0.5, "TGG": 13.2, 
    "CTT": 4.4, "CTC": 13.0, "CTA": 2.6, "CTG": 65.2, 
    "CCT": 8.1, "CCC": 29.5, "CCA": 5.1, "CCG": 20.7, 
    "CAT": 2.2, "CAC": 17.2, "CAA": 4.2, "CAG": 36.3, 
    "CGT": 4.9, "CGC": 34.9, "CGA": 2.0, "CGG": 11.2, 
    "ATT": 8.0, "ATC": 26.6, "ATA": 1.1, "ATG": 25.7, 
    "ACT": 5.2, "ACC": 27.7, "ACA": 4.1, "ACG": 15.9, 
    "AAT": 2.8, "AAC": 28.5, "AAA": 2.4, "AAG": 43.3, 
    "AGT": 2.6, "AGC": 22.8, "AGA": 0.7, "AGG": 2.7, 
    "GTT": 5.1, "GTC": 15.4, "GTA": 2.0, "GTG": 46.5, 
    "GCT": 16.7, "GCC": 54.6, "GCA": 10.6, "GCG": 44.4, 
    "GAT": 6.7, "GAC": 41.7, "GAA": 2.8, "GAG": 53.5, 
    "GGT": 9.5, "GGC": 62.0, "GGA": 5.0, "GGG": 9.7
}

codon_dic = {"GCT":"A","GCC":"A","GCA":"A","GCG":"A","TGT":"C","TGC":"C","GAA":"E","GAG":"E","GAT":"D","GAC":"D","GGT":"G","GGC":"G","GGA":"G","GGG":"G","TTT":"F","TTC":"F","ATT":"I","ATC":"I","ATA":"I","CAT":"H","CAC":"H","AAA":"K","AAG":"K","ATG":"M","CTT":"L","CTC":"L","CTA":"L","CTG":"L","TTA":"L","TTG":"L","AAT":"N","AAC":"N","CAA":"Q","CAG":"Q","CCT":"P","CCC":"P","CCA":"P","CCG":"P","TCT":"S","TCC":"S","TCA":"S","TCG":"S","AGT":"S","AGC":"S","CGT":"R","CGC":"R","CGA":"R","CGG":"R","AGA":"R","AGG":"R","ACT":"T","ACC":"T","ACA":"T","ACG":"T","TGG":"W","GTT":"V","GTC":"V","GTA":"V","GTG":"V","TAT":"Y","TAC":"Y"}

aa_dic = { "I" : ["ATT", "ATC", "ATA"], "L" : ["CTT", "CTC", "CTA", "CTG", "TTA", "TTG"], "V" : ["GTT", "GTC", "GTA", "GTG"], "F" : ["TTT", "TTC"], "M" : ["ATG"], "C" : ["TGT", "TGC"], "A" : ["GCT", "GCC", "GCA", "GCG"], "G" : ["GGT", "GGC", "GGA", "GGG"], "P" : ["CCT", "CCC", "CCA", "CCG"], "T" : ["ACT", "ACC", "ACA", "ACG"], "S" : ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"], "Y" : ["TAT", "TAC"], "W" : ["TGG"], "Q" : ["CAA", "CAG"], "N" : ["AAT", "AAC"], "H" : ["CAT", "CAC"], "E" : ["GAA", "GAG"], "D" : ["GAT", "GAC"], "K" : ["AAA", "AAG"], "R" : ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"] }

#print(len(codon_freq_table))

def normalize_codon_fre_table(codon_freq_table):
    new_table = {}
    N = 0.0
    for codon in codon_freq_table:
        N += codon_freq_table[codon]
    for codon in codon_freq_table:
        new_table[codon] = round( codon_freq_table[codon] / N * 1000.0, 1)
    return new_table

## make aa_codon_freq_dic
aa_codon_freq_dic = {}
for aa in aa_dic:
    freqs = []
    freq = 0.0
    for codon in aa_dic[aa]:
        freqs.append([codon_freq_table[codon], codon,])
        freq += codon_freq_table[codon]
    freqs.sort(reverse = True)
    aa_codon_freq_dic[aa] = freqs
    #print(aa, freq, freqs, aa_dic[aa])


import math
import random

#random.seed(24)

def translate(mRNA):
    seq = ""
    for i in range(0, len(mRNA),3):
        codon = mRNA[i:i+3]
        seq += codon_dic[codon]
    return seq

def codon_selection_prob(aa):
    population = []
    weights = []
    for freq, codon in aa_codon_freq_dic[aa]:
        population.append( codon )
        weights.append( freq )

    return random.choices(population, weights)[0]
    
def codon_optimization_score(mRNA):
    score = 0.0
    for i in range(0, len(mRNA),3):
        codon = mRNA[i:i+3]
        freq = codon_freq_table[codon]
        score += math.log10(freq)
    return score


def optimize_prob(seq):
    mRNA = ""
    for aa in seq:
        mRNA += codon_selection_prob(aa)
    return mRNA

def optimize_max(seq):
    mRNA = ""
    for aa in seq:
        mRNA += aa_codon_freq_dic[aa][0][1]
    return mRNA

def TestCode():
    print( codon_optimization_score("ATGGGAGGG") )
    print( codon_optimization_score("GGAATGGGG") )
    print( codon_optimization_score("ATGGGCGGC") )
    opt_mRNA = optimize_max("MGG")
    print(opt_mRNA)
    print( codon_optimization_score(opt_mRNA))

    prob_mRNA = optimize_prob("MGG")
    print(prob_mRNA)
    print( codon_optimization_score(prob_mRNA))

def codon_statistics(mRNA):
    major_codon = 0
    minor_codon = 0
    for i in range(0, len(mRNA),3):
        codon = mRNA[i:i+3]
        aa = codon_dic[codon]
        freqs = aa_codon_freq_dic[aa]
        if codon == freqs[0][1]:
            major_codon += 1
        else:
            minor_codon += 1
    return major_codon, minor_codon

import sys

def gc_simulation(codon_freq_table, N = 100000):
    population = []
    weights = []
    mRNA = ""
    for codon in codon_freq_table:
        population.append( codon )
        weights.append( codon_freq_table[codon] )

    for i in range(N):    
        mRNA += random.choices( population, weights )[0] 

    GC = 0.0
    for nt in mRNA:
        if nt in "GC":
            GC += 1
    return GC / len(mRNA) * 100.0


def init_gc_codons():
    HighGC_Codons = []
    LowGC_Codons = []
    for codon in codon_freq_table:
        GC = 0
        for nt in codon:
            if nt in "GC":
                GC += 1
        if GC >= 2:
            HighGC_Codons.append( codon )
        else:
            LowGC_Codons.append( codon )
    return HighGC_Codons, LowGC_Codons
    

def gc_content(mRNA):
    GC = 0.0
    for nt in mRNA:
        if nt in "GC":
            GC += 1
    GC_Content = GC / float( len(mRNA) )
    return GC_Content, GC, float( len(mRNA) )

def gc_adjustment(TargetGC):
    HighGC_Codons, LowGC_Codons = init_gc_codons()

    codon_usage_dic = {}
    Pre_N = 10000
    Post_N = 1000000
    population = []
    weights = []
    mRNA = ""
    for codon in codon_freq_table:
        population.append( codon )
        weights.append( codon_freq_table[codon] )

    for i in range(Pre_N):    
        codon = random.choices( population, weights )[0] 
        codon_usage_dic[codon] = codon_usage_dic.get(codon, 0) + 1
        mRNA += codon

    GC_Content, GC, total_len = gc_content( mRNA)

    #print(HighGC_Codons, len(HighGC_Codons))
    #print(LowGC_Codons, len(LowGC_Codons))
    #print(GC, total_len, TargetGC, GC_Content)


    for i in range(Post_N):
        #if i % 1000 == 0:
        #    print(GC, total_len, TargetGC, GC_Content)

        codon = random.choices( population, weights )[0] 
        #print(codon)
        if TargetGC >= GC_Content:
            if codon in LowGC_Codons: 
                #print("1")
                continue
            #print("2")
            codon_usage_dic[codon] = codon_usage_dic.get(codon, 0) + 1
            pass
        if TargetGC < GC_Content:
            if codon in HighGC_Codons: continue
            codon_usage_dic[codon] = codon_usage_dic.get(codon, 0) + 1
            pass

        for nt in codon:
            if nt in "GC":
                GC += 1
        total_len += 3
        GC_Content = GC / total_len

    return codon_usage_dic

def read_codon_table(filename):
    f = open(filename)
    codon_freq_table = eval(f.read())
    f.close()
    return codon_freq_table


## Load table from file
if len(sys.argv) > 3:
    codon_freq_table = read_codon_table(sys.argv[3])
    print("Load codon frequency table:", codon_freq_table)
    codon_freq_table = normalize_codon_fre_table(codon_freq_table)

if len(sys.argv) < 3:
    mRNA = "ATGGGAGGG"
    seq = "MGG"
    print( "options: score, max, prob, idt, twist, gc ,compare" )
    print( "example 1: python", sys.argv[0], "score", mRNA )
    print( "example 2: python", sys.argv[0], "max", seq )
    print( "example 3: python", sys.argv[0], "prob", seq )
    print( "example 4: python", sys.argv[0], "idt", seq )
    print( "example 5: python", sys.argv[0], "twist", seq )
    print( "example 6: python", sys.argv[0], "gc", "simulation" )
    print( "example 7: python", sys.argv[0], "gc", 0.6 )
    print( "example 8: python", sys.argv[0], "compare", seq )
    
else:
    if sys.argv[1] == "score":
        mRNA = sys.argv[2].upper()
        seq = translate(mRNA)
        print( "Translated seq:", seq )
        print( "Codon optimization score =", codon_optimization_score(mRNA) )
        print( "Number of major and minor codons =", codon_statistics(mRNA) )       
    if sys.argv[1] == "max":
        seq = sys.argv[2].upper()
        opt_mRNA = optimize_max(seq)
        print( opt_mRNA )
        print( "Codon optimization score =", codon_optimization_score(opt_mRNA) )
        print( "Number of major and minor codons =", codon_statistics(opt_mRNA) )    
        assert(seq == translate(opt_mRNA))   
    if sys.argv[1] == "prob":
        seq = sys.argv[2].upper()
        prob_mRNA = optimize_prob(seq)
        print( prob_mRNA )
        print( "Codon optimization score =", codon_optimization_score(prob_mRNA) )
        print( "Number of major and minor codons =", codon_statistics(prob_mRNA) )       
        assert(seq == translate(prob_mRNA))
    if sys.argv[1] == "idt":
        seq = sys.argv[2].upper()
        prob_mRNA = optimize_prob(seq)
        idt_mRNA = prob_mRNA.replace("GCCGCCGCCGCC", "GCCGCTGCCGCC")
        print( idt_mRNA )
        print( "Codon optimization score =", codon_optimization_score(idt_mRNA) )
        print( "Number of major and minor codons =", codon_statistics(idt_mRNA) )       
        assert(seq == translate(idt_mRNA))
    if sys.argv[1] == "twist":
        seq = sys.argv[2].upper()
        prob_mRNA = optimize_prob(seq)
        twist_mRNA = prob_mRNA.replace("GCCGCCGCCGCC", "GCCGCTGCCGCC")
        print( twist_mRNA )
        print( "Codon optimization score =", codon_optimization_score(twist_mRNA) )
        print( "Number of major and minor codons =", codon_statistics(twist_mRNA) )       
        assert(seq == translate(twist_mRNA))
    if sys.argv[1] == "gc":
        if sys.argv[2] == "simulation":
            print( "GC content simulation:", gc_simulation() )
        else:
            target_gc = float(sys.argv[2])
            adjusted_codon_table = gc_adjustment(target_gc)
            print( "Adjusted codon frequency table:", adjusted_codon_table )
            print( "GC content simulation:", gc_simulation(adjusted_codon_table) )
    if sys.argv[1] == "compare":
        codon_freq_table1 = normalize_codon_fre_table(read_codon_table( sys.argv[2] ))

        for codon in codon_freq_table:
            print( codon, "\t", codon_freq_table1[codon], "\t", codon_freq_table[codon] )
        
'''
I 35.7 [[26.6, 'ATC'], [8.0, 'ATT'], [1.1, 'ATA']] ['ATT', 'ATC', 'ATA']
L 89.8 [[65.2, 'CTG'], [13.0, 'CTC'], [4.4, 'CTT'], [4.0, 'TTG'], [2.6, 'CTA'], [0.6, 'TTA']] ['CTT', 'CTC', 'CTA', 'CTG', 'TTA', 'TTG']
V 69.0 [[46.5, 'GTG'], [15.4, 'GTC'], [5.1, 'GTT'], [2.0, 'GTA']] ['GTT', 'GTC', 'GTA', 'GTG']
F 32.1 [[27.1, 'TTC'], [5.0, 'TTT']] ['TTT', 'TTC']
M 25.7 [[25.7, 'ATG']] ['ATG']
C 14.5 [[13.1, 'TGC'], [1.4, 'TGT']] ['TGT', 'TGC']
A 126.29999999999998 [[54.6, 'GCC'], [44.4, 'GCG'], [16.7, 'GCT'], [10.6, 'GCA']] ['GCT', 'GCC', 'GCA', 'GCG']
G 86.2 [[62.0, 'GGC'], [9.7, 'GGG'], [9.5, 'GGT'], [5.0, 'GGA']] ['GGT', 'GGC', 'GGA', 'GGG']
P 63.400000000000006 [[29.5, 'CCC'], [20.7, 'CCG'], [8.1, 'CCT'], [5.1, 'CCA']] ['CCT', 'CCC', 'CCA', 'CCG']
T 52.9 [[27.7, 'ACC'], [15.9, 'ACG'], [5.2, 'ACT'], [4.1, 'ACA']] ['ACT', 'ACC', 'ACA', 'ACG']
S 65.5 [[22.8, 'AGC'], [16.1, 'TCG'], [16.1, 'TCC'], [4.7, 'TCT'], [3.2, 'TCA'], [2.6, 'AGT']] ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC']
Y 25.400000000000002 [[22.8, 'TAC'], [2.6, 'TAT']] ['TAT', 'TAC']
W 13.2 [[13.2, 'TGG']] ['TGG']
Q 40.5 [[36.3, 'CAG'], [4.2, 'CAA']] ['CAA', 'CAG']
N 31.3 [[28.5, 'AAC'], [2.8, 'AAT']] ['AAT', 'AAC']
H 19.4 [[17.2, 'CAC'], [2.2, 'CAT']] ['CAT', 'CAC']
E 56.3 [[53.5, 'GAG'], [2.8, 'GAA']] ['GAA', 'GAG']
D 48.400000000000006 [[41.7, 'GAC'], [6.7, 'GAT']] ['GAT', 'GAC']
K 45.699999999999996 [[43.3, 'AAG'], [2.4, 'AAA']] ['AAA', 'AAG']
R 56.400000000000006 [[34.9, 'CGC'], [11.2, 'CGG'], [4.9, 'CGT'], [2.7, 'AGG'], [2.0, 'CGA'], [0.7, 'AGA']] ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG']
'''