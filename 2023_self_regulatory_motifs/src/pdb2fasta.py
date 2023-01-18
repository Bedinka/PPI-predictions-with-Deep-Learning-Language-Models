#!/usr/bin/env python

from Bio import SeqIO
import os


d = './data/'
c = 0
l = len(os.listdir(d))
outf = open('./multifasta.fa', 'w')
# for each file in the directory...
for pdbfile in os.listdir(d):
    c+=1
    print('-------\n', c, '/', l, '\n-------')
    # open the pdb file and extract the sequence
    with open(d+pdbfile, 'rt') as f:
        seq = next(SeqIO.parse(f, 'pdb-atom'))
    # Add header to the fastas
    seq.id = pdbfile[:-4]
    seq.description = pdbfile[:-4]
    # write the sequence in the file in fasta format
    SeqIO.write(seq, outf, 'fasta')

