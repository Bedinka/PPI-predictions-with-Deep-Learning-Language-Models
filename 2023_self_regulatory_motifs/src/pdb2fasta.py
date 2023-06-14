#!/usr/bin/env python3

from Bio import SeqIO
import os
from inspect import getmembers


d = './data/'
c = 0
l = len(os.listdir(d))
outf = open('multifasta.fa', 'w')
# for each file in the directory...
for pdbfile in os.listdir(d):
    c+=1
    print('-------\n', c, '/', l, '\n-------')
    # open the pdb file and extract the sequence
    print(pdbfile, '\n-------')
    with open(os.path.join(d, pdbfile), 'rt') as f:
        for seq in (SeqIO.parse(f, 'pdb-seqres')):
            # Add header to the fastas
            if seq.id[:4] == 'XXXX':
                seq.id = pdbfile[:-4]
            # write the sequence in the file in fasta format
            SeqIO.write(seq, outf, 'fasta')
outf.close()