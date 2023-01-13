#!/usr/bin/env python

import sys
from Bio import SeqIO
import gzip
import os


d = './data/'
c = 0
# for each fine in the directory...
for pdbfile in os.listdir(d):
    c+=1
    print('-------\n', c, '\n-------')
    print(pdbfile)
    # open the pdb file and extract the sequence
    with open(d+pdbfile, 'rt') as f:
        seq = next(SeqIO.parse(f, 'pdb-atom'))
    # write the sequence in the directory ./fastas/ in fasta format
    with open('./fastas/'+pdbfile[:-4]+'.fa', 'w') as outf:
        SeqIO.write(seq, outf, 'fasta')
    print('./fastas/'+pdbfile[:-4]+'.fa')
