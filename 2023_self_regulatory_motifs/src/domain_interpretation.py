#!/usr/bin/env python

import sys
from Bio import SeqIO
from Bio import SearchIO
import pyhmmer

file = sys.argv[1]

with open(file, 'r') as f:
    for qresult in SearchIO.parse(f, 'hmmscan3-domtab'):
        print(qresult.id)