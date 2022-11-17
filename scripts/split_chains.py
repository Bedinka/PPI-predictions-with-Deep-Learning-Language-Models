#!/usr/bin/env python
import os

"""
Spliting pdb files by chain
"""

def splitPDBbyChain(filename):
    chains = {}
    
    f = open(f"pdb_new_files/{filename}/{filename}_clean.pdb")
    for line in f.readlines():  
        if line[:4] != "ATOM":
            continue
        chainID = line[21] 
        if chainID not in chains:
            chains[chainID] = [] 
        chains[chainID].append(line)
    print(chains.keys())

    [dirname, basename] = os.path.split(f"pdb_new_files/{filename}/")
    for chainID in chains: 
        chain_filename = f"pdb_new_files/{filename}/{filename}"+ basename[:-4] + "-" + chainID +".pdb"
        print(chain_filename)
        fout = open( chain_filename, "w" )
        for line in chains[chainID]:
            fout.write(line)
        fout.close()
    f.close()

for filename in pdb_id[:]:
    splitPDBbyChain(filename)
    print("\n")


