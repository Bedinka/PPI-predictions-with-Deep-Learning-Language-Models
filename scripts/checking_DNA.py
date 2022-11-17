#!/usr/bin/env python
import os 

"""
This script will screen for PDB structures with DNA molecules 
"""

import os

"""
This script will remove protein containing DNA from selected list
"""

pdb_DNA=[]
DNA_positive=[]

    

#Check for DNA strands in selected PDB IDs
for pdb in pdb_id:
    f=open(f"pdb_new_files/{pdb}/{pdb}.pdb","r")
    for lines in f.readlines():
        if "MOLECULE: DNA (5'-" in lines[:]:
            DNA_positive.append(pdb)
        else:
            continue


for pdbs in pdb_DNA:
    pdb_id.remove(pdbs)
    
print("Checking for protein contiaing DNA strans: ")
print(any(items in pdb_DNA for items in pdb_id))


