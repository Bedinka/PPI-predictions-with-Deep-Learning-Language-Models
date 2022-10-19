#!/usr/bin/env python

"""
This script will download PDB IDs from interaction list provided from Struct2Graph output data 
"""

import os 

pdb_id=[]
count=0

f=open("interactions_data_self.txt") 

for lines in f.readlines()[:300]:
    pdb_id.append(lines[:4])


os.mkdir("pdb_new_files/")


for pdbs in pdb_id:
    if os.path.exists(f"pdb_new_files/{pdbs}/"):
        continue
    else:
        os.mkdir(f"pdb_new_files/{pdbs}/")
        os.system(f"python pdb_fetch.py {pdbs} > pdb_new_files/{pdbs}/{pdbs}.pdb")
        count+=1
        print("{} {} downloaded".format(count,pdbs))
