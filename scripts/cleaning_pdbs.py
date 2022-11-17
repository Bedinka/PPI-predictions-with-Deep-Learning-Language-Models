#!/usr/bin/env python
import os

"""
Clean structures from HETATM, ANISOU, and UNK enteries
"""

pdb_noanisou=[]
pdb_nounk=[]

for pdbs in pdb_id:
    with open(f"pdb_new_files/{pdbs}/{pdbs}.pdb","r") as f:
        if "ANISOU" in f.read():
            os.system(f"grep -v 'ANISOU' pdb_new_files/{pdbs}/{pdbs}.pdb > pdb_new_files/{pdbs}/{pdbs}_noanisou.pdb")
            os.system(f"python pdb_cleaner.py pdb_new_files/{pdbs}/{pdbs}_noanisou.pdb > pdb_new_files/{pdbs}/{pdbs}_clean.pdb")

        else:
            os.system(f"python pdb_cleaner.py pdb_new_files/{pdbs}/{pdbs}.pdb > pdb_new_files/{pdbs}/{pdbs}_clean.pdb")


