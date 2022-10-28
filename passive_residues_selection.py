import os
import sys 
import pandas as pd
import glob

"""
This script will generate passive residues from selected active residues 
"""

def open_active_residues(x,y):
    global active_residues
    with open(f"{x}/active_residues_{y}.list","r") as p:
        active_residues = p.read()
    active_residues = ",".join(x for x in active_residues.split(" ") if "-" not in x)
    

#Extract PDB and chains IDs
pdbs = []
pdb_dict = {}

def get_pdbs():
    for i in os.listdir():
        if "." in i:
            continue
        else:
            pdbs.append(i)

        for i in pdbs:
            pdb_list = []
            sel = glob.glob(f"{i}/*.pdb")
            for a in sel:
                if "complex" in a:
                    continue
                else:
                    pdb_list.append( a.split(".")[0].split("/")[1] )
                    pdb_dict[i] = pdb_list

get_pdbs()

#Write passive_residue.sh file

with open("passive_residues.sh","w") as f:
    f.write("#!/bin/bash\n")
    for key, value in pdb_dict.items():
        pdb=key
        a=value[0]
        b=value[1]
        
        open_active_residues(pdb,a)
        f.write(f"./passive_from_active.py {pdb}/{a}.pdb {active_residues} > {pdb}/passive_residues_{a}.list\n")
        
        open_active_residues(pdb,b)
        f.write(f"./passive_from_active.py {pdb}/{b}.pdb {active_residues} > {pdb}/passive_residues_{b}.list\n")

    
