#!/usr/bin/env python
import os 
import glob

"""
generate list tbl file from combined active_passive
"""

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

with open("prep_ambig.sh","w") as f:
    f.write("#!/bin/bash\n")
    for key,value in pdb_dict.items():
            a=value[0]
            b=value[1]
            pdb=key

            f.write(f"./active-passive-to-ambig.py {pdb}/active_passive_{a}.list {pdb}/active_passive_{b}.list > {pdb}/active_passive_ambig.tbl\n")


    
    
