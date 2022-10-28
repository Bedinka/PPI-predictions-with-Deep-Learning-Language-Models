#!/usr/bin/env python
import os
import glob 

"""
generate run1.cns file 
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

with open("create_run1.csh","w") as f:
    f.write("#!/bin/csh\n")
    for key in pdb_dict.keys():
        f.write(f"cd {key}/\n")
        f.write(f"source ../../../../haddock_configure.csh\n")
        f.write("haddock2.4\n")
        f.write("cd ../\n")
        f.write("\n")
