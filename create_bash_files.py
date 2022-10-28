#!/usr/bin/env python 
import os 
import glob 

"""
Create bash files for docking  
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


start=0
end=10


while end <= 200:
    with open(f"docking_final_{start}_{end}.csh","w") as f:
        f.write("#!/bin/csh\n")
        for keys in list(pdb_dict.keys())[start:end]:
                f.write(f"source ../../../haddock_configure.csh\n")
                f.write(f"cd {keys}/run1/\n")
                f.write("haddock2.4\n")
                f.write(f"echo '{keys} has finished' \n")
                f.write("cd ../../ \n")
    start=end
    end+=10
