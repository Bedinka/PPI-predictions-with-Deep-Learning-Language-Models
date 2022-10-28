#!/usr/bin/env python
import os
import glob 

"""
Edit run.cns file 
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

with open("edit_run_cns.csh","w") as f:
    f.write("#!/bin/csh\n")
    for key in selected_pdbs:
        pdb = key
        a = value[0]
        b = value[1]
        
        f.write(f"cd {key}/run1/\n")
        f.write(f"""sed 's#/storage/scratch01/users/mshehata/haddock2.4-2022-01/cns/cns_solve_1.31-UU/intel-x86_64bit-linux/bin/cns#/storage/scratch01/users/mshehata/cns_solve_1.3/intel-x86_64bit-linux/bin/cns#; s/qsub/csh/g; s/structures_0=1000/structures_0=100/g; s/structures_1=200/structures_1=50/g; s/anastruc_1=200/anastruc_1=25/g' run.cns > run1.cns\n""")
        f.write("mv run1.cns run.cns\n")
        f.write("cd ../../\n")
        f.write("\n")
