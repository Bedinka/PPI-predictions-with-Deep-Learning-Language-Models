import glob
import os 
"""
combine the active and passive residue list
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

with open("docking_new/combination.sh","w") as f:
    f.write("#!/bin/bash\n")
    for key,value in selected_pdbs.items():
        pdb=key
        a=value[0]
        b=value[1]
        f.write(f"""cat {pdb}/active_residues_{a}.list empty.txt {pdb}/passive_residues_{a}.list > {pdb}/active_passive_{a}.list\n""")
        f.write(f"""cat {pdb}/active_residues_{b}.list empty.txt {pdb}/passive_residues_{b}.list > {pdb}/active_passive_{b}.list\n""")
        
        
        
        
