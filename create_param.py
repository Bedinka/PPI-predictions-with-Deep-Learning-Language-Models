#!/usr/bin/env python 
import os 
import glob 

"""
create run.param file 
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

for key,value in pdb_dict.items():
    pdb = key
    a = value[0]
    id_a = value[0].split("-")[1]
    b = value[1]
    id_b = value[1].split("-")[1]
    
    with open(f"{pdb}/run.param","w") as f:
        f.write(f"AMBIG_TBL=active_passive_ambig.tbl \nHADDOCK_DIR=/storage/scratch01/users/mshehata/haddock2.4-2022-01/ \nN_COMP=2 \nPDB_FILE1={a}.pdb \nPDB_FILE2={b}.pdb \nPROJECT_DIR=./ \nPROT_SEGID_1={id_a} \nPROT_SEGID_2={id_b} \nRUN_NUMBER=1")
        
    
