#!/usr/bin/env 
import os 

"""
This script will model missing atoms and side chains in PDB IDs that have been filtered.
Open babel and Ambertools should be installed.
"""


#Model missing atoms and side chains 
with open("docking/model_missing_side_chains.sh","w") as f:
    f.write("#!/bin/bash\n")
    for key,value in pdb_dict.items():
        with open("docking/{}/{}-{}.pdb".format(key,key,str(value[0]))) as pdb1:
            for lines in pdb1.readlines(1):
                a=lines[23:26]
        pdb1.close()      
        with open("docking/{}/{}-{}.pdb".format(key,key,str(value[1]))) as pdb2:
            for lines in pdb2.readlines(1):
                b=lines[23:26]
        pdb2.close()

        os.mkdir(f"docking/{key}/test/")
        
        f.write("pdb4amber -i docking/{}/{}-{}.pdb -o docking/{}/test/{}.pdb --add-missing-atoms\n".format(key,key,str(value[0]),key,str(value[0]))) 


        f.write("./renumber_pdb.py -i docking/{}/test/{}.pdb -s {} -r > docking/{}/test/{}-{}.pdb\n".format(key,str(value[0]),a,key,key,str(value[0])))

        f.write("obabel -ipdb docking/{}/test/{}-{}.pdb -opdb > docking/{}/{}-{}_noH.pdb\n".format(key,key,str(value[0]),key,key,str(value[0]))) 



        f.write("pdb4amber -i docking/{}/{}-{}.pdb -o docking/{}/test/{}.pdb --add-missing-atoms\n".format(key,key,str(value[1]),key,str(value[1]))) 



        f.write("./renumber_pdb.py -i docking/{}/test/{}.pdb -s {} -r > docking/{}/test/{}-{}.pdb\n".format(key,str(value[1]),b,key,key,str(value[1])))

        f.write("obabel -ipdb docking/{}/test/{}-{}.pdb -opdb > docking/{}/{}-{}_noH.pdb\n".format(key,key,str(value[1]),key,key,str(value[1]))) 

        
print("Single chain structures were removed: √")
