#!/usr/bin/env python 

"""
This script will prepare docked files for RMSD analysis 
"""


import pandas as pd
import os 

docked=[]
pdbs={}


df=pd.read_csv("selected_chains_for_docking.csv",sep=",")

#Selecting only docked models 
for pdb in os.listdir("docked_05/"):
	if ".txt" in pdb or ".sh" in pdb or ".csh" in pdb:
		continue
	else:
		if f"{}_run1_fcc_matrix.job{}" in os.listdir(f"docking_05/{pdb}/run1/"):
			docked.append(pdb)
		else:
			continue

#create dictionary with selected PDBIDs and their corresponding chains 

for items, rows in df.iterrows():
	if rows["PDB ID"] in docked:
		pdbs[ rows["PDB ID"]  ] = [ rows["PDB ID 1"]  , rows["PDB ID 2"]  ]
	else:
		continue  


#creating complex files for RMSD analysis
for key,value in pdbs.items():
    try:
        os.system(f"cat ../docking_05/{key}/{value[0]}_docking_final.pdb ../docking_05/{key}/{value[1]}_docking_final.pdb > ../docking_05/{key}/complex.pdb")
    except:
        print("An error has occured with {key}")


#modifying the ana_structures.csh file to perform the RMSD calculation on complex structure
try:
    os.system("sed 's#set refe=`head -n1 file.nam`#set refe=../../../complex.pdb#' ana_structures.csh > ana_struct_01.csh")
    os.system("sed 's#set refe=`head -n1 file.nam`#set refe=../../../../complex.pdb#' ana_structures.csh > ana_struct_water.csh")
except:
    raise Exception("An error occurred during the sed process")


#copy modified ana_structures to structure folder
for pdb in pdbs.keys():
    os.system(f"cp ana_struct_01.csh ../docking_05/{pdb}/run1/structures/it0")
    os.system(f"cp ana_struct_01.csh ../docking_05/{pdb}/run1/structures/it1")
    os.system(f"cp ana_struct_water.csh ../docking_05/{pdb}/run1/structures/it1/water")
