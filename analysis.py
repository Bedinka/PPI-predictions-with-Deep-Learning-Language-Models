#!/usr/bin/env python 
import os
import glob
import pandas as pd

pdbs = []
pdb_dict = {}

for i in os.listdir():
    if "." in i:
	continue
    else:
	pdbs.append(i)

for i in pdbs:
    pdb_dict[i] = ["{}/it0/structures_rmsd-sorted.dat".format(i), "{}/it1/structures_rmsd-sorted.dat".format(i), "{}/it1/water/structures_rmsd-sorted.dat".format(i)]

for key, value in pdb_dict.items():
    for i in range(len(value)):
        if "/it0/" in value[i]:
            df = pd.read_csv("{}".format(value),sep=" ")
            for items, rows in df.iterrows():
                it0_id.append(rows["#struc"])
                it0.append(rows["rmsd_all"])
        elif "/it1/" in value[i]:
            df1 = pd.read_csv("{}".format(value),sep=" ")
            for items, rows in df1.iterrows():
                it1_id.append(rows["#struc"])
                it1.append(rows["rmsd_all"]) 

for keys in pdb_dict.keys(): 
    dfw = pd.read_csv("{}/it1/water/structures_rmsd-sorted.dat".format(keys),sep=" ")
    for items, rows in dfw.iterrows():
        it1_water_id.append(rows["#struc"])
        it1_water.append(rows["rmsd_all"]) 

#Check if lists are empty 
try:
    print("Check it0:", len(it0) != 0)
except:
    print("it0 is empty")

try:
    print("Check it1:", len(it1) != 0)
except:
    print("it1 is empty")
  
try:  
    print("Check it1_water:", len(it1_water) !=0)
except:
    print("it1_water is empty")


#Compiling and saving data to csv files
compilation_it0 = {"it0_id" : it0_id,
        "it0" : it0}

compilation_it1 = {"it1_id" : it1_id,
        "it1" : it1}
        
compilation_it1_water = {"it1_water_id" : it1_water_id,
        "it1_water" : it1_water}

df_it0 = pd.DataFrame.from_dict(compilation_it0)
df_it1 = pd.DataFrame.from_dict(compilation_it1)
df_it1_water = pd.DataFrame.from_dict(compilation_it1_water)

df_it0.to_csv("it0.csv",sep=",",header=True,index=False)
df_it1.to_csv("it1.csv",sep=",",header=True,index=False)
df_it1_water.to_csv("it1_water.csv",sep=",",header=True,index=False)
