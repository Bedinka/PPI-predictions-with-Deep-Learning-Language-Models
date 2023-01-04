import numpy as np
import pandas as pd
import os
import sys
from Bio.PDB import *

pdb = []
dssp_type = ["H","B","E","G","I","T","S","L"]
data = {}

#Get PDB IDs from MCC dataframe 
with open("output_mcc.csv","r") as f:
    next(f)
    for lines in f.readlines():
        pdb.append(lines[:6])

#Open PDB file and get its sequence length 
def open_pdb(pdb_id,chain_id):
    global structure, seq_len
    parser = PDBParser()
    structure = parser.get_structure(f"{pdb_id}",f"{pdb_id}/{chain_id}.pdb")
    residues = []
    for res in structure.get_residues():
        if "W" in res.get_full_id()[3][0] or "H_" in res.get_full_id()[3][0]:
            continue
        else:
            residues.append(res.get_full_id()[3][1])
    seq_len = len(residues)

#Calculate structural properties of the structure
def dssp_analysis(pdb_id,chain_id):
    global H_count, B_count, E_count, G_count, I_count, T_count, S_count, L_count, dssp_profile
    dssp_profile = []
    model = structure[0]
    dssp_structure = DSSP(model,"{}".format(f"{pdb_id}/{chain_id}.pdb"))
    for values in dssp_structure.property_list:
        dssp_profile.append(values[2])
    H_count = round((len([str(x) for x in dssp_profile if x == "H"]) / seq_len) * 100 , 2)
    B_count = round((len([str(x) for x in dssp_profile if x == "B"]) / seq_len) * 100 , 2)
    E_count = round((len([str(x) for x in dssp_profile if x == "E"]) / seq_len) * 100 , 2)
    G_count = round((len([str(x) for x in dssp_profile if x == "G"]) / seq_len) * 100 , 2)
    I_count = round((len([str(x) for x in dssp_profile if x == "I"]) / seq_len) * 100 , 2)
    T_count = round((len([str(x) for x in dssp_profile if x == "T"]) / seq_len) * 100 , 2)
    S_count = round((len([str(x) for x in dssp_profile if x == "S"]) / seq_len) * 100 , 2)
    L_count = round((len([str(x) for x in dssp_profile if x == "-"]) / seq_len) * 100 , 2)
    data[chain_id] = [H_count, B_count, E_count, G_count, I_count, T_count, S_count, L_count]
    
#Extract DSSP information

for chain_id in pdb:
    pdb_id_ = chain_id.split("_")[0]
    open_pdb(pdb_id_, chain_id)
    dssp_analysis(pdb_id_,chain_id)
    
# #Save data into CSV file 
df = pd.DataFrame.from_dict(data).transpose()
df.reset_index(inplace=True)
df.columns = "PDB_ID, H_count, B_count, E_count, G_count, I_count, T_count, S_count, L_count".split(", ")
df.to_csv("DSSP_data.csv",index=None, sep=",")
