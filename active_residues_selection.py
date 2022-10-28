import pandas as pd
import os 
import glob

"""
This script select Active residues based on their FREESASA and ScanNet values 
"""

#Retrieve PDB IDs and chain IDs 
pdbs = []
pdb_dict = {}

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
        
        
no_scannet = []

for key,value in pdb_dict.items():
    pdb=key
    a=value[0]
    id_a = value[0].split("-")[1].split(".")[0]
    b=value[1]
    id_b = value[1].split("-")[1].split(".")[0]

    with open(f"{pdb}/active_residues_{a}.list","w") as f:
        #Load SASA and ScanNet datasets
        dfa=pd.read_csv(f"{pdb}/SASA_{id_a}.csv",sep=",")
        dfa_scan=pd.read_csv(f"{pdb}/{a}_single_ScanNet_interface_noMSA/predictions_{a}.csv",sep=",")
        dfa_scan["Residue Index"] = dfa_scan["Residue Index"].astype("int")
        
        if dfa_scan["Binding site probability"].mean() <= 0.3

        #Change 'Residue number' column name in dfa to 'Reside index'
        dfa.rename(columns={"Residue number":"Residue Index"},inplace=True)
        dfa["Residue Index"]=dfa["Residue Index"].astype("int")

        #Merge dataframes
        comb_a = pd.merge(dfa, dfa_scan, on="Residue Index")

        #Extracting active residue based on PBscore >=0.3 and relative SASA >= 0.4 '40%'
        comb_a = comb_a[(comb_a["Binding site probability"] >= 0.3) &
                     (comb_a["M.relative SASA"] >= 0.4)]

        ls_a=" ".join(str(num) for num in comb_a["Residue Index"].to_list())
        f.write(f"{ls_a}")

        comb_a.to_csv(f"{pdb}/SASA_SCANNET_{id_a}.csv",sep=",",header=True,index=False)


    with open(f"{pdb}/active_residues_{b}.list","w") as f:
        dfb=pd.read_csv(f"{pdb}/SASA_{id_b}.csv",sep=",")
        dfb_scan=pd.read_csv(f"{pdb}/{b}_single_ScanNet_interface_noMSA/predictions_{b}.csv",sep=",")
        dfb_scan["Residue Index"] = dfb_scan["Residue Index"].astype("int")

        #Change 'Residue number' column name in dfb to 'Reside index'
        dfb.rename(columns={"Residue number":"Residue Index"},inplace=True)
        dfb["Residue Index"]=dfb["Residue Index"].astype("int")

        #Merge dataframes
        comb_b = pd.merge(dfb, dfb_scan, on="Residue Index")

        #Extracting active residue based on PBscore >=0.3 and M.relative SASA >= 0.4 '40%'
        comb_b = comb_b[(comb_b["Binding site probability"] >= 0.3) &
                     (comb_b["M.relative SASA"] >= 0.4)]

        ls_b=" ".join(str(num) for num in comb_b["Residue Index"].to_list())
        f.write(f"{ls_b}")

        comb_b.to_csv(f"{pdb}/SASA_SCANNET_{id_b}.csv",sep=",",header=True,index=False)

