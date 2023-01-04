import math
import os
import pandas as pd 
import numpy as np 
from sklearn.metrics import confusion_matrix, matthews_corrcoef

"""
    This script will compute the Matthews correlation coefficient for
    interface and non-interface residues predicted by various predictors.
"""

predictors = [x.lstrip() for x in "scannet, scriber, csm_potential, masif".split(",")]

pdbs = [x for x in os.listdir() if x.endswith(".py") == False and x.endswith(".csv") == False and x.endswith(".sh") == False 
and x.endswith(".png") == False and len(os.listdir(x)) > 11]

data = {}

def open_dfs(pdb_id:str , chain_id: str):
    global df 
    df = pd.read_csv(f"{pdb_id}/Descriptors_{chain_id}.csv",sep=",")
    df.dropna(inplace=True)


def mcc_prep():
    """
    Select the descriptor columns  
        
            Parameters: 
                    pdb_id(str): PDB ID
                    chain_id(str): Chain ID
                    
            Returns:
                    descriptor lists 
    
    
    """
    global interface, scannet_prediction, scriber_prediction,  csm_potential_prediction, masif_prediction

    interface = np.array(df["Label"].to_list())
    scannet_prediction = np.array(df["Scannet"].to_list())
    scriber_prediction = np.array(df["Scriber"].to_list())
    csm_potential_prediction = np.array(df["csm_potential"].to_list())
    masif_prediction = np.array(df["masif"].to_list())
    
def mcc(pdb_id, pred):
    global cal, matrix
    TN, TP, FP, FN = 0, 0, 0, 0
    cal = 0
    matrix = confusion_matrix(interface,globals()[f"{pred}_prediction"])
    print(f"{pdb_id} matrix shape: ",matrix.shape)
    if matrix.shape[0] == 1:
        return (f"check {pdb}")
    else:
        TP = matrix[0][0]
        FP = matrix[0][1]
        FN = matrix[1][0]
        TN = matrix[1][1]

        cal = round( ( (TP*TN) - (FP*FN) ) / math.sqrt( (TP+FP) * (TP+FN) * (TN+FP) * (TN+FN) ),3)
    print(matrix, pred)
for pdb in pdbs: 
    desc = [x.split("/")[1] for x in glob.glob(f"{pdb}/Descriptor*")]
    for i in desc:
        chain_id = i.split(".")[0].split("_")[1]
        unique_id = pdb+"_"+chain_id
        open_dfs(pdb,chain_id)
        mcc_prep()
        mcc_values = []
        for a in predictors:
            mcc(pdb,a)
            mcc_values.append(cal)
        data[unique_id] = mcc_values
data
df_dic = pd.DataFrame.from_dict(data).transpose()
df_dic.reset_index(inplace=True)
df_dic.to_csv("test.csv",sep=",",index=None)
df_dic.to_csv("test.csv",sep=",",index=None)
df_dic.columns = "PDB_ID, scannet_MCC, scriber_MCC, csm_potential_MCC, masif_MCC".split(",")
df_dic.to_csv("output_mcc.csv",sep=",", index=False) 
